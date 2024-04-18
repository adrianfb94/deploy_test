/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Interpolate remapcon        First order conservative remapping
      Interpolate remapcon2       Second order conservative remapping
      Interpolate remapbil        Bilinear interpolation
      Interpolate remapbic        Bicubic interpolation
      Interpolate remapdis        Distance-weighted averaging
      Interpolate remapnn         Nearest neighbor remapping
      Interpolate remaplaf        Largest area fraction remapping
      Genweights  gencon          Generate first order conservative remap weights
      Genweights  gencon2         Generate second order conservative remap weights
      Genweights  genbil          Generate bilinear interpolation weights
      Genweights  genbic          Generate bicubic interpolation weights
      Genweights  gendis          Generate distance-weighted averaging weights
      Genweights  gennn           Generate nearest neighbor weights
      Genweights  genlaf          Generate largest area fraction weights
      Remap       remap           SCRIP grid remapping
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "remap.h"
#include "grid.h"


enum {REMAPCON, REMAPCON2, REMAPBIL, REMAPBIC, REMAPDIS, REMAPNN, REMAPLAF, REMAPSUM,
      GENCON, GENCON2, GENBIL, GENBIC, GENDIS, GENNN, GENLAF, REMAPXXX, REMAPYCON, GENYCON};

enum {HEAP_SORT, MERGE_SORT};

static
void get_map_type(int operfunc, int *map_type, int *submap_type, int *num_neighbors, int *remap_order)
{
  switch ( operfunc )
    {
    case REMAPYCON:
    case GENYCON:
      *map_type = MAP_TYPE_CONSERV_YAC;
      *remap_order = 1;
      break;
    case REMAPCON:
    case GENCON:
      *map_type = MAP_TYPE_CONSERV;
      *remap_order = 1;
      break;
    case REMAPCON2:
    case GENCON2:
      *map_type = MAP_TYPE_CONSERV;
      *remap_order = 2;
      break;
    case REMAPLAF:
    case GENLAF:
      *map_type = MAP_TYPE_CONSERV_YAC;
      *submap_type = SUBMAP_TYPE_LAF;
      break;
    case REMAPSUM:
      *map_type = MAP_TYPE_CONSERV;
      *submap_type = SUBMAP_TYPE_SUM;
      break;
    case REMAPBIL:
    case GENBIL:
      *map_type = MAP_TYPE_BILINEAR;
      break;
    case REMAPBIC:
    case GENBIC:
      *map_type = MAP_TYPE_BICUBIC;
      break;
    case REMAPDIS:
    case GENDIS:
      *map_type = MAP_TYPE_DISTWGT;
      if ( *num_neighbors == 0 ) *num_neighbors = 4;
      break;
    case REMAPNN:
    case GENNN:
      *map_type = MAP_TYPE_DISTWGT;
      *num_neighbors = 1;
      break;
    default:
      cdoAbort("Unknown mapping method");
      break;
    }
}

static
int maptype2operfunc(int map_type, int submap_type, int num_neighbors, int remap_order)
{
  int operfunc = -1;

  if      ( map_type == MAP_TYPE_CONSERV )     operfunc = (remap_order == 2) ? REMAPCON2 : REMAPCON;
  else if ( map_type == MAP_TYPE_CONSERV_YAC ) operfunc = (submap_type == SUBMAP_TYPE_LAF) ? REMAPLAF : REMAPYCON;
  else if ( map_type == MAP_TYPE_BILINEAR )    operfunc = REMAPBIL;
  else if ( map_type == MAP_TYPE_BICUBIC )     operfunc = REMAPBIC;
  else if ( map_type == MAP_TYPE_DISTWGT )     operfunc = (num_neighbors == 1) ? REMAPNN : REMAPDIS;
  else  cdoAbort("Unsupported mapping method (map_type = %d)", map_type);

  return operfunc;
} 

static
void print_remap_info(int operfunc, int remap_genweights, remapgrid_t *src_grid, remapgrid_t *tgt_grid, int nmiss)
{
  char line[256], tmpstr[256];

  line[0] = 0;

  if      ( operfunc == REMAPBIL  || operfunc == GENBIL  )  strcpy(line, "SCRIP bilinear");
  else if ( operfunc == REMAPBIC  || operfunc == GENBIC  )  strcpy(line, "SCRIP bicubic");
  else if ( operfunc == REMAPNN   || operfunc == GENNN   )  strcpy(line, "Nearest neighbor");
  else if ( operfunc == REMAPDIS  || operfunc == GENDIS  )  strcpy(line, "Distance-weighted average");
  else if ( operfunc == REMAPCON  || operfunc == GENCON  )  strcpy(line, "SCRIP first order conservative");
  else if ( operfunc == REMAPCON2 || operfunc == GENCON2 )  strcpy(line, "SCRIP second order conservative");
  else if ( operfunc == REMAPLAF  || operfunc == GENLAF  )  strcpy(line, "YAC largest area fraction");
  else if ( operfunc == REMAPYCON || operfunc == GENYCON )  strcpy(line, "YAC first order conservative");
  else                                                      strcpy(line, "Unknown");

  if ( remap_genweights )
    strcat(line, " weights from ");
  else
    strcat(line, " remapping from ");

  strcat(line, gridNamePtr(gridInqType(src_grid->gridID)));
  if ( src_grid->rank == 2 )
    snprintf(tmpstr, sizeof(tmpstr), " (%dx%d)", src_grid->dims[0], src_grid->dims[1]);
  else
    snprintf(tmpstr, sizeof(tmpstr), " (%d)", src_grid->dims[0]);
  strcat(line, tmpstr);
  strcat(line, " to ");
  strcat(line, gridNamePtr(gridInqType(tgt_grid->gridID)));
  if ( tgt_grid->rank == 2 )
    snprintf(tmpstr, sizeof(tmpstr), " (%dx%d)", tgt_grid->dims[0], tgt_grid->dims[1]);
  else
    snprintf(tmpstr, sizeof(tmpstr), " (%d)", tgt_grid->dims[0]);
  strcat(line, tmpstr);
  strcat(line, " grid");

  if ( nmiss > 0 )
    {
      snprintf(tmpstr, sizeof(tmpstr), ", with source mask (%d)", gridInqSize(src_grid->gridID)-nmiss);
      strcat(line, tmpstr);
    }

  cdoPrint(line);
}

static
void print_remap_warning(const char *remap_file, int operfunc, remapgrid_t *src_grid, int nmiss)
{
  char line[256];
  char tmpstr[256];

  line[0] = 0;
  (void)operfunc;
  /*
  if      ( operfunc == REMAPBIL  || operfunc == GENBIL  )  strcpy(line, "SCRIP bilinear");
  else if ( operfunc == REMAPBIC  || operfunc == GENBIC  )  strcpy(line, "SCRIP bicubic");
  else if ( operfunc == REMAPNN   || operfunc == GENNN   )  strcpy(line, "SCRIP nearest neighbor");
  else if ( operfunc == REMAPDIS  || operfunc == GENDIS  )  strcpy(line, "SCRIP distance-weighted average");
  else if ( operfunc == REMAPCON  || operfunc == GENCON  )  strcpy(line, "SCRIP first order conservative");
  else if ( operfunc == REMAPCON2 || operfunc == GENCON2 )  strcpy(line, "SCRIP second order conservative");
  else if ( operfunc == REMAPLAF  || operfunc == GENLAF  )  strcpy(line, "YAC largest area fraction");
  else if ( operfunc == REMAPYCON || operfunc == GENYCON )  strcpy(line, "YAC first order conservative");
  else                                                      strcpy(line, "Unknown");

  strcat(line, " remap weights from ");
  */
  strcat(line, "Remap weights from ");
  strcat(line, remap_file);
  strcat(line, " not used, ");
  strcat(line, gridNamePtr(gridInqType(src_grid->gridID)));
  if ( src_grid->rank == 2 )
    snprintf(tmpstr, sizeof(tmpstr), " (%dx%d)", src_grid->dims[0], src_grid->dims[1]);
  else
    snprintf(tmpstr, sizeof(tmpstr), " (%d)", src_grid->dims[0]);
  strcat(line, tmpstr);
  strcat(line, " grid");

  if ( nmiss > 0 )
    {
      snprintf(tmpstr, sizeof(tmpstr), " with mask (%d)", gridInqSize(src_grid->gridID)-nmiss);
      strcat(line, tmpstr);
    }

  strcat(line, " not found!");

  cdoWarning(line);
}


double remap_threshhold = 2;
double gridsearch_radius = 180;
int remap_test = 0;
int remap_order = 1;
int remap_non_global = FALSE;
int remap_num_srch_bins = 180;
static bool lremap_num_srch_bins = false;
static bool remap_extrapolate = false;
static bool lextrapolate = false;
int max_remaps = -1;
int sort_mode = HEAP_SORT;
double remap_frac_min = 0;
int remap_genweights = TRUE;

static
void get_remap_env(void)
{
  char *envstr;

  envstr = getenv("MAX_REMAPS");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival > 0 )
	{
	  max_remaps = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set MAX_REMAPS to %d", max_remaps);
	}
    }

  envstr = getenv("REMAP_MAX_ITER");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_set_int(REMAP_MAX_ITER, ival);
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_MAX_ITER to %d", ival);
	}
    }
  /*
  envstr = getenv("REMAP_ORDER");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_order = ival;
	  if ( remap_order == 0 ) remap_order = 1;
	  if ( remap_order != 1 && remap_order != 2 )
	    cdoAbort("REMAP_ORDER must be 1 or 2");
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_ORDER to %d", remap_order);
	}
    }
  */
  envstr = getenv("REMAP_TEST");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_test = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_TEST to %d", remap_test);
	}
    }

#if defined(_OPENMP)
  sort_mode = (ompNumThreads == 1) ? HEAP_SORT : MERGE_SORT;
#endif

  envstr = getenv("REMAP_SORT_MODE");
  if ( envstr )
    {
      if      ( strcmp(envstr, "heap")  == 0 ) sort_mode = HEAP_SORT;
      else if ( strcmp(envstr, "merge") == 0 ) sort_mode = MERGE_SORT;

      if ( cdoVerbose )
	{
	  if      ( sort_mode == HEAP_SORT )
	    cdoPrint("Set sort_mode to HEAP_SORT");
	  else if ( sort_mode == MERGE_SORT )
	    cdoPrint("Set sort_mode to MERGE_SORT");
	}
    }

  envstr = getenv("REMAP_THRESHHOLD");
  if ( envstr )
    {
      double fval = atof(envstr);
      if ( fval > 0 )
	{
	  remap_threshhold = fval;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_THRESHHOLD to %g", remap_threshhold);
	}
    }

  remap_set_threshhold(remap_threshhold);

  envstr = getenv("CDO_REMAP_RADIUS");
  if ( envstr )
    {
      double fval = radius_str_to_deg(envstr);
      if ( fval < 0 || fval > 180 ) cdoAbort("%s=%g out of bounds (0-180 deg)!", "CDO_REMAP_RADIUS", fval);
      gridsearch_radius = fval;
      if ( cdoVerbose ) cdoPrint("Set CDO_REMAP_RADIUS to %g", gridsearch_radius);
    }

  envstr = getenv("CDO_GRIDSEARCH_RADIUS");
  if ( envstr )
    {
      double fval = radius_str_to_deg(envstr);
      if ( fval < 0 || fval > 180 ) cdoAbort("%s=%g out of bounds (0-180 deg)!", "CDO_GRIDSEARCH_RADIUS", fval);
      gridsearch_radius = fval;
      if ( cdoVerbose ) cdoPrint("Set CDO_GRIDSEARCH_RADIUS to %g", gridsearch_radius);
    }
  
  if ( cdoVerbose )
    cdoPrint("remap_radius = %g deg", gridsearch_radius);

  envstr = getenv("REMAP_AREA_MIN");
  if ( envstr )
    {
      double fval = atof(envstr);
      if ( fval > 0 )
	{
	  remap_frac_min = fval;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_AREA_MIN to %g", remap_frac_min);
	}
    }

  envstr = getenv("REMAP_NUM_SRCH_BINS");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival > 0 )
	{
	  remap_num_srch_bins = ival;
	  lremap_num_srch_bins = true;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_NUM_SRCH_BINS to %d", remap_num_srch_bins);
	}
    }

  envstr = getenv("REMAP_NON_GLOBAL");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival >= 0 )
	{
	  remap_non_global = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_NON_GLOBAL to %d", remap_non_global);
	}
    }

  int remap_store_link_fast = TRUE;

  envstr = getenv("REMAP_STORE_LINK_FAST");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival >= 0 )
	{
	  remap_store_link_fast = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set REMAP_STORE_LINK_FAST to %d", remap_store_link_fast);
	}
    }

  remap_set_int(REMAP_STORE_LINK_FAST, remap_store_link_fast);

  envstr = getenv("REMAP_EXTRAPOLATE");
  if ( envstr )
    {
      if ( *envstr )
	{
	  if ( memcmp(envstr, "ON", 2) == 0 || memcmp(envstr, "on", 2) == 0 )
	    {
	      lextrapolate = true;
	      remap_extrapolate = true;
	    }
	  else if ( memcmp(envstr, "OFF", 3) == 0 || memcmp(envstr, "off", 3) == 0 )
	    {
	      lextrapolate = true;
	      remap_extrapolate = false;
	    }
	  else
	    cdoWarning("Environment variable REMAP_EXTRAPOLATE has wrong value!");

	  if ( cdoVerbose )
	    {
	      if ( remap_extrapolate )
		cdoPrint("Extrapolation enabled!");
	      else
		cdoPrint("Extrapolation disabled!");
	    }
	}
    }

  envstr = getenv("CDO_REMAP_GENWEIGHTS");
  if ( envstr )
    {
      if ( *envstr )
	{
	  if ( memcmp(envstr, "ON", 2) == 0 || memcmp(envstr, "on", 2) == 0 )
	    {
	      remap_genweights = TRUE;
	    }
	  else if ( memcmp(envstr, "OFF", 3) == 0 || memcmp(envstr, "off", 3) == 0 )
	    {
	      remap_genweights = FALSE;
	    }
	  else
	    cdoWarning("Environment variable CDO_REMAP_GENWEIGHTS has wrong value!");

	  if ( cdoVerbose )
	    {
	      if ( remap_genweights == TRUE )
		cdoPrint("Generation of weights enabled!");
	      else
		cdoPrint("Generation of weights disabled!");
	    }
	}
    }
}

static
void set_halo_to_missval(size_t nx, size_t ny, double *array, double missval)
{
  for ( size_t j = 0; j < ny+4; j++ ) array[j*(nx+4)+0]      = missval;
  for ( size_t j = 0; j < ny+4; j++ ) array[j*(nx+4)+1]      = missval;
  for ( size_t j = 0; j < ny+4; j++ ) array[j*(nx+4)+nx+2]   = missval;
  for ( size_t j = 0; j < ny+4; j++ ) array[j*(nx+4)+nx+3]   = missval;
  for ( size_t i = 0; i < nx+4; i++ ) array[     0*(nx+4)+i] = missval;
  for ( size_t i = 0; i < nx+4; i++ ) array[     1*(nx+4)+i] = missval;
  for ( size_t i = 0; i < nx+4; i++ ) array[(ny+2)*(nx+4)+i] = missval;
  for ( size_t i = 0; i < nx+4; i++ ) array[(ny+3)*(nx+4)+i] = missval;
}

static
bool is_global_grid(int gridID)
{
  bool global_grid = true;
  bool non_global = remap_non_global || !gridIsCircular(gridID);

  int gridtype = gridInqType(gridID);
  int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
  if ( (projtype == CDI_PROJ_RLL)  ||
       (projtype == CDI_PROJ_LAEA) ||
       (projtype == CDI_PROJ_SINU) ||
       (projtype == CDI_PROJ_LCC)  ||
       (gridtype == GRID_LONLAT && non_global) ||
       (gridtype == GRID_CURVILINEAR && non_global) ) global_grid = false;

  return global_grid;
}

static
void scale_gridbox_area(size_t gridsize, const double *restrict array1, size_t gridsize2, double *restrict array2, const double *restrict grid2_area)
{
  double array1sum = 0;
  for ( size_t i = 0; i < gridsize; i++ ) array1sum += array1[i];

  double array2sum = 0;
  for ( size_t i = 0; i < gridsize2; i++ ) array2sum += grid2_area[i];

  for ( size_t i = 0; i < gridsize2; i++ ) array2[i] = grid2_area[i]/array2sum*array1sum;

  static bool lgridboxinfo = true;
  if ( lgridboxinfo )
    {
      cdoPrint("gridbox_area replaced and scaled to %g", array1sum);
      lgridboxinfo = false;
    }
}

static
int set_remapgrids(int filetype, int vlistID, int ngrids, std::vector<bool>& remapgrids)
{
  int index;
  for ( index = 0; index < ngrids; index++ )
    {
      remapgrids[index] = true;

      int gridID = vlistGrid(vlistID, index);
      int gridtype = gridInqType(gridID);
      int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
      bool lproj4param = (gridtype == GRID_PROJECTION) && grid_has_proj4param(gridID);
  
      if ( gridtype != GRID_LONLAT      &&
           !lproj4param                 &&
           projtype != CDI_PROJ_RLL     &&
           projtype != CDI_PROJ_LAEA    &&
           projtype != CDI_PROJ_SINU    &&
           projtype != CDI_PROJ_LCC     &&
	   gridtype != GRID_GAUSSIAN    &&
	   gridtype != GRID_GME         &&
	   gridtype != GRID_CURVILINEAR &&
	   gridtype != GRID_UNSTRUCTURED )
	{
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    {
	      if ( !cdoRegulargrid && filetype == CDI_FILETYPE_GRB )
		cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular Gaussian grid!", gridNamePtr(gridtype));
	      else
		cdoAbort("Unsupported grid type: %s, use CDO operator -setgridtype,regular to convert reduced to regular Gaussian grid!", gridNamePtr(gridtype));
	    }
	  else if ( gridtype == GRID_GENERIC && gridInqSize(gridID) <= 2 )
            {
              remapgrids[index] = false;
            }
	  else
            {
              char varname[CDI_MAX_NAME];
              int nvars = vlistNvars(vlistID);
              for ( int varID = 0; varID < nvars; ++varID )
                if ( gridID == vlistInqVarGrid(vlistID, varID) )
                  {
                    vlistInqVarName(vlistID, varID, varname);
                    break;
                  }
              cdoAbort("Unsupported %s coordinates (Variable: %s)!", gridNamePtr(gridtype), varname);
            }
	}
    }

  for ( index = 0; index < ngrids; index++ )
    if ( remapgrids[index] ) break;

  if ( index == ngrids ) cdoAbort("No remappable grid found!");

  return index;
}

static
int set_max_remaps(int vlistID)
{
  int max_remaps = 0;

  const int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; index++ )
    {
      const int zaxisID = vlistZaxis(vlistID, index);
      const int zaxissize = zaxisInqSize(zaxisID);
      if ( zaxissize > max_remaps ) max_remaps = zaxissize;
    }
  
  const int nvars = vlistNvars(vlistID);
  if ( nvars > max_remaps ) max_remaps = nvars;

  max_remaps++;

  if ( cdoVerbose ) cdoPrint("Set max_remaps to %d", max_remaps);

  return max_remaps;
}

static
int get_norm_opt(void)
{
  int norm_opt = NORM_OPT_FRACAREA;

  char *envstr = getenv("CDO_REMAP_NORMALIZE_OPT"); // obsolate
  if ( envstr && *envstr )
    {
      if      ( memcmp(envstr, "frac", 4) == 0 ) norm_opt = NORM_OPT_FRACAREA;
      else if ( memcmp(envstr, "dest", 4) == 0 ) norm_opt = NORM_OPT_DESTAREA;
      else if ( memcmp(envstr, "none", 4) == 0 ) norm_opt = NORM_OPT_NONE;
      else cdoWarning("CDO_REMAP_NORMALIZE_OPT=%s unsupported!", envstr);
    }

  envstr = getenv("CDO_REMAP_NORM");
  if ( envstr && *envstr )
    {
      if      ( memcmp(envstr, "frac", 4) == 0 ) norm_opt = NORM_OPT_FRACAREA;
      else if ( memcmp(envstr, "dest", 4) == 0 ) norm_opt = NORM_OPT_DESTAREA;
      else if ( memcmp(envstr, "none", 4) == 0 ) norm_opt = NORM_OPT_NONE;
      else cdoWarning("CDO_REMAP_NORM=%s unsupported!", envstr);
    }

  if ( cdoVerbose )
    {
      if      ( norm_opt == NORM_OPT_FRACAREA ) cdoPrint("Normalization option: frac");
      else if ( norm_opt == NORM_OPT_DESTAREA ) cdoPrint("Normalization option: dest");
      else                                      cdoPrint("Normalization option: none");
    }

  return norm_opt;
}

static
void remap_normalize(int norm_opt, size_t gridsize, double *array, double missval, remapgrid_t *tgt_grid)
{
  /* used only to check the result of remapcon */
  double grid_err;

  if ( norm_opt == NORM_OPT_NONE )
    {
      for ( size_t i = 0; i < gridsize; i++ )
	{
	  if ( !DBL_IS_EQUAL(array[i], missval) )
	    {
	      grid_err = tgt_grid->cell_frac[i]*tgt_grid->cell_area[i];

	      if ( fabs(grid_err) > 0 )
		array[i] /= grid_err;
	      else
		array[i] = missval;
	    }
	}
    }
  else if ( norm_opt == NORM_OPT_DESTAREA )
    {
      for ( size_t i = 0; i < gridsize; i++ )
	{
	  if ( !DBL_IS_EQUAL(array[i], missval) )
	    {
	      if ( fabs(tgt_grid->cell_frac[i]) > 0 )
		array[i] /= tgt_grid->cell_frac[i];
	      else
		array[i] = missval;
	    }
	}
    }
}

static
void remap_set_frac_min(size_t gridsize, double *array, double missval, remapgrid_t *tgt_grid)
{
  if ( remap_frac_min > 0 )
    {
      for ( size_t i = 0; i < gridsize; i++ )
	{
	  //printf("%zu %g %g\n", i, remaps[r].tgt_grid.cell_frac[i], remaps[r].tgt_grid.cell_area[i]);
	  if ( tgt_grid->cell_frac[i] < remap_frac_min ) array[i] = missval;
	}
    }
}


int timer_remap, timer_remap_init, timer_remap_sort;
int timer_remap_bil, timer_remap_bic, timer_remap_dis, timer_remap_con, timer_remap_con_l1, timer_remap_con_l2;

static
void init_remap_timer(void)
{
  timer_remap        = timer_new("remap");
  timer_remap_init   = timer_new("remap init");
  timer_remap_sort   = timer_new("remap sort");
  timer_remap_bil    = timer_new("remap bil");
  timer_remap_bic    = timer_new("remap bic");
  timer_remap_dis    = timer_new("remap dis");
  timer_remap_con    = timer_new("remap con");
  timer_remap_con_l1 = timer_new("remap con loop1");
  timer_remap_con_l2 = timer_new("remap con loop2");
}

static
void links_per_value(remapvars_t *remapvars)
{
  long lpv = -1;

  remapvars->links_per_value = lpv;
  return;
  /*
  size_t num_links = remapvars->num_links;

  if ( num_links > 0 )
    {
      lpv = 1;
      const size_t *restrict dst_add = remapvars->tgt_cell_add;
      size_t n = 0;
      size_t ival = dst_add[n];
      for ( n = 1; n < num_links; ++n )
        if ( dst_add[n] == ival ) lpv++;
        else break;

    printf("lpv %zu\n", lpv);

      if ( num_links%lpv != 0 ) lpv = 0;

      n++;
      if ( n < num_links )
        {
          size_t lpv2 = 0;
          size_t ival2 = dst_add[n];
          for ( ; n < num_links; ++n )
            if ( dst_add[n] == ival2 ) lpv2++;
            else if ( lpv == lpv2 )
              {
                lpv2 = 0;
                ival2 = dst_add[n];
              }
            else
              {
                lpv = 0;
                break;
              }
        }
      
  printf("lpv %zu\n", lpv);

      if ( lpv == 1 )
        {
          for ( size_t n = 1; n < num_links; ++n )
            {
              if ( dst_add[n] == dst_add[n-1] )
                {
                  lpv = 0;
                  break;
                }
            }
        }
      else if ( lpv > 1 )
        {
          for ( size_t n = 1; n < num_links/lpv; ++n )
            {
              ival = dst_add[n*lpv];
              for ( size_t k = 1; k < lpv; ++k )
                {
                  if ( dst_add[n*lpv+k] != ival )
                    {
                      lpv = 0;
                      break;
                    }
                }
              if ( lpv == 0 ) break;
            }
        }
    }

  printf("lpv %zu\n", lpv);

  remapvars->links_per_value = lpv;
  */
}

static
void sort_remap_add(remapvars_t *remapvars)
{
  if ( cdoTimer ) timer_start(timer_remap_sort);
  if ( sort_mode == MERGE_SORT )
    { /* 
      ** use a combination of the old sort_add and a split and merge approach.
      ** The chunk size is determined by MERGE_SORT_LIMIT_SIZE in remaplib.c. 
      ** OpenMP parallelism is supported
      */   
      sort_iter(remapvars->num_links, remapvars->num_wts,
		remapvars->tgt_cell_add, remapvars->src_cell_add,
		remapvars->wts, ompNumThreads);
    }
  else
    { /* use a pure heap sort without any support of parallelism */
      sort_add(remapvars->num_links, remapvars->num_wts,
	       remapvars->tgt_cell_add, remapvars->src_cell_add,
	       remapvars->wts);
    }
  if ( cdoTimer ) timer_stop(timer_remap_sort);
}


void *Remap(void *argument)
{
  int streamID2 = -1;
  int nrecs;
  int varID, levelID;
  size_t gridsize, gridsize2;
  int gridID1 = -1, gridID2;
  int nmiss1, nmiss2;
  size_t i, j;
  int r = -1;
  int nremaps = 0;
  int norm_opt = NORM_OPT_NONE;
  int map_type = -1;
  int submap_type = SUBMAP_TYPE_NONE;
  int num_neighbors = 0;
  int need_gradiants = FALSE;
  char varname[CDI_MAX_NAME];
  double missval;
  remap_t *remaps = NULL;
  char *remap_file = NULL;

  if ( cdoTimer ) init_remap_timer();

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("remapcon",     REMAPCON,     0, NULL);
  cdoOperatorAdd("remapcon2",    REMAPCON2,    0, NULL);
  cdoOperatorAdd("remapbil",     REMAPBIL,     0, NULL);
  cdoOperatorAdd("remapbic",     REMAPBIC,     0, NULL);
  cdoOperatorAdd("remapdis",     REMAPDIS,     0, NULL);
  cdoOperatorAdd("remapnn",      REMAPNN,      0, NULL);
  cdoOperatorAdd("remaplaf",     REMAPLAF,     0, NULL);
  cdoOperatorAdd("remapsum",     REMAPSUM,     0, NULL);
  cdoOperatorAdd("gencon",       GENCON,       1, NULL);
  cdoOperatorAdd("gencon2",      GENCON2,      1, NULL);
  cdoOperatorAdd("genbil",       GENBIL,       1, NULL);
  cdoOperatorAdd("genbic",       GENBIC,       1, NULL);
  cdoOperatorAdd("gendis",       GENDIS,       1, NULL);
  cdoOperatorAdd("gennn",        GENNN,        1, NULL);
  cdoOperatorAdd("genlaf",       GENLAF,       1, NULL);
  cdoOperatorAdd("remap",        REMAPXXX,     0, NULL);
  cdoOperatorAdd("remapycon",    REMAPYCON,    0, NULL);
  cdoOperatorAdd("genycon",      GENYCON,      1, NULL);
  // clang-format on

  int operatorID   = cdoOperatorID();
  int operfunc     = cdoOperatorF1(operatorID);
  int lwrite_remap = cdoOperatorF2(operatorID);
  int lremapxxx    = operfunc == REMAPXXX;

  remap_set_int(REMAP_WRITE_REMAP, lwrite_remap);

  if ( operfunc == REMAPDIS || operfunc == GENDIS ||
       operfunc == REMAPNN  || operfunc == GENNN )
    remap_extrapolate = true;

  get_remap_env();

  if ( cdoVerbose )
    {
      if ( remap_extrapolate )
	cdoPrint("Extrapolation enabled!");
      else
	cdoPrint("Extrapolation disabled!");
    }

  // open stream before calling cdoDefineGrid!!!
  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  if ( lremapxxx )
    {
      operatorInputArg("grid description file or name, remap weights file (SCRIP NetCDF)");
      operatorCheckArgc(2);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
      remap_file = operatorArgv()[1];
    }
  else
    {
      operatorInputArg("grid description file or name");
      if ( operfunc == REMAPDIS && operatorArgc() == 2 )
        {
          gridID2 = cdoDefineGrid(operatorArgv()[0]);
          int inum = parameter2int(operatorArgv()[1]);
          //  if ( inum < 1 || inum > 9 ) cdoAbort("Number of nearest neighbors out of range (1-9)!", inum);
          num_neighbors = inum;
        }
      else
        {
          operatorCheckArgc(1);
          gridID2 = cdoDefineGrid(operatorArgv()[0]);
        }
    }

  if ( gridInqType(gridID2) == GRID_GENERIC ) cdoAbort("Unsupported target grid type (generic)!");

  int filetype = pstreamInqFiletype(streamID1);

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  std::vector<bool> remapgrids(ngrids);
  int index = set_remapgrids(filetype, vlistID1, ngrids, remapgrids);
  gridID1 = vlistGrid(vlistID1, index);

  for ( index = 0; index < ngrids; index++ )
    if ( remapgrids[index] )
      vlistChangeGridIndex(vlistID2, index, gridID2);

  if ( max_remaps == -1 ) max_remaps = set_max_remaps(vlistID1);

  if ( max_remaps > 0 )
    {
      remaps = (remap_t*) Malloc(max_remaps*sizeof(remap_t));
      for ( r = 0; r < max_remaps; r++ )
	{
	  remaps[r].nused    = 0;
	  remaps[r].gridID   = -1;
	  remaps[r].gridsize = 0;
	  remaps[r].nmiss    = 0;
	}
    }

  if ( lwrite_remap || lremapxxx )
    remap_genweights = TRUE;

  if ( lremapxxx )
    {
      size_t gridsize2;

      read_remap_scrip(remap_file, gridID1, gridID2, &map_type, &submap_type, &num_neighbors,
		       &remap_order, &remaps[0].src_grid, &remaps[0].tgt_grid, &remaps[0].vars);

      if ( remaps[0].vars.links_per_value == 0 ) links_per_value(&remaps[0].vars);
            
      nremaps = 1;
      gridsize = remaps[0].src_grid.size;
      remaps[0].gridID   = gridID1;
      remaps[0].gridsize = gridInqSize(gridID1);

      if ( map_type == MAP_TYPE_DISTWGT && !lextrapolate ) remap_extrapolate = true;
      if ( gridIsCircular(gridID1)      && !lextrapolate ) remap_extrapolate = true;

      if ( map_type == MAP_TYPE_DISTWGT && !remap_extrapolate && gridInqSize(gridID1) > 1 && !is_global_grid(gridID1) )
	{
	  remaps[0].gridsize += 4*(gridInqXsize(gridID1)+2) + 4*(gridInqYsize(gridID1)+2);
	  remaps[0].src_grid.non_global = true;
	}

      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].src_grid.nvgp;

      if ( gridsize != remaps[0].gridsize )
	cdoAbort("Size of source grid and weights from %s differ!", remap_file);

      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].src_grid.size;

      for ( i = 0; i < gridsize; i++ )
        if ( remaps[0].src_grid.mask[i] == FALSE )
          remaps[0].nmiss++;

      gridsize2 = gridInqSize(gridID2);
      if ( gridInqType(gridID2) == GRID_GME )
	{
	  remaps[0].tgt_grid.nvgp = gridInqSize(gridID2);
	  remaps[0].tgt_grid.vgpm = (int*) Realloc(remaps[0].tgt_grid.vgpm, gridInqSize(gridID2)*sizeof(int));
	  int gridID2_gme = gridToUnstructured(gridID2, 1);
	  gridInqMaskGME(gridID2_gme, remaps[0].tgt_grid.vgpm);
	  size_t isize = 0;
	  for ( size_t i = 0; i < gridsize2; ++i )
	    if ( remaps[0].tgt_grid.vgpm[i] ) isize++;
	  gridsize2 = isize;
	}
      /*
      printf("grid2 %zu %d %zu\n", gridsize2, remaps[0].tgt_grid.nvgp, remaps[0].tgt_grid.size);
      */
      if ( remaps[0].tgt_grid.size != gridsize2 )
	cdoAbort("Size of target grid and weights from %s differ!", remap_file);

      operfunc = maptype2operfunc(map_type, submap_type, num_neighbors, remap_order);

      if ( remap_test ) reorder_links(&remaps[0].vars);
    }
  else
    {
      get_map_type(operfunc, &map_type, &submap_type, &num_neighbors, &remap_order);
    }

  if ( remap_genweights == FALSE &&
       map_type != MAP_TYPE_BILINEAR && map_type != MAP_TYPE_BICUBIC &&
       map_type != MAP_TYPE_DISTWGT  && map_type != MAP_TYPE_CONSERV_YAC )
    remap_genweights = TRUE;

  remap_set_int(REMAP_GENWEIGHTS, remap_genweights);

  if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSERV_YAC ) norm_opt = get_norm_opt();

  size_t grid1sizemax = vlistGridsizeMax(vlistID1);

  if ( map_type == MAP_TYPE_BICUBIC ) need_gradiants = TRUE;
  if ( map_type == MAP_TYPE_CONSERV && remap_order == 2 )
    {
      if ( cdoVerbose ) cdoPrint("Second order remapping");
      need_gradiants = TRUE;
    }
  else
    remap_order = 1;

  double *grad1_lat    = need_gradiants ? (double*) Malloc(grid1sizemax*sizeof(double)) : NULL;
  double *grad1_lon    = need_gradiants ? (double*) Malloc(grid1sizemax*sizeof(double)) : NULL;
  double *grad1_latlon = need_gradiants ? (double*) Malloc(grid1sizemax*sizeof(double)) : NULL;

  double *array1 = (double*) Malloc(grid1sizemax*sizeof(double));
  int *imask = (int*) Malloc(grid1sizemax*sizeof(int));

  gridsize = gridInqSize(gridID2);
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

  if ( ! lwrite_remap )
    {
      streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
      pstreamDefVlist(streamID2, vlistID2);
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      if ( ! lwrite_remap ) 
	pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss1);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);

	  if ( !remapgrids[vlistGridIndex(vlistID1, gridID1)] )
	    {
	      if ( lwrite_remap ) continue;
	      else
		{
		  nmiss2 = nmiss1;
		  *array2 = *array1;
		  goto SKIPVAR;
		}
	    }

	  if ( map_type != MAP_TYPE_CONSERV && map_type != MAP_TYPE_CONSERV_YAC && 
	       gridInqType(gridID1) == GRID_GME && gridInqType(gridID2) == GRID_GME )
	    cdoAbort("Only conservative remapping is available to remap between GME grids!");

	  missval = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(gridID1);

	  if ( gridIsCircular(gridID1) && !lextrapolate ) remap_extrapolate = true;
	  if ( map_type == MAP_TYPE_DISTWGT && !remap_extrapolate && gridInqSize(gridID1) > 1 && !is_global_grid(gridID1) )
	    {
	      long nx = gridInqXsize(gridID1);
	      long ny = gridInqYsize(gridID1);
	      size_t gridsize_new = gridsize + 4*(nx+2) + 4*(ny+2);
	      if ( gridsize_new > grid1sizemax )
		{
		  grid1sizemax = gridsize_new;
		  array1 = (double*) Realloc(array1, grid1sizemax*sizeof(double));
		  imask  = (int*) Realloc(imask, grid1sizemax*sizeof(int));

		  if ( need_gradiants )
		    {
		      grad1_lat    = (double*) Realloc(grad1_lat, grid1sizemax*sizeof(double));
		      grad1_lon    = (double*) Realloc(grad1_lon, grid1sizemax*sizeof(double));
		      grad1_latlon = (double*) Realloc(grad1_latlon, grid1sizemax*sizeof(double));
		    }
		}
	      
	      for ( long j = ny-1; j >= 0; j-- )
		for ( long i = nx-1; i >= 0; i-- )
		  array1[(j+2)*(nx+4)+i+2] = array1[j*nx+i];

	      set_halo_to_missval(nx, ny, array1, missval);

	      gridsize = gridsize_new;
	      nmiss1 += 4*(nx+2) + 4*(ny+2);
	    }

	  for ( i = 0; i < gridsize; i++ )
            imask[i] = DBL_IS_EQUAL(array1[i], missval) ? FALSE : TRUE;

	  for ( r = nremaps-1; r >= 0; r-- )
	    {
	      if ( gridID1 == remaps[r].gridID && nmiss1 == remaps[r].nmiss )
		{
		  if ( memcmp(imask, remaps[r].src_grid.mask, remaps[r].src_grid.size*sizeof(int)) == 0 )
                    {
                      remaps[r].nused++;
                      break;
                    }
                }
	    }

	  if ( cdoVerbose && r >= 0 ) cdoPrint("Using remap %d", r);

	  if ( r < 0 )
	    {
	      if ( nremaps < max_remaps )
		{
		  r = nremaps;
		  nremaps++;
		}
	      else
		{
                  int n0 = 0;
                  if ( max_remaps > 1 && remaps[0].nused > remaps[1].nused ) n0 = 1;
                  remapVarsFree(&remaps[n0].vars);
                  remapGridFree(&remaps[n0].src_grid);
                  remapGridFree(&remaps[n0].tgt_grid);
                  for ( r = n0+1; r < nremaps; r++ ) memcpy(&remaps[r-1], &remaps[r], sizeof(remap_t));
                  r = nremaps - 1;
                  remaps[r].nused    = 0;
                  remaps[r].gridID   = -1;
                  remaps[r].gridsize = 0;
                  remaps[r].nmiss    = 0;
		}

	      if ( remaps[r].gridID != gridID1 )
		{
		  if ( gridIsCircular(gridID1) && !lextrapolate ) remap_extrapolate = true;
		  remaps[r].src_grid.non_global = false;
		  if ( map_type == MAP_TYPE_DISTWGT && !remap_extrapolate && gridInqSize(gridID1) > 1 && !is_global_grid(gridID1) )
		    {
		      remaps[r].src_grid.non_global = true;
		    }
		  /*
		    remaps[r].src_grid.luse_cell_area = FALSE;
		    remaps[r].tgt_grid.luse_cell_area = FALSE;
		  */
		  if ( gridInqType(gridID1) != GRID_UNSTRUCTURED && lremap_num_srch_bins == false )
		    {
		      if ( !remap_extrapolate && map_type == MAP_TYPE_DISTWGT )
			{
			  remap_num_srch_bins = 1;
			}
		      else
			{
			  int maxbins = 720;
			  int ysize1 = gridInqYsize(gridID1);
			  remap_num_srch_bins = ysize1/2 + ysize1%2;
			  if ( remap_num_srch_bins > maxbins ) remap_num_srch_bins = maxbins;
			  if ( remap_num_srch_bins < 1 )       remap_num_srch_bins = 1;
			}
		    }

		  remap_set_int(REMAP_NUM_SRCH_BINS, remap_num_srch_bins);

		  remaps[r].vars.norm_opt = norm_opt;
		  remaps[r].vars.pinit = false;
		  
		  if ( (map_type == MAP_TYPE_BILINEAR || map_type == MAP_TYPE_BICUBIC) &&
		       (gridInqType(gridID1) == GRID_GME || gridInqType(gridID1) == GRID_UNSTRUCTURED) )
		    cdoAbort("Bilinear/bicubic interpolation doesn't support unstructured source grids!");

		  /* initialize grid information for both grids */
		  if ( cdoTimer ) timer_start(timer_remap_init);
		  remap_grids_init(map_type, remap_extrapolate, gridID1, &remaps[r].src_grid, gridID2, &remaps[r].tgt_grid);
		  if ( cdoTimer ) timer_stop(timer_remap_init);
		}

	      remaps[r].gridID = gridID1;
	      remaps[r].nmiss  = nmiss1;

	      if ( gridInqType(gridID1) == GRID_GME )
		{
		  j = 0;
		  for ( i = 0; i < gridsize; i++ )
		    if ( remaps[r].src_grid.vgpm[i] ) imask[j++] = imask[i];
		}

	      memcpy(remaps[r].src_grid.mask, imask, remaps[r].src_grid.size*sizeof(int));

	      if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSERV_YAC )
		{
		  memset(remaps[r].src_grid.cell_area, 0, remaps[r].src_grid.size*sizeof(double));
		  memset(remaps[r].src_grid.cell_frac, 0, remaps[r].src_grid.size*sizeof(double));
		  memset(remaps[r].tgt_grid.cell_area, 0, remaps[r].tgt_grid.size*sizeof(double));
		}
	      memset(remaps[r].tgt_grid.cell_frac, 0, remaps[r].tgt_grid.size*sizeof(double));

	      /* initialize some remapping variables */
	      if ( cdoTimer ) timer_start(timer_remap_init);
	      remap_vars_init(map_type, remaps[r].src_grid.size, remaps[r].tgt_grid.size, &remaps[r].vars);
	      if ( cdoTimer ) timer_stop(timer_remap_init);

              print_remap_info(operfunc, remap_genweights, &remaps[r].src_grid, &remaps[r].tgt_grid, nmiss1);

	      if ( remap_genweights )
		{
		  if      ( map_type == MAP_TYPE_CONSERV     ) scrip_remap_conserv_weights(&remaps[r].src_grid, &remaps[r].tgt_grid, &remaps[r].vars);
		  else if ( map_type == MAP_TYPE_BILINEAR    ) scrip_remap_bilinear_weights(&remaps[r].src_grid, &remaps[r].tgt_grid, &remaps[r].vars);
		  else if ( map_type == MAP_TYPE_BICUBIC     ) scrip_remap_bicubic_weights(&remaps[r].src_grid, &remaps[r].tgt_grid, &remaps[r].vars);
		  else if ( map_type == MAP_TYPE_DISTWGT     ) remap_distwgt_weights(num_neighbors, &remaps[r].src_grid, &remaps[r].tgt_grid, &remaps[r].vars);
		  else if ( map_type == MAP_TYPE_CONSERV_YAC ) remap_conserv_weights(&remaps[r].src_grid, &remaps[r].tgt_grid, &remaps[r].vars);

		  if ( map_type == MAP_TYPE_CONSERV && remaps[r].vars.num_links != remaps[r].vars.max_links )
		    resize_remap_vars(&remaps[r].vars, remaps[r].vars.num_links-remaps[r].vars.max_links);
		  
		  if ( remaps[r].vars.sort_add ) sort_remap_add(&remaps[r].vars);
		  if ( remaps[r].vars.links_per_value == -1 ) links_per_value(&remaps[r].vars);

		  if ( lwrite_remap ) goto WRITE_REMAP;

		  if ( remap_test ) reorder_links(&remaps[r].vars);
		}
	    }

	  if ( gridInqType(gridID1) == GRID_GME )
	    {
	      j = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( remaps[r].src_grid.vgpm[i] ) array1[j++] = array1[i];
	    }
	  
	  if ( remap_genweights )
	    {
              remaps[r].nused++;

	      if ( need_gradiants )
		{
		  if ( remaps[r].src_grid.rank != 2 && remap_order == 2 )
		    cdoAbort("Second order remapping is not only available for unstructured grids!");

		  remap_gradients(remaps[r].src_grid, array1, grad1_lat, grad1_lon, grad1_latlon);
		}

	      if ( operfunc == REMAPLAF )
		remap_laf(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
			  remaps[r].vars.num_wts, remaps[r].vars.tgt_cell_add, remaps[r].vars.src_cell_add, array1);
	      else if ( operfunc == REMAPSUM )
		remap_sum(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
			  remaps[r].vars.num_wts, remaps[r].vars.tgt_cell_add, remaps[r].vars.src_cell_add, array1);
	      else
		remap(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		      remaps[r].vars.num_wts, remaps[r].vars.tgt_cell_add, remaps[r].vars.src_cell_add,
		      array1, grad1_lat, grad1_lon, grad1_latlon, remaps[r].vars.links, remaps[r].vars.links_per_value);
	    }
	  else
	    {
	      if      ( map_type == MAP_TYPE_BILINEAR    ) scrip_remap_bilinear(&remaps[r].src_grid, &remaps[r].tgt_grid, array1, array2, missval);
	      else if ( map_type == MAP_TYPE_BICUBIC     ) scrip_remap_bicubic(&remaps[r].src_grid, &remaps[r].tgt_grid, array1, array2, missval);
              else if ( map_type == MAP_TYPE_DISTWGT     ) remap_distwgt(num_neighbors, &remaps[r].src_grid, &remaps[r].tgt_grid, array1, array2, missval);
	      else if ( map_type == MAP_TYPE_CONSERV_YAC ) remap_conserv(&remaps[r].src_grid, &remaps[r].tgt_grid, array1, array2, missval);
	    }

	  gridsize2 = gridInqSize(gridID2);

	  if ( operfunc == REMAPCON || operfunc == REMAPCON2 || operfunc == REMAPYCON )
	    {
	      /* used only to check the result of remapcon */
	      if ( 0 ) remap_normalize(remaps[r].vars.norm_opt, gridsize2, array2, missval, &remaps[r].tgt_grid);

	      remap_set_frac_min(gridsize2, array2, missval, &remaps[r].tgt_grid);
	    }

	  if ( operfunc == REMAPSUM )
	    {
	      for ( i = 0; i < gridsize; i++ )
		printf("1 %zd %g %g %g %g\n", i, array1[i], remaps[r].src_grid.cell_frac[i], remaps[r].src_grid.cell_area[i],remaps[r].src_grid.cell_frac[i]);
	      double array1sum = 0;
	      for ( i = 0; i < gridsize; i++ )
		array1sum += remaps[r].src_grid.cell_area[i];

	      for ( i = 0; i < gridsize2; i++ )
		printf("2 %zd %g %g %g %g\n", i, array2[i], remaps[r].tgt_grid.cell_frac[i],remaps[r].tgt_grid.cell_area[i],remaps[r].tgt_grid.cell_frac[i]);
	      double array2sum = 0;
	      for ( i = 0; i < gridsize2; i++ )
		array2sum += remaps[r].tgt_grid.cell_area[i];

	      printf("array1sum %g, array2sum %g\n", array1sum, array2sum);
	    }

	  vlistInqVarName(vlistID1, varID, varname);
	  if ( operfunc == REMAPCON || operfunc == REMAPCON2 || operfunc == REMAPYCON )
	    if ( strcmp(varname, "gridbox_area") == 0 )
	      {
		scale_gridbox_area(gridsize, array1, gridsize2, array2, remaps[r].tgt_grid.cell_area);
	      }

	  /* calculate some statistics */
	  if ( cdoVerbose )
	    remap_stat(remap_order, remaps[r].src_grid, remaps[r].tgt_grid, remaps[r].vars, array1, array2, missval);

	  if ( gridInqType(gridID2) == GRID_GME )
	    {
	      int nd, ni, ni2, ni3;
 	      gridInqParamGME(gridID2, &nd, &ni, &ni2, &ni3);
	      j = remaps[r].tgt_grid.size;

	      for ( i = gridsize2; i > 0 ; i-- )
		if ( remaps[r].tgt_grid.vgpm[i-1] ) array2[i-1] = array2[--j];

	      gme_grid_restore(array2, ni, nd);
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize2; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss2++;

	SKIPVAR:

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, array2, nmiss2);
	}

      tsID++;
    }

  pstreamClose(streamID2);

 WRITE_REMAP:
 
  if ( lwrite_remap ) 
    write_remap_scrip(cdoStreamName(1)->args, map_type, submap_type, num_neighbors, remap_order,
		      remaps[r].src_grid, remaps[r].tgt_grid, remaps[r].vars);

  pstreamClose(streamID1);

  if ( imask )  Free(imask);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  if ( grad1_latlon ) Free(grad1_latlon);
  if ( grad1_lon ) Free(grad1_lon);
  if ( grad1_lat ) Free(grad1_lat);

  if ( max_remaps > 0 )
    {
      if ( lremapxxx && remap_genweights && remaps[0].nused == 0 )
        print_remap_warning(remap_file, operfunc, &remaps[0].src_grid, remaps[0].nmiss);
      
      for ( r = 0; r < nremaps; r++ )
	{
	  remapVarsFree(&remaps[r].vars);
	  remapGridFree(&remaps[r].src_grid);
	  remapGridFree(&remaps[r].tgt_grid);
	}
      
      if ( remaps ) Free(remaps);
    }

  cdoFinish();

  return 0;
}
