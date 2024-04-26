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

   Output field with grid cell center or cell bounds for plotting with GMT

    - outputcenter
    - outputbounds
    - outputboundscpt
    - outputvector
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"
#include "color.h"

double intlin(double x, double y1, double x1, double y2, double x2);

static
int check_ncorner(int ncorner, const double *lon_bounds, const double *lat_bounds)
{
  int ncorner_new = ncorner;
  int k;

  for ( k = ncorner-1; k > 0; --k )
    if ( IS_NOT_EQUAL(lon_bounds[k], lon_bounds[k-1]) ||
	 IS_NOT_EQUAL(lat_bounds[k], lat_bounds[k-1]) ) break;

  if ( k < ncorner-1 ) ncorner_new = k+1;

  return ncorner_new;
}


void make_cyclic(double *array1, double *array2, int nlon, int nlat)
{
  int ij1, ij2;

  for ( int j = 0; j < nlat; ++j )
    for ( int i = 0; i < nlon; ++i )
      {
        ij1 = j*nlon+i;
        ij2 = j*(nlon+1)+i;
        array2[ij2] = array1[ij1];
      }

  for ( int j = 0; j < nlat; ++j )
    {
      ij2 = j*(nlon+1);
      array2[ij2+nlon] = array2[ij2];
    }
}

static
void array_stat(int ngp, double *restrict array, double missval, double *minval, double *maxval, double *meanval)
{
  double rmin = DBL_MAX;
  double rmax = -DBL_MAX;
  double rmean = 0;
  
  int nvals = 0;
  for ( int i = 0; i < ngp; i++ )
    {
      if ( !DBL_IS_EQUAL(array[i], missval) )
        {
          if ( array[i] < rmin ) rmin = array[i];
          if ( array[i] > rmax ) rmax = array[i];
          rmean += array[i];
          nvals++;
        }
    }

  if ( IS_EQUAL(rmin,  DBL_MAX) ) rmin = missval;
  if ( IS_EQUAL(rmax, -DBL_MAX) ) rmax = missval;

  if ( nvals ) rmean /= nvals;
  else         rmean = missval;

  *minval = rmin;
  *maxval = rmax;
  *meanval = rmean;
}

static
void output_vrml(int nlon, int nlat, int ngp, double *restrict array, double missval, CPT *cpt)
{
  double minval, maxval, meanval;
  array_stat(ngp, array, missval, &minval, &maxval, &meanval);

  double dx = 10./nlon;

  printf("Viewpoint {\n");
  printf("  description \"viewpoint1\"\n");
  printf("  orientation 0 0 1 0\n");
  printf("  position 0.0 0.0 10.0\n");
  printf("}\n");
  printf("\n");
  printf("Background {\n");
  printf("  skyColor [\n");
  printf("    0.0 0.1 0.8,\n");
  printf("    0.0 0.5 1.0,\n");
  printf("    1.0 1.0 1.0\n");
  printf("  ]\n");
  printf("  skyAngle [0.785, 1.571]\n");
  printf("\n");
  printf("  groundColor [\n");
  printf("    0.0 0.0 0.0,\n");
  printf("    0.3 0.3 0.3,\n");
  printf("    0.5 0.5 0.5\n");
  printf("  ]\n");
  printf("  groundAngle [0.785, 1.571]\n");
  printf("}\n");
  printf("\n");
  printf("Transform {\n");
  printf("  children [\n");
  printf("    Shape {\n");
  printf("      appearance Appearance {\n");
  printf("        material Material {}\n");
  printf("      }\n");
  printf("      geometry ElevationGrid {\n");
  printf("        colorPerVertex TRUE\n");
  printf("        solid FALSE\n");
  printf("        xDimension %d\n", nlon);
  printf("        zDimension %d\n", nlat);
  printf("        xSpacing %g\n", dx);
  printf("        zSpacing %g\n", dx);
  printf("        color Color {\n");
  printf("          color [\n");
  for ( int j = nlat-1; j >= 0 ; --j )
    for ( int i = 0; i < nlon; ++i )
      {
        int r = 0, g = 0, b = 0;
        double val = array[j*nlon+i];
        
        if ( !DBL_IS_EQUAL(val, missval) )
          {
            int n;
            for ( n = 0; n < cpt->ncolors; n++ )
              if ( val > cpt->lut[n].z_low && val <= cpt->lut[n].z_high ) break;
            
            if ( n == cpt->ncolors )
              {
                r = cpt->bfn[0].rgb[0];  g = cpt->bfn[0].rgb[1];  b = cpt->bfn[0].rgb[2];
              }
            else
              {
                //  r = cpt->lut[n].rgb_high[0];  g = cpt->lut[n].rgb_high[1];  b = cpt->lut[n].rgb_high[2];
                r = intlin(val, cpt->lut[n].rgb_low[0], cpt->lut[n].z_low, cpt->lut[n].rgb_high[0], cpt->lut[n].z_high);
                g = intlin(val, cpt->lut[n].rgb_low[1], cpt->lut[n].z_low, cpt->lut[n].rgb_high[1], cpt->lut[n].z_high);
                b = intlin(val, cpt->lut[n].rgb_low[2], cpt->lut[n].z_low, cpt->lut[n].rgb_high[2], cpt->lut[n].z_high);
              }
          }
        else
          {
            r = cpt->bfn[2].rgb[0];  g = cpt->bfn[2].rgb[1];  b = cpt->bfn[2].rgb[2]; 
          }
        printf(" %.3g %.3g %.3g,\n", r/255., g/255., b/255.);
      }
  printf("          ]\n");
  printf("        }\n");
  printf("        height [\n");

  for ( int j = nlat-1; j >= 0 ; --j )
    for ( int i = 0; i < nlon; ++i )
      printf("%g,\n", array[j*nlon+i]);

  printf("        ]\n");
  printf("      }\n");
  printf("    }\n");
  printf("  ]\n");
  printf("  translation -5 0 %g\n", -5.*nlat/nlon);
  printf("  rotation 0.0 0.0 0.0 0.0\n");
  printf("  scale 1.0 %g 1.0\n", 0.5/(maxval-minval));
  printf("}\n");
}


void *Outputgmt(void *argument)
{
  int varID0;
  int gridsize2 = 0;
  int nrecs;
  int levelID;
  int nmiss;
  int ninc = 1;
  bool lzon = false, lmer = false, lhov = false;
  bool lgrid_gen_bounds = false, luse_grid_corner = false;
  char varname[CDI_MAX_NAME];
  double *array2 = NULL;
  double *uf = NULL, *vf = NULL, *alpha = NULL, *auv = NULL;
  double *grid_center_lat2 = NULL, *grid_center_lon2 = NULL;
  double *grid_corner_lat = NULL, *grid_corner_lon = NULL;
  int *grid_mask = NULL;
  CPT cpt;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];	  

  cdoInitialize(argument);

  // clang-format off
  int OUTPUTCENTER    = cdoOperatorAdd("gmtxyz",          0, 0, NULL);
  int OUTPUTCENTER2   = cdoOperatorAdd("outputcenter2",   0, 0, NULL);
  int OUTPUTCENTERCPT = cdoOperatorAdd("outputcentercpt", 0, 0, NULL);
  int OUTPUTBOUNDS    = cdoOperatorAdd("gmtcells",        0, 0, NULL);
  int OUTPUTBOUNDSCPT = cdoOperatorAdd("outputboundscpt", 0, 0, NULL);
  int OUTPUTVECTOR    = cdoOperatorAdd("outputvector",    0, 0, NULL);
  int OUTPUTTRI       = cdoOperatorAdd("outputtri",       0, 0, NULL);
  int OUTPUTVRML      = cdoOperatorAdd("outputvrml",      0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTVECTOR )
    {
      operatorInputArg("increment");
      operatorCheckArgc(1);
      ninc = parameter2int(operatorArgv()[0]);
      if ( ninc < 1 ) cdoAbort("Increment must be greater than 0!");
    }

  if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
    luse_grid_corner = true;

  if ( operatorID == OUTPUTCENTERCPT || operatorID == OUTPUTBOUNDSCPT || operatorID == OUTPUTVRML )
    {
       operatorCheckArgc(1);
       char *cpt_file = operatorArgv()[0];

       FILE *cpt_fp = fopen(cpt_file, "r");
       if ( cpt_fp == NULL )
	cdoAbort("Open failed on color palette table %s", cpt_file);

      int status = cptRead(cpt_fp, &cpt);
      if ( status != 0 )
	cdoAbort("Error during read of color palette table %s", cpt_file);
      
      if ( cdoVerbose ) cptWrite(stderr, cpt);
    }

  int streamID = pstreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  int varID = 0;
  vlistInqVarName(vlistID, varID, varname);
  int code    = vlistInqVarCode(vlistID, varID);
  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  double missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 1);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    {
      gridID = gridToCurvilinear(gridID, 1);
      lgrid_gen_bounds = true;
    }

  int gridsize = gridInqSize(gridID);
  int nlon     = gridInqXsize(gridID);
  int nlat     = gridInqYsize(gridID);
  int nlev     = zaxisInqSize(zaxisID);

  if ( gridInqMaskGME(gridID, NULL) )
    {
      grid_mask = (int*) Malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
    }

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED )
    {
      if ( nlon == 1 && nlat  > 1 && nlev == 1 ) lhov = true;
      if ( nlon == 1 && nlat  > 1 && nlev  > 1 ) lzon = true;
      if ( nlon  > 1 && nlat == 1 && nlev  > 1 ) lmer = true;
    }
  else
    {
      nlat = 1;
    }

  if ( cdoVerbose && lhov ) cdoPrint("Process hovmoeller data");
  if ( cdoVerbose && lzon ) cdoPrint("Process zonal data");
  if ( cdoVerbose && lmer ) cdoPrint("Process meridional data");
  /*
  if ( lzon || lmer ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available for zonal/meridional data!");
    }
  */
  if ( lhov ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available hovmoeller data!");
    }

  int ncorner = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  bool grid_is_circular = gridIsCircular(gridID);

  double *grid_center_lat = (double*) Malloc(gridsize*sizeof(double));
  double *grid_center_lon = (double*) Malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");

  int nvals = gridsize;
  double *plon = grid_center_lon;
  double *plat = grid_center_lat;

  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      gridsize2 = nlat*(nlon+1);

      grid_center_lat2 = (double*) Malloc(gridsize2*sizeof(double));
      grid_center_lon2 = (double*) Malloc(gridsize2*sizeof(double));

      make_cyclic(grid_center_lat, grid_center_lat2, nlon, nlat);
      make_cyclic(grid_center_lon, grid_center_lon2, nlon, nlat);

      for ( int j = 0; j < nlat; ++j )
	{
	  int ij2 = j*(nlon+1);
	  grid_center_lon2[ij2+nlon] += 360;
	}

      nvals = gridsize2;
      plon = grid_center_lon2;
      plat = grid_center_lat2;
    }

  double *zaxis_center_lev = (double*) Malloc(nlev*sizeof(double));
  double *zaxis_lower_lev  = (double*) Malloc(nlev*sizeof(double));
  double *zaxis_upper_lev  = (double*) Malloc(nlev*sizeof(double));

  cdoZaxisInqLevels(zaxisID, zaxis_center_lev);

  if ( luse_grid_corner )
    {
      if ( ncorner == 0 ) cdoAbort("grid corner missing!");
      int nalloc = ncorner*gridsize;
      grid_corner_lat = (double*) Realloc(grid_corner_lat, nalloc*sizeof(double));
      grid_corner_lon = (double*) Realloc(grid_corner_lon, nalloc*sizeof(double));

      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqYbounds(gridID, grid_corner_lat);
	  gridInqXbounds(gridID, grid_corner_lon);
	}
      else
	{
	  if ( lgrid_gen_bounds )
	    {
	      char xunitstr[CDI_MAX_NAME];
	      char yunitstr[CDI_MAX_NAME];
	      gridInqXunits(gridID, xunitstr);
	      gridInqYunits(gridID, yunitstr);
	      if ( ! lzon ) grid_cell_center_to_bounds_X2D(xunitstr, nlon, nlat, grid_center_lon, grid_corner_lon, 0);
	      if ( ! lmer ) grid_cell_center_to_bounds_Y2D(yunitstr, nlon, nlat, grid_center_lat, grid_corner_lat);
	    }
	  else
	    cdoAbort("Grid corner missing!");
	}


      /* Note: using units from latitude instead from bounds */
      grid_to_degree(units, ncorner*gridsize, grid_corner_lon, "grid corner lon");
      grid_to_degree(units, ncorner*gridsize, grid_corner_lat, "grid corner lat");

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	{
	  zaxisInqLbounds(zaxisID, zaxis_lower_lev);
	  zaxisInqUbounds(zaxisID, zaxis_upper_lev);
	}
      else
	{
	  zaxis_lower_lev[0] = zaxis_center_lev[0];
	  for ( int i = 1; i < nlev; ++i )
	    zaxis_lower_lev[i] = 0.5*(zaxis_center_lev[i] + zaxis_center_lev[i-1]);

	  zaxis_upper_lev[nlev-1] = zaxis_center_lev[nlev-1];
	  for ( int i = 0; i < nlev-1; ++i )
	    zaxis_upper_lev[i] = zaxis_lower_lev[i+1];

	  if ( cdoVerbose )
	    for ( int i = 0; i < nlev; ++i )
	      fprintf(stderr, "level: %d %g %g %g\n",
		     i+1, zaxis_lower_lev[i], zaxis_center_lev[i], zaxis_upper_lev[i]);
	}
    }

  double *array = (double*) Malloc(gridsize*sizeof(double));
  double *parray = array;
						
  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      array2 = (double*) Malloc(nlat*(nlon+1)*sizeof(double));
      parray = array2;
    }

  if ( operatorID == OUTPUTVECTOR )
    {
      uf    = (double*) Malloc(gridsize*sizeof(double));
      vf    = (double*) Malloc(gridsize*sizeof(double));
      alpha = (double*) Malloc(gridsize*sizeof(double));
      auv   = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      if ( tsID == 0 && operatorID != OUTPUTTRI )
	{
	  if ( operatorID == OUTPUTVRML )
	    printf("#VRML V2.0 utf8\n\n");
#if defined(VERSION)
	  fprintf(stdout, "# Generated by CDO version %s\n", VERSION);
	  fprintf(stdout, "#\n");
#endif
	  fprintf(stdout, "# Operator = %s\n", cdoOperatorName(operatorID));
	  if      ( lhov )  fprintf(stdout, "# Mode     = hovmoeller\n");
	  else if ( lzon )  fprintf(stdout, "# Mode     = zonal\n");
	  else if ( lmer )  fprintf(stdout, "# Mode     = meridional\n");
	  else              fprintf(stdout, "# Mode     = horizonal\n");

	  if ( operatorID == OUTPUTVECTOR )
	    fprintf(stdout, "# Increment = %d\n", ninc);
	  fprintf(stdout, "#\n");
	  fprintf(stdout, "# File  = %s\n", cdoStreamName(0)->args);
	  fprintf(stdout, "# Date  = %s\n", vdatestr);
	  fprintf(stdout, "# Time  = %s\n", vtimestr);
	  fprintf(stdout, "# Name  = %s\n", varname);
	  fprintf(stdout, "# Code  = %d\n", code);
	}

      varID0 = varID;

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID, &varID, &levelID);

	  if ( varID != varID0 ) continue;
	  if ( recID > 0 && !lzon && !lmer ) continue;

	  pstreamReadRecord(streamID, array, &nmiss);

	  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
	    make_cyclic(array, array2, nlon, nlat);

	  double level = zaxis_center_lev[levelID];

	  if ( (tsID == 0 || lzon || lmer) && operatorID != OUTPUTTRI )
	    fprintf(stdout, "# Level = %g\n", level);
	  if ( lhov )
	    fprintf(stdout, "# Timestep = %d\n", tsID+1);

	  if ( operatorID != OUTPUTTRI ) fprintf(stdout, "#\n");

	  if ( operatorID == OUTPUTCENTER || operatorID == OUTPUTCENTER2 || operatorID == OUTPUTCENTERCPT )
	    {
              if ( cdoVerbose )
                {
                  double minval, maxval, meanval;
                  array_stat(gridsize, array, missval, &minval, &maxval, &meanval);
                  double range = maxval - minval;
                  fprintf(stderr, "makecpt -T%g/%g/%g -Crainbow > gmt.cpt\n", minval, maxval, range/20);
                  fprintf(stderr, "pscontour -K -JQ0/10i -Rd -I -Cgmt.cpt data.gmt > gmtplot.ps\n");
                  fprintf(stderr, "pscoast -O -J -R -Dc -W -B40g20 >> gmtplot.ps\n");
                }

	      for ( int i = 0; i < nvals; i++ )
		{
		  if ( grid_mask && grid_mask[i] == 0 ) continue;

		  if ( operatorID == OUTPUTCENTER )
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else if ( lhov )
			fprintf(stdout, " %d  %g  %g\n", tsID+1, grid_center_lat[i], array[i]);
		      else
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i]);
		    }
		  else if ( operatorID == OUTPUTCENTER2 )
		    {
		      fprintf(stdout, " %g  %g  %g\n", plon[i], plat[i], parray[i]);
		    }
		  else
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else
			fprintf(stdout, " %g  %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
		    }
		}
	      fprintf(stdout, "#\n");
	    }
	  else if ( operatorID == OUTPUTTRI )
	    {
	      int c1, c2, c3;
	      int ip1;
	      if ( gridInqType(gridID) != GRID_CURVILINEAR ) cdoAbort("Unsupported grid!");

	      int mlon = nlon-1;
	      /* if ( gridIsCircular(gridID) ) mlon = nlon; */
	      for ( int j = 0; j < nlat-1; ++j )
                for ( int i = 0; i < mlon; ++i )
                  {
                    ip1 = i+1;
                    if ( i == nlon-1 ) ip1 = 0;
                    c1 = (j)*nlon+ip1;
                    c2 = (j)*nlon+i;
                    c3 = (j+1)*nlon+i;
                    fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
                    c1 = (j)*nlon+i+1;
                    c2 = (j+1)*nlon+i;
                    c3 = (j+1)*nlon+ip1;
                    fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
                  }
	    }
	  else if ( operatorID == OUTPUTVECTOR )
	    {
	      if ( nrecs < 2 ) cdoAbort("Too few fields!");

	      memcpy(uf, array, gridsize*sizeof(double));
	      pstreamInqRecord(streamID, &varID, &levelID);
	      pstreamReadRecord(streamID, vf, &nmiss);

	      for ( int j = 0; j < nlat; j += ninc )
		for ( int i = 0; i < nlon; i += ninc )
		  {
		    /* compute length of velocity vector */
		    auv[IX2D(j,i,nlon)] = sqrt(uf[IX2D(j,i,nlon)]*uf[IX2D(j,i,nlon)] + 
					       vf[IX2D(j,i,nlon)]*vf[IX2D(j,i,nlon)]);

		    alpha[IX2D(j,i,nlon)] = atan2(vf[IX2D(j,i,nlon)],uf[IX2D(j,i,nlon)]);
		    alpha[IX2D(j,i,nlon)] = 90. - alpha[IX2D(j,i,nlon)]*RAD2DEG;

		    if ( alpha[IX2D(j,i,nlon)] <   0 ) alpha[IX2D(j,i,nlon)] += 360;
		    if ( alpha[IX2D(j,i,nlon)] > 360 ) alpha[IX2D(j,i,nlon)] -= 360;

		    if ( fabs(auv[IX2D(j,i,nlon)]) > 0 )
		      fprintf(stdout, " %g  %g  %g  %g\n",
			      grid_center_lon[IX2D(j,i,nlon)], grid_center_lat[IX2D(j,i,nlon)],
			      alpha[IX2D(j,i,nlon)], auv[IX2D(j,i,nlon)]);
		  }

	      fprintf(stdout, "#\n");
	      break;
	    }
	  else if ( operatorID == OUTPUTVRML )
	    {
              output_vrml(nlon, nlat, gridsize, array, missval, &cpt);
            }
	  else if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	    {
              if ( cdoVerbose )
                {
                  double minval, maxval, meanval;
                  array_stat(gridsize, array, missval, &minval, &maxval, &meanval);
                  double range = maxval - minval;
                  fprintf(stderr, "makecpt -T%g/%g/%g -Crainbow > gmt.cpt\n", minval, maxval, range/20);
                  fprintf(stderr, "psxy -K -JQ0/10i -Rd -L -Cgmt.cpt -m data.gmt > gmtplot.ps\n");
                  // fprintf(stderr, "psxy -K -Jx0.028id -Rd -L -Cgmt.cpt -m data.gmt > gmtplot.ps\n");
                  // fprintf(stderr, "psxy -K -JN0/10i -Rd -L -Cgmt.cpt -m data.gmt > gmtplot.ps\n");
                  fprintf(stderr, "pscoast -O -J -R -Dc -W -B40g20 >> gmtplot.ps\n");
                  fprintf(stderr, "ps2pdf gmtplot.ps\n");
                }

	      for ( int i = 0; i < gridsize; i++ )
		{
		  if ( grid_mask && grid_mask[i] == 0 ) continue;

		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    fprintf(stdout, "> -Z%g", array[i]);
		  else
		    fprintf(stdout, "> -ZNaN");

		  if ( operatorID == OUTPUTBOUNDSCPT )
		    {
		      int r = 0, g = 0, b = 0, n;

		      if ( !DBL_IS_EQUAL(array[i], missval) )
			{
			  for ( n = 0; n < cpt.ncolors; n++ )
			    if ( array[i] > cpt.lut[n].z_low && array[i] <= cpt.lut[n].z_high ) break;

			  if ( n == cpt.ncolors )
			    {
			      r = cpt.bfn[0].rgb[0];  g = cpt.bfn[0].rgb[1];  b = cpt.bfn[0].rgb[2];
			    }
			  else
			    {
			      r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
			    }
			}
		      else
			{
			  r = cpt.bfn[2].rgb[0];  g = cpt.bfn[2].rgb[1];  b = cpt.bfn[2].rgb[2]; 
			}

		      fprintf(stdout, " -G%d/%d/%d", r, g, b);
		    }

		  fprintf(stdout, "\n");

		  if ( lzon )
		    {
		      double xlev[4], xlat[4];
		      double levmin = zaxis_lower_lev[levelID];
		      double levmax = zaxis_upper_lev[levelID];
		      double latmin = grid_corner_lat[i*4];
		      double latmax = grid_corner_lat[i*4];
		      for ( int ic = 1; ic < 4; ic++ )
			{
			  if ( grid_corner_lat[i*4+ic] < latmin ) latmin = grid_corner_lat[i*4+ic];
			  if ( grid_corner_lat[i*4+ic] > latmax ) latmax = grid_corner_lat[i*4+ic];
			}
		      xlev[0] = levmin;
		      xlev[1] = levmax;
		      xlev[2] = levmax;
		      xlev[3] = levmin;
		      xlat[0] = latmin;
		      xlat[1] = latmin;
		      xlat[2] = latmax;
		      xlat[3] = latmax;
		      for ( int ic = 0; ic < 4; ic++ )
			fprintf(stdout, "   %g  %g\n", xlat[ic], xlev[ic]);
		      fprintf(stdout, "   %g  %g\n", xlat[0], xlev[0]);
		    }
		  else if ( lmer )
		    {
		      double xlev[4], xlon[4];
		      double levmin = zaxis_lower_lev[levelID];
		      double levmax = zaxis_upper_lev[levelID];
		      double lonmin = grid_corner_lon[i*4];
		      double lonmax = grid_corner_lon[i*4];
		      for ( int ic = 1; ic < 4; ic++ )
			{
			  if ( grid_corner_lon[i*4+ic] < lonmin ) lonmin = grid_corner_lon[i*4+ic];
			  if ( grid_corner_lon[i*4+ic] > lonmax ) lonmax = grid_corner_lon[i*4+ic];
			}
		      xlev[0] = levmin;
		      xlev[1] = levmin;
		      xlev[2] = levmax;
		      xlev[3] = levmax;
		      xlon[0] = lonmin;
		      xlon[1] = lonmax;
		      xlon[2] = lonmax;
		      xlon[3] = lonmin;
		      for ( int ic = 0; ic < 4; ic++ )
			fprintf(stdout, "   %g  %g\n", xlon[ic], xlev[ic]);
		      fprintf(stdout, "   %g  %g\n", xlon[0], xlev[0]);
		    }
		  else if ( lhov )
		    {
		      cdoAbort("Implementation for hovmoeller data missing!");
		    }
		  else
		    {
		      const double *lon_bounds = grid_corner_lon+i*ncorner;
		      const double *lat_bounds = grid_corner_lat+i*ncorner;
		      int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

		      for ( int ic = 0; ic < ncorner_new; ic++ )
			fprintf(stdout, "   %g  %g\n", lon_bounds[ic], lat_bounds[ic]);
		      fprintf(stdout, "   %g  %g\n", lon_bounds[0], lat_bounds[0]);
		    }
		}
	      fprintf(stdout, "\n");
	    }
	}

      if ( ! lhov ) break;

      tsID++;
    }

  pstreamClose(streamID);

  if ( array  ) Free(array);
  if ( array2 ) Free(array2);
  if ( grid_mask ) Free(grid_mask);
  if ( grid_center_lon ) Free(grid_center_lon);
  if ( grid_center_lat ) Free(grid_center_lat);
  if ( grid_center_lon2 ) Free(grid_center_lon2);
  if ( grid_center_lat2 ) Free(grid_center_lat2);
  if ( grid_corner_lon ) Free(grid_corner_lon);
  if ( grid_corner_lat ) Free(grid_corner_lat);

  Free(zaxis_center_lev);
  Free(zaxis_lower_lev);
  Free(zaxis_upper_lev);

  cdoFinish();

  return 0;
}
