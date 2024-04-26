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

#include <cdi.h>
#include "cdo_int.h"


double cdoZaxisInqLevel(int zaxisID, int levelID)
{
  int zaxistype  = zaxisInqType(zaxisID);
  double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) :
                 (zaxistype == ZAXIS_SURFACE) ? 0 : levelID+1;
  return level;
}


int cdoZaxisInqLevels(int zaxisID, double *levels)
{
  int size = zaxisInqLevels(zaxisID, NULL);

  if ( levels )
    {
      if ( size )
        zaxisInqLevels(zaxisID, levels);
      else
        {
          size = zaxisInqSize(zaxisID);
          if ( size == 1 && zaxisInqType(zaxisID) == ZAXIS_SURFACE )
            levels[0] = 0;
          else
            for ( int i = 0; i < size; ++i ) levels[i] = i+1;
        }
    }

  return size;
}

static
void compare_lat_reg2d(int ysize, int gridID1, int gridID2)
{
  if ( ysize > 1 )
    {      
      double *yvals1 = (double*) Malloc(ysize*sizeof(double));
      double *yvals2 = (double*) Malloc(ysize*sizeof(double));

      gridInqYvals(gridID1, yvals1);
      gridInqYvals(gridID2, yvals2);
		
      if ( IS_EQUAL(yvals1[0], yvals2[ysize-1]) &&
	   IS_EQUAL(yvals1[ysize-1], yvals2[0]) )
	{
	  if ( yvals1[0] > yvals2[0] )
	    cdoWarning("Latitude orientation differ! First grid: N->S; second grid: S->N");
	  else
	    cdoWarning("Latitude orientation differ! First grid: S->N; second grid: N->S");
	}
      else
	{
	  for ( int i = 0; i < ysize; ++i )
	    if ( fabs(yvals1[i] - yvals2[i]) > 3.e-5 )
	      {
		cdoWarning("Grid latitudes differ!");
		break;
	      }
	}

      Free(yvals1);
      Free(yvals2);
    }
}

static
void compare_lon_reg2d(int xsize, int gridID1, int gridID2)
{
  if ( xsize > 1 )
    {
      double *xvals1 = (double*) Malloc(xsize*sizeof(double));
      double *xvals2 = (double*) Malloc(xsize*sizeof(double));

      gridInqXvals(gridID1, xvals1);
      gridInqXvals(gridID2, xvals2);
		  
      for ( int i = 0; i < xsize; ++i )
	if ( fabs(xvals1[i] - xvals2[i]) > 3.e-5 )
	  {
	    cdoWarning("Grid longitudes differ!");
	    break;
	  }
      
      Free(xvals1);
      Free(xvals2);
    }
}

static
void compare_grid_unstructured(int gridID1, int gridID2)
{
  if ( gridInqXvals(gridID1, NULL) && gridInqXvals(gridID1, NULL) == gridInqXvals(gridID2, NULL) &&
       gridInqYvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) == gridInqYvals(gridID2, NULL) )
    {
      int gridsize = gridInqSize(gridID1);
      
      double *xvals1 = (double*) Malloc(gridsize*sizeof(double));
      double *xvals2 = (double*) Malloc(gridsize*sizeof(double));
      double *yvals1 = (double*) Malloc(gridsize*sizeof(double));
      double *yvals2 = (double*) Malloc(gridsize*sizeof(double));

      gridInqXvals(gridID1, xvals1);
      gridInqXvals(gridID2, xvals2);
      gridInqYvals(gridID1, yvals1);
      gridInqYvals(gridID2, yvals2);

      int inc = gridsize > 10000 ? gridsize/1000 : 1;
      for ( int i = 0; i < gridsize; i += inc )
        if ( fabs(xvals1[i] - xvals2[i]) > 2.e-5 || fabs(yvals1[i] - yvals2[i]) > 2.e-5 )
          {
            // printf("%d %g %g %g %g %g %g\n", i, xvals1[i], xvals2[i], yvals1[i], yvals2[i], xvals1[i] - xvals2[i], yvals1[i] - yvals2[i]);
            cdoWarning("Geographic location of some grid points differ!");
            break;
          }
      
      Free(xvals1);
      Free(xvals2);
      Free(yvals1); 
      Free(yvals2); 
    }
}


void cdoCompareGrids(int gridID1, int gridID2)
{
  /* compare grids of first variable */

  if ( gridInqType(gridID1) == gridInqType(gridID2) )
    {
      if ( gridInqType(gridID1) == GRID_GAUSSIAN || gridInqType(gridID1) == GRID_LONLAT )
	{
	  int xsize = gridInqXsize(gridID1);
	  int ysize = gridInqYsize(gridID1);
		
	  if ( ysize == gridInqYsize(gridID2) )
	    compare_lat_reg2d(ysize, gridID1, gridID2);
	  else
	    cdoWarning("ysize of input grids differ!");
		
	  if ( xsize == gridInqXsize(gridID2) )
	    compare_lon_reg2d(xsize, gridID1, gridID2);
	  else
	    cdoWarning("xsize of input grids differ!");
	}
      else if ( gridInqType(gridID1) == GRID_CURVILINEAR || gridInqType(gridID1) == GRID_UNSTRUCTURED )
	{
	  compare_grid_unstructured(gridID1, gridID2);
	}
    }
  else if ( gridInqSize(gridID1) > 1 )
    {
      cdoWarning("Grids have different types! First grid: %s; second grid: %s",
		 gridNamePtr(gridInqType(gridID1)), gridNamePtr(gridInqType(gridID2)));
    }
}

static
int cmpnames(const void *s1, const void *s2)
{
  const char *name1 = (const char *)s1;
  const char *name2 = (const char *)s2;

  return strcmp(name1, name2);
}


void vlistCompare(int vlistID1, int vlistID2, int flag)
{
  int varID;
  bool lchecknames = false;

  int nvars = vlistNvars(vlistID1);

  if ( nvars != vlistNvars(vlistID2) )
    cdoAbort("Input streams have different number of variables per timestep!");

  if ( vlistNrecs(vlistID1) != vlistNrecs(vlistID2) )
    cdoAbort("Input streams have different number of %s per timestep!", nvars==1 ? "layers" : "records");

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( nvars > 1 )
	{
	  if ( flag & CMP_NAME )
	    {
	      char name1[CDI_MAX_NAME], name2[CDI_MAX_NAME];
	      vlistInqVarName(vlistID1, varID, name1);
	      vlistInqVarName(vlistID2, varID, name2);
	      strtolower(name1);
	      strtolower(name2);
	      if ( strcmp(name1, name2) != 0 )
		{
		  cdoWarning("Input streams have different parameter names!");
		  lchecknames = true;
		  flag -= CMP_NAME;
		  //    break;
		}
	    }
	}

      if ( flag & CMP_GRIDSIZE )
	{
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	    cdoAbort("Grid size of the input parameters do not match!");
	}
      
      if ( flag & CMP_NLEVEL )
	{
          int zaxisID1 = vlistInqVarZaxis(vlistID1, varID);
          int zaxisID2 = vlistInqVarZaxis(vlistID2, varID);
          if ( zaxisID1 != zaxisID2 )
            {
              int nlev1 = zaxisInqSize(zaxisID1);
              int nlev2 = zaxisInqSize(zaxisID2);
              if ( nlev1 != nlev2 )
                cdoAbort("Number of levels of the input parameters do not match!");

              double *lev1 = (double*) Malloc(nlev1*sizeof(double));
              double *lev2 = (double*) Malloc(nlev1*sizeof(double));
              cdoZaxisInqLevels(zaxisID1, lev1);
              cdoZaxisInqLevels(zaxisID2, lev2);

              bool ldiffer = false;
              for ( int i = 0; i < nlev1; ++i )
                if ( IS_NOT_EQUAL(lev1[i], lev2[i]) )
                  { ldiffer = true; break; }
              if ( ldiffer )
                {
                  ldiffer = false;
                  for ( int i = 0; i < nlev1; ++i )
                    if ( IS_NOT_EQUAL(lev1[i], lev2[nlev1-1-i]) )
                      { ldiffer = true; break; }

                  if ( ldiffer )
                    cdoWarning("Input parameters have different levels!");
                  else
                    cdoWarning("Z-axis orientation differ!");
                  break;
                }
              
              Free(lev1);
              Free(lev2);
            }
        }
    }

  if ( flag & CMP_GRID )
    {
      int gridID1 = vlistInqVarGrid(vlistID1, 0);
      int gridID2 = vlistInqVarGrid(vlistID2, 0);
      cdoCompareGrids(gridID1, gridID2);
    }

  if ( lchecknames )
    {
      NEW_2D(char, names1, nvars, CDI_MAX_NAME);
      NEW_2D(char, names2, nvars, CDI_MAX_NAME);

      for ( varID = 0; varID < nvars; varID++ )
	vlistInqVarName(vlistID1, varID, names1[varID]);
      for ( varID = 0; varID < nvars; varID++ )
	vlistInqVarName(vlistID2, varID, names2[varID]);

      qsort(names1[0], nvars, CDI_MAX_NAME, cmpnames);
      qsort(names2[0], nvars, CDI_MAX_NAME, cmpnames);

      for ( varID = 0; varID < nvars; varID++ )
	if ( strcmp(names1[varID], names2[varID]) != 0 ) break;

      if ( varID == nvars )
	cdoPrint("Use CDO option --sortname to sort the parameter by name (NetCDF only)!");

      DELETE_2D(names1);
      DELETE_2D(names2);
    }
}


int vlistCompareX(int vlistID1, int vlistID2, int flag)
{
  int nvars = vlistNvars(vlistID1);
  int nvars2 = vlistNvars(vlistID2);
  int nlevels2 = zaxisInqSize(vlistInqVarZaxis(vlistID2, 0));

  if ( nvars2 != 1 )
    cdoAbort("Internal problem, vlistCompareX() called with unexpected vlistID2 argument!");

  for ( int varID = 0; varID < nvars; varID++ )
    {
      if ( flag & CMP_GRIDSIZE )
	{
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, 0)) )
	    cdoAbort("Grid size of the input parameters do not match!");
	}
      
      if ( flag & CMP_NLEVEL )
	{
	  if ( (zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
                nlevels2) && nlevels2 > 1 )
	    cdoAbort("Number of levels of the input parameters do not match!");
	}
    }

  if ( flag & CMP_GRID )
    {
      int gridID1 = vlistInqVarGrid(vlistID1, 0);
      int gridID2 = vlistInqVarGrid(vlistID2, 0);
      cdoCompareGrids(gridID1, gridID2);
    }

  return nlevels2;
}


int vlistIsSzipped(int vlistID)
{
  int lszip = FALSE;
  int nvars = vlistNvars(vlistID);

  for ( int varID = 0; varID < nvars; varID++ )
    {						
      int comptype = vlistInqVarCompType(vlistID, varID);
      if ( comptype == CDI_COMPRESS_SZIP )
	{
	  lszip = TRUE;
	  break;
	}
    }      

  return lszip;
}


int vlistInqNWPV(int vlistID, int varID)
{
  int nwpv; // number of words per value; real:1  complex:2

  if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_CPX32 || 
       vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_CPX64 )
    nwpv = 2;
  else
    nwpv = 1;

  return nwpv;
}


int vlist_check_gridsize(int vlistID)
{
  bool lerror = false;
  int ngrids = vlistNgrids(vlistID);
  int gridID = vlistGrid(vlistID, 0);
  int ngp    = gridInqSize(gridID);

  /* check gridsize */
  for ( int index = 0; index < ngrids; ++index )
    {
      gridID = vlistGrid(vlistID, index);
      if ( ngp != gridInqSize(gridID) )
	{
	  lerror = true;
	  break;
	}
    }

  if ( lerror )
    {
      cdoPrint("This operator requires all variables on the same horizontal grid.");
      cdoPrint("Horizontal grids found:");
      for ( int index = 0; index < ngrids; ++index )
	{
	  gridID = vlistGrid(vlistID, index);
	  cdoPrint("  grid=%d  type=%s  points=%d", index+1, gridNamePtr(gridInqType(gridID)), gridInqSize(gridID));
	}
      cdoAbort("The input stream contains variables on different horizontal grids!");
    }

  return ngp;
}


double *vlist_read_vct(int vlistID, int *rzaxisIDh, int *rnvct, int *rnhlev, int *rnhlevf, int *rnhlevh)
{
  double *vct = NULL;
  int zaxisIDh = -1;
  int nhlev = 0, nhlevf = 0, nhlevh = 0;
  int nvct = 0;
  
  bool lhavevct = false;
  int nzaxis = vlistNzaxis(vlistID);
  for ( int iz = 0; iz < nzaxis; ++iz )
    {
      // bool mono_level = false;
      bool mono_level = true;
      int zaxisID = vlistZaxis(vlistID, iz);
      int nlevel  = zaxisInqSize(zaxisID);
      int zaxistype = zaxisInqType(zaxisID);

      if ( cdoVerbose )
        cdoPrint("ZAXIS_HYBRID = %d ZAXIS_HYBRID_HALF=%d nlevel=%d mono_level=%d",
                 zaxisInqType(zaxisID) == ZAXIS_HYBRID,  zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF, nlevel, mono_level);

      if ( (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && nlevel > 1 && !mono_level )
	{
	  int l;
	  double *level = (double*) Malloc(nlevel*sizeof(double));
	  cdoZaxisInqLevels(zaxisID, level);
	  for ( l = 0; l < nlevel; l++ )
	    {
	      if ( (l+1) != (int) (level[l]+0.5) ) break;
	    }
	  if ( l == nlevel ) mono_level = true; 
	  Free(level);
	}

      if ( (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && nlevel > 1 && mono_level )
	{
	  nvct = zaxisInqVctSize(zaxisID);
	  if ( nlevel == (nvct/2 - 1) )
	    {
	      if ( lhavevct == false )
		{
		  lhavevct = true;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev;
		  nhlevh   = nhlevf + 1;
	      
		  vct = (double*) Malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);
                  if ( cdoVerbose )
                    cdoPrint("Detected half-level model definition : nlevel == (nvct/2 - 1) (nlevel: %d, nvct: %d, nhlevf: %d, nhlevh: %d) ",
                             nlevel, nvct,nhlevf, nhlevh);  
		}
	    }
	  else if ( nlevel == (nvct/2) )
	    {
	      if ( lhavevct == false )
		{
		  lhavevct = true;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev - 1;
		  nhlevh   = nhlev;
	      
		  vct = (double*) Malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);
                  if ( cdoVerbose )
                    cdoPrint("Detected full-level model definition : nlevel == (nvct/2) (nlevel: %d, nvct: %d, nhlevf: %d, nhlevh: %d) ",
                             nlevel, nvct,nhlevf, nhlevh);  
		}
	    }
	  else if ( nlevel == (nvct - 4 - 1) )
	    {
	      if ( lhavevct == false )
		{
		  int vctsize;
		  int voff = 4;

		  double *rvct = (double*) Malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID,rvct);

		  if ( (int)(rvct[0]+0.5) == 100000 && rvct[voff] < rvct[voff+1] )
		    {
		      lhavevct = true;
		      zaxisIDh = zaxisID;
		      nhlev    = nlevel;
		      nhlevf   = nhlev;
		      nhlevh   = nhlev + 1;

		      vctsize = 2*nhlevh;
		      vct = (double*) Malloc(vctsize*sizeof(double));

		      /* calculate VCT for LM */

		      for ( int i = 0; i < vctsize/2; i++ )
			{
			  if ( rvct[voff+i] >= rvct[voff] && rvct[voff+i] <= rvct[3] )
			    {
			      vct[i] = rvct[0]*rvct[voff+i];
			      vct[vctsize/2+i] = 0;
			    }
			  else
			    {
			      vct[i] = (rvct[0]*rvct[3]*(1-rvct[voff+i]))/(1-rvct[3]);
			      vct[vctsize/2+i] = (rvct[voff+i]-rvct[3])/(1-rvct[3]);
			    }
			}
		      
		      if ( cdoVerbose )
			{
			  for ( int i = 0; i < vctsize/2; i++ )
			    fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
			}
		    }
		  Free(rvct);
		}
	    }
	}
    }

  *rzaxisIDh = zaxisIDh;
  *rnvct   = nvct;
  *rnhlev  = nhlev;
  *rnhlevf = nhlevf;
  *rnhlevh = nhlevh;
  
  return vct;
}


void vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2)
{
  int nvct0 = 0;
  double *vct = NULL;

  int nzaxis = vlistNzaxis(vlistID1);
  for ( int i = 0; i < nzaxis; ++i )
    {
      int zaxisID = vlistZaxis(vlistID1, i);
      int nlevel  = zaxisInqSize(zaxisID);

      if ( zaxisID == zaxisID1 && nlevel > 1 )
	{
	  int nvct = zaxisInqVctSize(zaxisID);
          if ( vct == NULL )
            {
              nvct0 = nvct;
              vct = (double*) Malloc(nvct*sizeof(double));
              zaxisInqVct(zaxisID, vct);

              vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
            }
          else
            {
              if ( nvct0 == nvct && memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
                vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
	    }
	}
    }

  if ( vct ) Free(vct);
}


int vlist_get_psvarid(int vlistID, int zaxisID)
{
  int psvarid = -1;
  char name[CDI_MAX_NAME];
  char psname[CDI_MAX_NAME]; psname[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_PSNAME, CDI_MAX_NAME, psname);

  if ( psname[0] )
    {
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; ++varID )
        {
          vlistInqVarName(vlistID, varID, name);
          if ( strcmp(name, psname) == 0 )
            {
              psvarid = varID;
              break;
            }
        }
      if ( cdoVerbose && psvarid == -1 )
        cdoWarning("Surface pressure variable not found - %s", psname);
    }

  return psvarid;
}
