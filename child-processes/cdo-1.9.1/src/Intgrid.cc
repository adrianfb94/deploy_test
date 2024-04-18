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

      Intgrid    interpolate     PINGO grid interpolation
      Intgrid    intgridbil      Bilinear grid interpolation
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"
#include "grid.h"


int genThinoutGrid(int gridID1, int xinc, int yinc)
{
  int gridtype = gridInqType(gridID1);
  int nlon1 = gridInqXsize(gridID1);
  int nlat1 = gridInqYsize(gridID1);

  int nlon2 = nlon1/xinc;
  int nlat2 = nlat1/yinc;
  if ( nlon1%xinc ) nlon2++;
  if ( nlat1%yinc ) nlat2++;
  int gridsize2 = nlon2*nlat2;

  int gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      double *xvals1 = (double*) Malloc(nlon1*sizeof(double));
      double *yvals1 = (double*) Malloc(nlat1*sizeof(double));
      double *xvals2 = (double*) Malloc(nlon2*sizeof(double));
      double *yvals2 = (double*) Malloc(nlat2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      int olat = 0;
      for ( int ilat = 0; ilat < nlat1; ilat+=yinc )
	{
	  yvals2[olat] = yvals1[ilat];
	  olat++;
	}

      int olon = 0;
      for ( int ilon = 0; ilon < nlon1; ilon+=xinc )
	{
	  xvals2[olon] = xvals1[ilon];
	  olon++;
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      Free(xvals1);
      Free(yvals1);
      Free(xvals2);
      Free(yvals2);
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}


int genBoxavgGrid(int gridID1, int xinc, int yinc)
{
  int i, j, i1;

  int gridtype = gridInqType(gridID1);
  int nlon1 = gridInqXsize(gridID1);
  int nlat1 = gridInqYsize(gridID1);

  int nlon2 = nlon1/xinc;
  int nlat2 = nlat1/yinc;
  if ( nlon1%xinc ) nlon2++;
  if ( nlat1%yinc ) nlat2++;
  int gridsize2 = nlon2*nlat2;

  int gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      double *grid1_corner_lon = NULL, *grid1_corner_lat = NULL;
      double *grid2_corner_lon = NULL, *grid2_corner_lat = NULL;
      double *xvals1 = (double*) Malloc(nlon1*sizeof(double));
      double *yvals1 = (double*) Malloc(nlat1*sizeof(double));
      double *xvals2 = (double*) Malloc(nlon2*sizeof(double));
      double *yvals2 = (double*) Malloc(nlat2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridInqYbounds(gridID1, NULL) && gridInqXbounds(gridID1, NULL) )
	{
	  grid1_corner_lon = (double*) Malloc(2*nlon1*sizeof(double));
	  grid1_corner_lat = (double*) Malloc(2*nlat1*sizeof(double));
	  grid2_corner_lon = (double*) Malloc(2*nlon2*sizeof(double));
	  grid2_corner_lat = (double*) Malloc(2*nlat2*sizeof(double));
	  gridInqXbounds(gridID1, grid1_corner_lon);
	  gridInqYbounds(gridID1, grid1_corner_lat);
	}

      j = 0;
      for ( i = 0; i < nlon1; i += xinc )
	{
	  i1 = i+(xinc-1);
	  if ( i1 >= nlon1-1 ) i1 = nlon1-1; 
	  xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i])/2;
	  if ( grid2_corner_lon )
	    {
	      grid2_corner_lon[2*j] = grid1_corner_lon[2*i];
	      grid2_corner_lon[2*j+1] = grid1_corner_lon[2*i1+1];
	    }
	  j++;
	}
      j = 0;
      for ( i = 0; i < nlat1; i += yinc )
	{
	  i1 = i+(yinc-1);
	  if ( i1 >= nlat1-1 ) i1 = nlat1-1; 
	  yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i])/2;
	  if ( grid2_corner_lat )
	    {
	      grid2_corner_lat[2*j] = grid1_corner_lat[2*i];
	      grid2_corner_lat[2*j+1] = grid1_corner_lat[2*i1+1];
	    }
	  j++;
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      Free(xvals1);
      Free(yvals1);
      Free(xvals2);
      Free(yvals2);

      if ( grid2_corner_lon && grid2_corner_lat )
	{
	  gridDefNvertex(gridID2, 2);
	  gridDefXbounds(gridID2, grid2_corner_lon);
	  gridDefYbounds(gridID2, grid2_corner_lat);

	  Free(grid2_corner_lon);
	  Free(grid2_corner_lat);
	}
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}

static
void boxavg(field_type *field1, field_type *field2, int xinc, int yinc)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1  = field1->ptr;
  double *array2  = field2->ptr;
  double missval = field1->missval;

  int nlon1 = gridInqXsize(gridID1);
  int nlat1 = gridInqYsize(gridID1);

  int nlon2 = gridInqXsize(gridID2);
  int nlat2 = gridInqYsize(gridID2);

  double **xfield1 = (double **) Malloc(nlat1*sizeof(double *));

  for ( int ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;

  double **xfield2 = (double **) Malloc(nlat2 * sizeof(double *));

  for ( int ilat = 0; ilat < nlat2; ilat++ )
    xfield2[ilat] = array2 + ilat*nlon2;

  for ( int ilat = 0; ilat < nlat2; ilat++ )
    for ( int ilon = 0; ilon < nlon2; ilon++ )
      {
	xfield2[ilat][ilon] = 0;

	int in = 0;
	for ( int j = 0; j < yinc; ++j )
	  {
	    int jj = ilat*yinc+j;
	    if ( jj >= nlat1 ) break;
	    for ( int i = 0; i < xinc; ++i )
	      {
		int ii = ilon*xinc+i;
		if ( ii >= nlon1 ) break;
		in++;
		xfield2[ilat][ilon] += xfield1[jj][ii];
	      }
	  }
	xfield2[ilat][ilon] /= in;
      }

  int nmiss = 0;
  for ( int i = 0; i < nlat2*nlon2; i++ )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

  field2->nmiss = nmiss;

  Free(xfield2);
  Free(xfield1);
}

static
void thinout(field_type *field1, field_type *field2, int xinc, int yinc)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1  = field1->ptr;
  double *array2  = field2->ptr;
  double missval = field1->missval;

  int nlon1 = gridInqXsize(gridID1);
  int nlat1 = gridInqYsize(gridID1);

  int nlon2 = gridInqXsize(gridID2);
  int nlat2 = gridInqYsize(gridID2);

  double **xfield1 = (double **) Malloc(nlat1*sizeof(double *));

  for ( int ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;

  double **xfield2 = (double **) Malloc(nlat2*sizeof(double *));

  for ( int ilat = 0; ilat < nlat2; ilat++ )
    xfield2[ilat] = array2 + ilat*nlon2;

  int olat = 0;
  for ( int ilat = 0; ilat < nlat1; ilat+=yinc )
    {
      int olon = 0;
      for ( int ilon = 0; ilon < nlon1; ilon+=xinc )
	{
	  xfield2[olat][olon] = xfield1[ilat][ilon];
	  olon++;
	}
      olat++;
    }

  int nmiss = 0;
  for ( int i = 0; i < nlat2*nlon2; i++ )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
  
  field2->nmiss = nmiss;

  Free(xfield2);
  Free(xfield1);
}



void *Intgrid(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridID1 = -1, gridID2 = -1;
  int nmiss;
  int xinc = 0, yinc = 0;
  double missval;

  cdoInitialize(argument);

  // clang-format off
  int INTGRIDBIL  = cdoOperatorAdd("intgridbil",  0, 0, NULL);
  int INTPOINT    = cdoOperatorAdd("intpoint",    0, 0, NULL);
  int INTERPOLATE = cdoOperatorAdd("interpolate", 0, 0, NULL);
  int BOXAVG      = cdoOperatorAdd("boxavg",      0, 0, NULL);
  int THINOUT     = cdoOperatorAdd("thinout",     0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  // open stream before calling cdoDefineGrid!!!
  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  if ( operatorID == INTGRIDBIL || operatorID == INTERPOLATE )
    {
      operatorInputArg("grid description file or name");
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == INTPOINT )
    {
      operatorInputArg("longitude and latitude");
      operatorCheckArgc(2);
      double slon = parameter2double(operatorArgv()[0]);
      double slat = parameter2double(operatorArgv()[1]);
      gridID2 = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
      gridDefXvals(gridID2, &slon);
      gridDefYvals(gridID2, &slat);
    }
  else if ( operatorID == THINOUT || operatorID == BOXAVG )
    {
      operatorInputArg("xinc, yinc");
      operatorCheckArgc(2);
      xinc = parameter2int(operatorArgv()[0]);
      yinc = parameter2int(operatorArgv()[1]);
    }

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  for ( int index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( operatorID == BOXAVG || operatorID == THINOUT )
	{
	  if ( index == 0 )
	    {
	      if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN 
		   /* && gridInqType(gridID1) != GRID_CURVILINEAR */ )
		cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

	      if ( operatorID == BOXAVG )
		gridID2 = genBoxavgGrid(gridID1, xinc, yinc);
	      else
		gridID2 = genThinoutGrid(gridID1, xinc, yinc);
	    }
	  else
	    cdoAbort("Too many different grids!");
	}
      else
	{
          bool ldistgen = false;
          if ( grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2) ) ldistgen = true;
          
	  if ( !ldistgen && gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN )
	    cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );
	}

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  field1.grid    = gridID1;
	  field1.nmiss   = nmiss;
	  field1.missval = missval;
	  field1.ptr     = array1;
	  field2.grid    = gridID2;
	  field2.ptr     = array2;
	  field2.nmiss   = 0;

	  if ( operatorID == INTGRIDBIL || operatorID == INTPOINT )
	    intgridbil(&field1, &field2);
	  else if ( operatorID == INTERPOLATE )
	    interpolate(&field1, &field2);
	  else if ( operatorID == BOXAVG )
	    boxavg(&field1, &field2, xinc, yinc);
	  else if ( operatorID == THINOUT )
	    thinout(&field1, &field2, xinc, yinc);

	  nmiss = field2.nmiss;

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  cdoFinish();

  return 0;
}
