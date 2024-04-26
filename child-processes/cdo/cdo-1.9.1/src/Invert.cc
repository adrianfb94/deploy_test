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

      Invert     invertlat       Invert latitude
      Invert     invertlon       Invert longitude
      Invert     invertlatdes    Invert latitude description
      Invert     invertlondes    Invert longitude description
      Invert     invertlatdata   Invert latitude data
      Invert     invertlondata   Invert longitude data
*/

#include <cdi.h>
#include "cdo_int.h"
#include "pstream.h"


static
void invertLonDes(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for ( int index = 0; index < ngrids; index++ )
    {
      int gridID1 = vlistGrid(vlistID, index);
      int gridID2 = gridDuplicate(gridID1);

      int gridtype = gridInqType(gridID1);

      if ( !(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION ||
             gridtype == GRID_LONLAT  || gridtype == GRID_CURVILINEAR) )
	cdoAbort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      if ( gridInqXvals(gridID1, NULL) )
	{
	  int nlon = gridInqXsize(gridID1);
	  int nlat = gridInqYsize(gridID1);
	  int size = (gridtype == GRID_CURVILINEAR) ? nlon*nlat : nlon;

	  double *xv1 = (double*) Malloc(size*sizeof(double));
	  double *xv2 = (double*) Malloc(size*sizeof(double));

	  gridInqXvals(gridID1, xv1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      for ( int ilat = 0; ilat < nlat; ilat++ )
		for ( int ilon = 0; ilon < nlon; ilon++ )
		  xv2[ilat*nlon + nlon-ilon-1] = xv1[ilat*nlon + ilon];
	    }
	  else
	    {
	      for ( int ilon = 0; ilon < nlon; ilon++ )
		xv2[nlon-ilon-1] = xv1[ilon];
	    }

	  gridDefXvals(gridID2, xv2);

	  if ( xv2 ) Free(xv2);
	  if ( xv1 ) Free(xv1);
	}

      if ( gridInqXbounds(gridID1, NULL) )
	{
	  int nlon = gridInqXsize(gridID1);
	  int nlat = gridInqYsize(gridID1);
	  int nv   = gridInqNvertex(gridID1);
	  int size = (gridtype == GRID_CURVILINEAR) ? nv*nlon*nlat : nv*nlon;

	  double *xb1 = (double*) Malloc(size*sizeof(double));
	  double *xb2 = (double*) Malloc(size*sizeof(double));

	  gridInqXbounds(gridID1, xb1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      for ( int ilat = 0; ilat < nlat; ilat++ )
		for ( int ilon = 0; ilon < nlon; ilon++ )
		  for ( int iv = 0; iv < nv; iv++ )
		    xb2[ilat*nlon*nv + (nlon-ilon-1)*nv + iv] = xb1[ilat*nlon*nv + ilon*nv + iv];
	    }
	  else
	    {
		for ( int ilon = 0; ilon < nlon; ilon++ )
		  {
		    xb2[nlon*2-ilon*2-1] = xb1[ilon*2];
		    xb2[nlon*2-ilon*2-2] = xb1[ilon*2+1];
		  }
	    }

	  gridDefXbounds(gridID2, xb2);

	  if ( xb2 ) Free(xb2);
	  if ( xb1 ) Free(xb1);
	}

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static
void invertLatCoord(int gridID)
{
  int gridtype = gridInqType(gridID);

  if ( gridInqYvals(gridID, NULL) )
    {
      int nlon = gridInqXsize(gridID);
      int nlat = gridInqYsize(gridID);
      int size = (gridtype == GRID_CURVILINEAR) ? nlon*nlat : nlat;

      double *yv1 = (double*) Malloc(size*sizeof(double));
      double *yv2 = (double*) Malloc(size*sizeof(double));

      if ( gridtype == GRID_CURVILINEAR )
        {
          gridInqXvals(gridID, yv1);
          
          for ( int ilat = 0; ilat < nlat; ilat++ )
            for ( int ilon = 0; ilon < nlon; ilon++ )
              yv2[(nlat-ilat-1)*nlon + ilon] = yv1[ilat*nlon + ilon];

          gridDefXvals(gridID, yv2);
          
          gridInqYvals(gridID, yv1);

          for ( int ilat = 0; ilat < nlat; ilat++ )
            for ( int ilon = 0; ilon < nlon; ilon++ )
              yv2[(nlat-ilat-1)*nlon + ilon] = yv1[ilat*nlon + ilon];

          gridDefYvals(gridID, yv2);
        }
      else
        {
          gridInqYvals(gridID, yv1);

          for ( int ilat = 0; ilat < nlat; ilat++ )
            yv2[nlat-ilat-1] = yv1[ilat];

          gridDefYvals(gridID, yv2);
        }

      if ( yv2 ) Free(yv2);
      if ( yv1 ) Free(yv1);
    }

  if ( gridInqYbounds(gridID, NULL) )
    {
      int nlon = gridInqXsize(gridID);
      int nlat = gridInqYsize(gridID);
      int nv   = gridInqNvertex(gridID);
      int size = (gridtype == GRID_CURVILINEAR) ? nv*nlon*nlat : nv*nlat;

      double *yb1 = (double*) Malloc(size*sizeof(double));
      double *yb2 = (double*) Malloc(size*sizeof(double));

      gridInqYbounds(gridID, yb1);

      if ( gridtype == GRID_CURVILINEAR )
        {
          for ( int ilat = 0; ilat < nlat; ilat++ )
            for ( int ilon = 0; ilon < nlon; ilon++ )
              for ( int iv = 0; iv < nv; iv++ )
                yb2[(nlat-ilat-1)*nlon*nv + ilon*nv + iv] = yb1[ilat*nlon*nv + ilon*nv + iv];
        }
      else
        {
          for ( int ilat = 0; ilat < nlat; ilat++ )
            {
              yb2[nlat*2-ilat*2-1] = yb1[ilat*2];
              yb2[nlat*2-ilat*2-2] = yb1[ilat*2+1];
            }
        }

      gridDefYbounds(gridID, yb2);

      if ( yb2 ) Free(yb2);
      if ( yb1 ) Free(yb1);
    }
}

static
void invertLatDes(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for ( int index = 0; index < ngrids; index++ )
    {
      int gridID1 = vlistGrid(vlistID, index);
      int gridID2 = gridDuplicate(gridID1);

      int gridtype = gridInqType(gridID1);

      if ( !(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION ||
             gridtype == GRID_LONLAT  || gridtype == GRID_CURVILINEAR) )
	cdoAbort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      invertLatCoord(gridID2);

      int projID = gridInqProj(gridID2);
      if ( projID != CDI_UNDEFID ) invertLatCoord(projID);

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static
void invertLonData(double *array1, double *array2, int gridID1)
{
  int nlon = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  if ( nlat > 0 )
    {
      double **field1 = (double **) Malloc(nlat*sizeof(double *));
      double **field2 = (double **) Malloc(nlat*sizeof(double *));
  
      for ( int ilat = 0; ilat < nlat; ilat++ )
	{
	  field1[ilat] = array1 + ilat*nlon;
	  field2[ilat] = array2 + ilat*nlon;
	}

      for ( int ilat = 0; ilat < nlat; ilat++ )
	for ( int ilon = 0; ilon < nlon; ilon++ )
	  field2[ilat][nlon-ilon-1] = field1[ilat][ilon];
  
      if ( field1 ) Free(field1);
      if ( field2 ) Free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}

static
void invertLatData(double *array1, double *array2, int gridID1)
{
  int nlon = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  if ( nlat > 0 )
    {
      double **field1 = (double **) Malloc(nlat*sizeof(double *));
      double **field2 = (double **) Malloc(nlat*sizeof(double *));
  
      for ( int ilat = 0; ilat < nlat; ilat++ )
	{
	  field1[ilat] = array1 + ilat*nlon;
	  field2[ilat] = array2 + ilat*nlon;
	}

      for ( int ilat = 0; ilat < nlat; ilat++ )
	memcpy(field2[nlat-ilat-1], field1[ilat], nlon*sizeof(double));
      
      if ( field1 ) Free(field1);
      if ( field2 ) Free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}


void *Invert(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridID1;
  int nmiss;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("invertlat",     func_all, func_lat, NULL);
  cdoOperatorAdd("invertlon",     func_all, func_lon, NULL);
  cdoOperatorAdd("invertlatdes",  func_hrd, func_lat, NULL);
  cdoOperatorAdd("invertlondes",  func_hrd, func_lon, NULL);
  cdoOperatorAdd("invertlatdata", func_fld, func_lat, NULL);
  cdoOperatorAdd("invertlondata", func_fld, func_lon, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc1 = cdoOperatorF1(operatorID);
  int operfunc2 = cdoOperatorF2(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc1 == func_all || operfunc1 == func_hrd )
    {
      if ( operfunc2 == func_lat ) invertLatDes(vlistID2);
      else                         invertLonDes(vlistID2);
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss);

	  pstreamDefRecord(streamID2, varID, levelID);

	  if ( operfunc1 == func_all || operfunc1 == func_fld )
	    {
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      if ( operfunc2 == func_lat )
		invertLatData(array1, array2, gridID1);
	      else
		invertLonData(array1, array2, gridID1);

	      pstreamWriteRecord(streamID2, array2, nmiss);     
	    }
	  else
	    {
	      pstreamWriteRecord(streamID2, array1, nmiss);     
	    }
	}
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);

  cdoFinish();

  return 0;
}
