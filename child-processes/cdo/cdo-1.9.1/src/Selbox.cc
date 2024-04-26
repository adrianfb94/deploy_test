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

      Selbox     sellonlatbox    Select lon/lat box
      Selbox     selindexbox     Select index box
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


static
void correct_xvals(long nlon, long inc, double *xvals)
{
  if ( IS_EQUAL(xvals[0], xvals[(nlon-1)*inc]) ) xvals[(nlon-1)*inc] += 360;

  if ( xvals[0] > xvals[(nlon-1)*inc] )
    for ( int i = 0; i < nlon; i++ )
      if ( xvals[i*inc] >= 180 ) xvals[i*inc] -= 360;

  for ( int i = 0; i < nlon; i++ )
    {
      if ( xvals[i*inc] < -180 ) xvals[i*inc] += 360;
      if ( xvals[i*inc] >  360 ) xvals[i*inc] -= 360;
    }

  if ( xvals[0] > xvals[(nlon-1)*inc] )
    for ( int i = 1; i < nlon; i++ )
      if ( xvals[i*inc] < xvals[(i-1)*inc] ) xvals[i*inc] += 360;
}

static
void gengridxyvals(int gridtype, int gridID1, int gridID2, int nlon, int nlat, int nlon2, int nlat2,
                   int lat1, int lat2, int lon11, int lon12, int lon21, int lon22, const char *xunits)
{
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;

  int lxvals = gridInqXvals(gridID1, NULL);
  int lyvals = gridInqYvals(gridID1, NULL);

  if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lxvals && lyvals )
	{
	  xvals1 = (double*) Malloc(nlon*nlat*sizeof(double));
	  yvals1 = (double*) Malloc(nlon*nlat*sizeof(double));
	  xvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));
	}
    }
  else
    {
      if ( lxvals ) xvals1 = (double*) Malloc(nlon*sizeof(double));
      if ( lyvals ) yvals1 = (double*) Malloc(nlat*sizeof(double));
      if ( lxvals ) xvals2 = (double*) Malloc(nlon2*sizeof(double));
      if ( lyvals ) yvals2 = (double*) Malloc(nlat2*sizeof(double));
    }

  double *pxvals2 = xvals2;
  double *pyvals2 = yvals2;

  if ( xvals1 ) gridInqXvals(gridID1, xvals1);
  if ( yvals1 ) gridInqYvals(gridID1, yvals1);

  if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lxvals && lyvals )
	for ( int ilat = lat1; ilat <= lat2; ilat++ )
	  {
	    for ( int ilon = lon21; ilon <= lon22; ilon++ )
	      {
		*pxvals2++ = xvals1[ilat*nlon + ilon];
		*pyvals2++ = yvals1[ilat*nlon + ilon];
	      }
	    for ( int ilon = lon11; ilon <= lon12; ilon++ )
	      {
		*pxvals2++ = xvals1[ilat*nlon + ilon];
		*pyvals2++ = yvals1[ilat*nlon + ilon];
	      }
	  }
    }
  else
    {
      if ( lxvals )
	{
	  for ( int i = lon21; i <= lon22; i++ ) *pxvals2++ = xvals1[i];
	  for ( int i = lon11; i <= lon12; i++ ) *pxvals2++ = xvals1[i];
	  if ( xunits && strncmp(xunits, "degree", 6) == 0 ) correct_xvals(nlon2, 1, xvals2);
	}
      
      if ( lyvals ) for ( int i = lat1;  i <= lat2;  i++ ) *pyvals2++ = yvals1[i];
    }
  /*
    for ( int i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
    for ( int i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
  */
  if ( xvals2 ) gridDefXvals(gridID2, xvals2);
  if ( yvals2 ) gridDefYvals(gridID2, yvals2);

  if ( xvals1 ) Free(xvals1);
  if ( yvals1 ) Free(yvals1);
  if ( xvals2 ) Free(xvals2);
  if ( yvals2 ) Free(yvals2);
}

static
int gengrid(int gridID1, int lat1, int lat2, int lon11, int lon12, int lon21, int lon22)
{
  int nlon = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  int nlon21 = lon12 - lon11 + 1;
  int nlon22 = lon22 - lon21 + 1;
  int nlon2 = nlon21 + nlon22;
  int nlat2 = lat2 - lat1 + 1;

  if ( cdoVerbose ) cdoPrint("nlon=%d nlat=%d", nlon2, nlat2);

  int gridtype = gridInqType(gridID1);

  int gridID2 = gridCreate(gridtype, nlon2*nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefNP(gridID2, gridInqNP(gridID1));
  gridDefDatatype(gridID2, gridInqDatatype(gridID1));

  grid_copy_attributes(gridID1, gridID2);

  if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridID1, gridID2);

  char xunits[CDI_MAX_NAME]; xunits[0] = 0;
  char yunits[CDI_MAX_NAME]; yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  gengridxyvals(gridtype, gridID1, gridID2, nlon, nlat, nlon2, nlat2,
                lat1, lat2, lon11, lon12, lon21, lon22, xunits);

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      double *xbounds1 = NULL, *ybounds1 = NULL;
      double *xbounds2 = NULL, *ybounds2 = NULL;

      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = (double*) Malloc(4*nlon*nlat*sizeof(double));
	  ybounds1 = (double*) Malloc(4*nlon*nlat*sizeof(double));
	  xbounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xbounds1 = (double*) Malloc(2*nlon*sizeof(double));
	  ybounds1 = (double*) Malloc(2*nlat*sizeof(double));
	  xbounds2 = (double*) Malloc(2*nlon2*sizeof(double));
	  ybounds2 = (double*) Malloc(2*nlat2*sizeof(double));
	}

      double *pxbounds2 = xbounds2;
      double *pybounds2 = ybounds2;

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  gridDefNvertex(gridID2, 4);
	  for ( int ilat = lat1; ilat <= lat2; ilat++ )
	    {
	      for ( int ilon = 4*lon21; ilon < 4*(lon22+1); ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon + ilon];
		}
	      for ( int ilon = 4*lon11; ilon < 4*(lon12+1); ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon + ilon];
		}
	    }
	}
      else
	{
	  gridDefNvertex(gridID2, 2);
	  for ( int i = 2*lon21; i < 2*(lon22+1); i++ ) *pxbounds2++ = xbounds1[i];
	  for ( int i = 2*lon11; i < 2*(lon12+1); i++ ) *pxbounds2++ = xbounds1[i];
	  for ( int i = 2*lat1;  i < 2*(lat2+1);  i++ ) *pybounds2++ = ybounds1[i];

	  if ( strncmp(xunits, "degree", 6) == 0 )
	    {
	      correct_xvals(nlon2, 2, xbounds2);
	      correct_xvals(nlon2, 2, xbounds2+1);
	    }
	}

      gridDefXbounds(gridID2, xbounds2);
      gridDefYbounds(gridID2, ybounds2);

      Free(xbounds1);
      Free(ybounds1);
      Free(xbounds2);
      Free(ybounds2);
    }

  int projID1 = gridInqProj(gridID1);
  if ( projID1 != CDI_UNDEFID && gridInqType(projID1) == GRID_PROJECTION )
    {
      int projID2 = gridCreate(GRID_PROJECTION, nlon2*nlat2);
      gridDefXsize(projID2, nlon2);
      gridDefYsize(projID2, nlat2);

      grid_copy_attributes(projID1, projID2);
      grid_copy_mapping(projID1, projID2);

      gengridxyvals(GRID_PROJECTION, projID1, projID2, nlon, nlat, nlon2, nlat2,
                    lat1, lat2, lon11, lon12, lon21, lon22, NULL);

      gridDefProj(gridID2, projID2);
    }

  return gridID2;
}


int gengridcell(int gridID1, int gridsize2, int *cellidx)
{
  int gridID2 = -1;
  int gridtype  = gridInqType(gridID1);
  int gridsize1 = gridInqSize(gridID1);
  int prec      = gridInqDatatype(gridID1);

  if ( gridtype == GRID_CURVILINEAR ) gridtype = GRID_UNSTRUCTURED;

  if ( gridtype == GRID_UNSTRUCTURED )
    gridID2 = gridCreate(gridtype, gridsize2);
  else
    return gridID2;

  gridDefDatatype(gridID2, prec);

  grid_copy_attributes(gridID1, gridID2);

  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      double *xvals1 = (double*) Malloc(gridsize1*sizeof(double));
      double *yvals1 = (double*) Malloc(gridsize1*sizeof(double));
      double *xvals2 = (double*) Malloc(gridsize2*sizeof(double));
      double *yvals2 = (double*) Malloc(gridsize2*sizeof(double));

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      for ( int i = 0; i < gridsize2; ++i )
	{
	  xvals2[i] = xvals1[cellidx[i]];
	  yvals2[i] = yvals1[cellidx[i]];
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      Free(xvals1);
      Free(yvals1);
      Free(xvals2);
      Free(yvals2);
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      int nv = gridInqNvertex(gridID1);

      double *xbounds1 = (double*) Malloc(nv*gridsize1*sizeof(double));
      double *ybounds1 = (double*) Malloc(nv*gridsize1*sizeof(double));
      double *xbounds2 = (double*) Malloc(nv*gridsize2*sizeof(double));
      double *ybounds2 = (double*) Malloc(nv*gridsize2*sizeof(double));

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      gridDefNvertex(gridID2, nv);

      for ( int i = 0; i < gridsize2; ++i )
	{
	  for ( int k = 0; k < nv; ++k )
	    {
	      xbounds2[i*nv+k] = xbounds1[cellidx[i]*nv+k];
	      ybounds2[i*nv+k] = ybounds1[cellidx[i]*nv+k];
	    }
	}

      gridDefXbounds(gridID2, xbounds2);
      gridDefYbounds(gridID2, ybounds2);

      Free(xbounds1);
      Free(ybounds1);
      Free(xbounds2);
      Free(ybounds2);
    }

  return gridID2;
}


void genlonlatbox_reg(int gridID, double xlon1, double xlon2, double xlat1, double xlat2,
		      int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  int nlon = gridInqXsize(gridID);
  int nlat = gridInqYsize(gridID);

  double *xvals = (double*) Malloc(nlon*sizeof(double));
  double *yvals = (double*) Malloc(nlat*sizeof(double));

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  gridInqXunits(gridID, xunits);
  gridInqYunits(gridID, yunits);

  double xfact = 1, yfact = 1;
  if ( strncmp(xunits, "radian", 6) == 0 ) xfact = RAD2DEG;
  if ( strncmp(yunits, "radian", 6) == 0 ) yfact = RAD2DEG;

  for ( int ilat = 0; ilat < nlat; ilat++ ) yvals[ilat] *= yfact;
  for ( int ilon = 0; ilon < nlon; ilon++ ) xvals[ilon] *= xfact;

  if ( IS_NOT_EQUAL(xlon1, xlon2) )
    {
      xlon2 -= 360 * floor ((xlon2 - xlon1) / 360);
      if ( IS_EQUAL(xlon1, xlon2) ) xlon2 += 360;
    }
  else
    {
      xlon2 += 0.00001;
    }

  xlon2 -= 360 * floor ((xlon1 - xvals[0]) / 360);
  xlon1 -= 360 * floor ((xlon1 - xvals[0]) / 360);

  // while ( nlon == 1 || (xvals[nlon-1] - xvals[0]) >= 360 ) nlon--;

  for ( *lon21 = 0; *lon21 < nlon && xvals[*lon21] < xlon1; (*lon21)++ );
  for ( *lon22 = *lon21; *lon22 < nlon && xvals[*lon22] < xlon2; (*lon22)++ );

  if ( *lon22 >= nlon || xvals[*lon22] > xlon2 ) (*lon22)--;

  xlon1 -= 360;
  xlon2 -= 360;

  for ( *lon11 = 0; xvals[*lon11] < xlon1; (*lon11)++ );
  for ( *lon12 = *lon11; *lon12 < nlon && xvals[*lon12] < xlon2; (*lon12)++ );
  
  // (*lon12)--;
  if ( *lon12 >= nlon || xvals[*lon12] > xlon2 ) (*lon12)--;
  if ( *lon12 >= 0 && IS_EQUAL(xvals[*lon12], xvals[*lon21]) ) (*lon12)--;

  if ( *lon12 - *lon11 + 1 + *lon22 - *lon21 + 1 <= 0 )
    cdoAbort("Longitudinal dimension is too small!");
  
  if ( yvals[0] > yvals[nlat - 1] )
    {
      if ( xlat1 > xlat2 )
	{
	  for ( *lat1 = 0; *lat1 < nlat && yvals[*lat1] > xlat1; (*lat1)++ );
	  for ( *lat2 = nlat - 1; *lat2 && yvals[*lat2] < xlat2; (*lat2)-- );
	}
      else
	{
	  for ( *lat1 = 0; *lat1 < nlat && yvals[*lat1] > xlat2; (*lat1)++ );
	  for ( *lat2 = nlat - 1; *lat2 && yvals[*lat2] < xlat1; (*lat2)-- );
	}
    }
  else
    {
      if ( xlat1 < xlat2 )
	{
	  for ( *lat1 = 0; *lat1 < nlat && yvals[*lat1] < xlat1; (*lat1)++ );
	  for ( *lat2 = nlat - 1; *lat2 && yvals[*lat2] > xlat2; (*lat2)-- );
	}
      else
	{
	  for ( *lat1 = 0; *lat1 < nlat && yvals[*lat1] < xlat2; (*lat1)++ );
	  for ( *lat2 = nlat - 1; *lat2 && yvals[*lat2] > xlat1; (*lat2)-- );
	}
    }

  Free(xvals);
  Free(yvals);

  if ( *lat2 - *lat1 + 1 <= 0 )
    cdoAbort("Latitudinal dimension is too small!");
}

static
void genlonlatbox_curv(int gridID, double xlon1, double xlon2, double xlat1, double xlat2,
                       int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{  
  int nlon = gridInqXsize(gridID);
  int nlat = gridInqYsize(gridID);
  int gridsize = nlon*nlat;

  int grid_is_circular = gridIsCircular(gridID);

  double *xvals = (double *) Malloc(gridsize*sizeof(double));
  double *yvals = (double *) Malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  gridInqXunits(gridID, xunits);
  gridInqYunits(gridID, yunits);

  double xfact = 1, yfact = 1;
  if ( strncmp(xunits, "radian", 6) == 0 ) xfact = RAD2DEG;
  if ( strncmp(yunits, "radian", 6) == 0 ) yfact = RAD2DEG;

  if ( xlon1 > xlon2 ) 
    cdoAbort("The second longitude have to be greater than the first one!");
  /*
  if ( xlon2 < 180. )
    {
      xlon1 += 360;
      xlon2 += 360;
    }
  */
  if ( xlat1 > xlat2 )
    {
      double xtemp = xlat1;
      xlat1 = xlat2;
      xlat2 = xtemp;
    }

  *lat1 = nlat-1;
  *lat2 = 0;
  *lon11 = 0;
  *lon12 = -1;
  /*
  *lon11 = nlon-1;
  *lon12 = 0;
  */
  *lon21 = nlon-1;
  *lon22 = 0;

  bool lp2 = false;
  double xfirst, xlast, ylast;

  if ( grid_is_circular )
    {
      for ( int ilat = 0; ilat < nlat; ilat++ )
        {
          xlast = xfact * xvals[ilat*nlon + nlon-1];
          ylast = yfact * yvals[ilat*nlon + nlon-1];
          if ( ylast >= xlat1 && ylast <= xlat2 )
            if ( xlon1 <= xlast && xlon2 > xlast && (xlon2-xlon1) < 360 )
              {
                *lon11 = nlon-1;
                *lon12 = 0;
                lp2 = true;
                break;
              }
        }
    }

  for ( int ilat = 0; ilat < nlat; ilat++ )
    for ( int ilon = 0; ilon < nlon; ilon++ )
      {
        double xval = xfact * xvals[ilat*nlon + ilon];
        double yval = yfact * yvals[ilat*nlon + ilon];
        if ( yval >= xlat1 && yval <= xlat2 )
          {
            if ( lp2 )
              {
                xfirst = xfact * xvals[ilat*nlon];
                if ( xfirst < xlon1 ) xfirst = xlon1;
                  
                xlast = xfact * xvals[ilat*nlon + nlon-1];
                if ( xlast > xlon2 ) xlast = xlon2;

                if ( xval >= xlon1 && xval <= xlast )
                  {
                    if ( ilon < *lon21 ) *lon21 = ilon;
                    if ( ilon > *lon22 ) *lon22 = ilon;
                    if ( ilat < *lat1 ) *lat1 = ilat;
                    if ( ilat > *lat2 ) *lat2 = ilat;
                  }
                else if ( xval >= xfirst && xval <= xlon2 )
                  {
                    if ( ilon < *lon11 ) *lon11 = ilon;
                    if ( ilon > *lon12 ) *lon12 = ilon;
                    if ( ilat < *lat1 ) *lat1 = ilat;
                    if ( ilat > *lat2 ) *lat2 = ilat;
                  }
              }
            else
              {
                if ( ((xval     >= xlon1 && xval     <= xlon2) ||
                      (xval-360 >= xlon1 && xval-360 <= xlon2) ||
                      (xval+360 >= xlon1 && xval+360 <= xlon2)) )
                  {
                    if ( ilon < *lon21 ) *lon21 = ilon;
                    if ( ilon > *lon22 ) *lon22 = ilon;
                    if ( ilat < *lat1 ) *lat1 = ilat;
                    if ( ilat > *lat2 ) *lat2 = ilat;
                  }
              }
            /*
            if ( xval >= xlon1 && xval <= xlon2 )
              {
                if ( ilon < *lon21 ) *lon21 = ilon;
                if ( ilon > *lon22 ) *lon22 = ilon;
                if ( ilat < *lat1 ) *lat1 = ilat;
                if ( ilat > *lat2 ) *lat2 = ilat;
              }
            else if ( xval >= xlon1-360 && xval <= xlon2-360 )
              {
                if ( ilon < *lon11 ) *lon11 = ilon;
                if ( ilon > *lon12 ) *lon12 = ilon;
                if ( ilat < *lat1 ) *lat1 = ilat;
                if ( ilat > *lat2 ) *lat2 = ilat;
              }
            */
          }
      }
  /*
  while ( *lon12 >= *lon21 ) (*lon12)--;
  if ( *lon12 <= *lon11 ) { *lon11 = nlon-1; *lon12 = 0; }
  */
  if ( *lon12 == 0 && *lon11 > 0 ) *lon11 = -1;

  if ( *lat2 - *lat1 + 1 <= 0 )
    cdoAbort("Latitudinal dimension is too small!");

  Free(xvals);
  Free(yvals);
}


void getlonlatparams(int argc_offset, double *xlon1, double *xlon2, double *xlat1, double *xlat2)
{
  bool lset = false;

  int nargc = operatorArgc() - argc_offset;
  if ( nargc == 1 )
    {
      const char *gridname = operatorArgv()[argc_offset+0];
      if ( strcmp(gridname, "europe") == 0 )
        {
          *xlon1 = -30;
          *xlon2 =  60;
          *xlat1 =  30;
          *xlat2 =  80;
          lset = true;
        }
    }

  if ( !lset )
    {
      operatorCheckArgc(argc_offset+4);

      *xlon1 = parameter2double(operatorArgv()[argc_offset+0]);
      *xlon2 = parameter2double(operatorArgv()[argc_offset+1]);
      *xlat1 = parameter2double(operatorArgv()[argc_offset+2]);
      *xlat2 = parameter2double(operatorArgv()[argc_offset+3]);
    }
}


void genlonlatbox(int argc_offset, int gridID, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  double xlon1 = 0, xlon2 = 0, xlat1 = 0, xlat2 = 0;
  getlonlatparams(argc_offset, &xlon1, &xlon2, &xlat1, &xlat2);

  int gridtype = gridInqType(gridID);

  if ( gridtype == GRID_CURVILINEAR )
    genlonlatbox_curv(gridID, xlon1, xlon2, xlat1, xlat2, lat1, lat2, lon11, lon12, lon21, lon22);
  else
    genlonlatbox_reg(gridID, xlon1, xlon2, xlat1, xlat2, lat1, lat2, lon11, lon12, lon21, lon22);
}

static
int genlonlatgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  genlonlatbox(0, gridID1, lat1, lat2, lon11, lon12, lon21, lon22);

  int gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return gridID2;
}

static
int gencellgrid(int gridID1, int *gridsize2, int **cellidx)
{
  int cellinc = 4096;
  
  int argc_offset = 0;
  operatorCheckArgc(argc_offset+4);

  double xlon1 = parameter2double(operatorArgv()[argc_offset+0]);
  double xlon2 = parameter2double(operatorArgv()[argc_offset+1]);
  double xlat1 = parameter2double(operatorArgv()[argc_offset+2]);
  double xlat2 = parameter2double(operatorArgv()[argc_offset+3]);

  double x;
  if ( xlon1 >= xlon2 ) { x = xlon1; xlon1 = xlon2; xlon2 = x; }
  if ( xlat1 >= xlat2 ) { x = xlat1; xlat1 = xlat2; xlat2 = x; }

  int gridtype = gridInqType(gridID1);
  int gridsize1 = gridInqSize(gridID1);

  if ( gridtype != GRID_UNSTRUCTURED ) cdoAbort("Internal problem, wrong grid type!");

  double *xvals = (double*) Malloc(gridsize1*sizeof(double));
  double *yvals = (double*) Malloc(gridsize1*sizeof(double));
  gridInqXvals(gridID1, xvals);
  gridInqYvals(gridID1, yvals);

  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  gridInqXunits(gridID1, xunits);
  gridInqYunits(gridID1, yunits);

  double xfact = 1, yfact = 1;
  if ( strncmp(xunits, "radian", 6) == 0 ) xfact = RAD2DEG;
  if ( strncmp(yunits, "radian", 6) == 0 ) yfact = RAD2DEG;

  /* find gridsize2 */
  *cellidx = NULL;
  int maxcell = 0;
  int nvals = 0;
  for ( int i = 0; i < gridsize1; ++i )
    {
      double xval = xvals[i]*xfact;
      double yval = yvals[i]*yfact;
      if ( yval >= xlat1 && yval <= xlat2 )
	if ( (xval     >= xlon1 && xval     <= xlon2) ||
	     (xval+360 >= xlon1 && xval+360 <= xlon2) ||
	     (xval-360 >= xlon1 && xval-360 <= xlon2)  )
	  {
	    nvals++;
	    if ( nvals > maxcell )
	      {
		maxcell += cellinc;
		*cellidx = (int*) Realloc(*cellidx, maxcell*sizeof(int));
	      }
	    (*cellidx)[nvals-1] = i;
	  }
    }

  if ( nvals == 0 ) cdoAbort("No grid points found!");

  *gridsize2 = nvals;

  Free(xvals);
  Free(yvals);

  int gridID2 = gengridcell(gridID1, *gridsize2, *cellidx);

  return gridID2;
}


void genindexbox(int argc_offset, int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  operatorCheckArgc(argc_offset+4);

  *lon11 = parameter2int(operatorArgv()[argc_offset+0]);
  *lon12 = parameter2int(operatorArgv()[argc_offset+1]);
  *lat1  = parameter2int(operatorArgv()[argc_offset+2]);
  *lat2  = parameter2int(operatorArgv()[argc_offset+3]);

  if ( *lat1 > *lat2 )
    {
      int temp = *lat1;
      *lat1 = *lat2;
      *lat2 = temp;
    }

  int nlon = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  if ( *lat1 < 1 )
    {
      cdoWarning("First latitude index out of range, set to 1!");
      *lat1 = 1;
    }
  if ( *lat1 > nlat )
    {
      cdoWarning("First latitude index out of range, set to %d!", nlat);
      *lat1 = nlat;
    }
  if ( *lat2 < 1 )
    {
      cdoWarning("First latitude index out of range, set to 1!");
      *lat2 = 1;
    }
  if ( *lat2 > nlat )
    {
      cdoWarning("Last latitude index out of range, set to %d!", nlat);
      *lat2 = nlat;
    }
  if ( *lon11 < 1 )
    {
      cdoWarning("First longitude index out of range, set to 1!");
      *lon11 = 1;
    }
  if ( *lon11 > nlon )
    {
      cdoWarning("First longitude index out of range, set to %d!", nlon);
      *lon11 = nlon;
    }
  if ( *lon12 < 1 )
    {
      cdoWarning("Last longitude index out of range, set to 1!");
      *lon12 = 1;
    }
  if ( *lon12 > nlon+1 )
    {
      cdoWarning("Last longitude index out of range, set to %d!", nlon);
      *lon12 = nlon;
    }

  (*lon11)--;
  (*lon12)--;
  (*lat1)--;
  (*lat2)--;

  if ( *lon11 > *lon12 )
    {
      *lon21 = *lon11;
      *lon22 = nlon - 1;
      *lon11 = 0;
    }
  else
    {
      if ( *lon12 > nlon-1 )
	{
	  *lon21 = *lon11;
	  *lon22 = nlon - 1;
	  *lon11 = 0;
	  *lon12 = 0;
	}
      else
	{
	  *lon21 = 0;
	  *lon22 = -1;
	}
    }
}

static
int genindexgrid(int gridID1, int *lat1, int *lat2, int *lon11, int *lon12, int *lon21, int *lon22)
{
  genindexbox(0, gridID1, lat1, lat2, lon11, lon12, lon21, lon22);

  int gridID2 = -1;
  if ( gridInqType(gridID1) == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_LCC )
    gridID2 = cdo_define_subgrid_grid(gridID1, *lon11, *lon12, *lat1, *lat2);
  else
    gridID2 = gengrid(gridID1, *lat1, *lat2, *lon11, *lon12, *lon21, *lon22);

  return gridID2;
}

static
void window(int nwpv, double *array1, int gridID1, double *array2,
	    long lat1, long lat2, long lon11, long lon12, long lon21, long lon22)
{
  long ilat, ilon;
  long nlon = gridInqXsize(gridID1);

  if ( nwpv == 2 )
    {
      for ( ilat = lat1; ilat <= lat2; ilat++ )
	{
	  for ( ilon = lon21; ilon <= lon22; ilon++ )
	    {
	      *array2++ = array1[ilat*nlon*2 + ilon*2];
	      *array2++ = array1[ilat*nlon*2 + ilon*2+1];
	    }
	  for ( ilon = lon11; ilon <= lon12; ilon++ )
	    {
	      *array2++ = array1[ilat*nlon*2 + ilon*2];
	      *array2++ = array1[ilat*nlon*2 + ilon*2+1];
	    }
	}
    }
  else
    {
      for ( ilat = lat1; ilat <= lat2; ilat++ )
	{
	  for ( ilon = lon21; ilon <= lon22; ilon++ )
	    *array2++ = array1[ilat*nlon + ilon];
	  for ( ilon = lon11; ilon <= lon12; ilon++ )
	    *array2++ = array1[ilat*nlon + ilon];
	}
    }
}

static
void window_cell(int nwpv, double *array1, int gridID1, double *array2, long gridsize2, int *cellidx)
{
  if ( nwpv == 2 )
    {
      for ( long i = 0; i < gridsize2; ++i )
	{
	  array2[i*2]   = array1[cellidx[i]*2];
	  array2[i*2+1] = array1[cellidx[i]*2+1];
	}
    }
  else
    {
      for ( long i = 0; i < gridsize2; ++i )
	array2[i] = array1[cellidx[i]];
    }
}


void *Selbox(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridID1 = -1, gridID2;
  int index, gridtype = -1;
  int nmiss;
  double missval;
  typedef struct {
    int gridID1, gridID2;
    int *cellidx, nvals;
    int lat1, lat2, lon11, lon12, lon21, lon22; 
  } sbox_t;

  cdoInitialize(argument);

  // clang-format off
  int SELLONLATBOX = cdoOperatorAdd("sellonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude");
  int SELINDEXBOX  = cdoOperatorAdd("selindexbox",  0, 0, "index of first and last longitude and index of first and last latitude");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  bool *vars = (bool*) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;

  int ngrids = vlistNgrids(vlistID1);
  sbox_t *sbox = (sbox_t *) Malloc(ngrids*sizeof(sbox_t));

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
           (gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL) ||
           (gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_LCC) ||
	   (operatorID == SELINDEXBOX && gridtype == GRID_GENERIC && 
	    gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) ||
	   (operatorID == SELLONLATBOX && gridtype == GRID_UNSTRUCTURED) )
	{
	  if ( operatorID == SELLONLATBOX )
	    {
	      int gridsize = gridInqSize(gridID1);
	      if ( gridsize == 1 ) continue;

	      if ( gridtype == GRID_UNSTRUCTURED )
		gridID2 = gencellgrid(gridID1, &sbox[index].nvals, &sbox[index].cellidx);
	      else
		gridID2 = genlonlatgrid(gridID1, &sbox[index].lat1, &sbox[index].lat2, &sbox[index].lon11, 
					&sbox[index].lon12, &sbox[index].lon21, &sbox[index].lon22);
	    }
	  else
	    gridID2 = genindexgrid(gridID1, &sbox[index].lat1, &sbox[index].lat2, &sbox[index].lon11, 
				   &sbox[index].lon12, &sbox[index].lon21, &sbox[index].lon22);
	  
	  sbox[index].gridID1 = gridID1;
	  sbox[index].gridID2 = gridID2;

	  vlistChangeGridIndex(vlistID2, index, gridID2);

	  for ( varID = 0; varID < nvars; varID++ )
	    if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	      vars[varID] = true;
	}
      else if ( gridtype == GRID_GENERIC && gridInqXsize(gridID1) <= 1 && gridInqYsize(gridID1) <=1 )
	{
	}
      else
	{
	  cdoPrint("Unsupported grid type: %s", gridNamePtr(gridtype));
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoPrint("Use option -R to convert Gaussian reduced grid to a regular grid!");
	  cdoAbort("Unsupported grid type!");
	}
    }

  if ( cdoVerbose )
    {
      if ( gridtype != GRID_UNSTRUCTURED )
	{
	  cdoPrint("box1 - idx1,idx2,idy1,idy2: %d,%d,%d,%d",
		   sbox[0].lon21+1, sbox[0].lon22+1, sbox[0].lat1+1, sbox[0].lat2+1);
	  cdoPrint("box2 - idx1,idx2,idy1,idy2: %d,%d,%d,%d",
		   sbox[0].lon11+1, sbox[0].lon12+1, sbox[0].lat1+1, sbox[0].lat2+1);
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    if ( vars[varID] ) break;

  if ( varID >= nvars ) cdoWarning("No variables selected!");

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  int gridsize2 = vlistGridsizeMax(vlistID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double*) Malloc(gridsize2*sizeof(double));

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

	  if ( vars[varID] )
	    {
              // number of words per value; real:1  complex:2
	      int nwpv = vlistInqNWPV(vlistID1, varID);
	      
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      for ( index = 0; index < ngrids; index++ )
		if ( gridID1 == sbox[index].gridID1 ) break;

	      if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

	      gridsize2 = gridInqSize(sbox[index].gridID2);

	      if ( operatorID == SELLONLATBOX && gridtype == GRID_UNSTRUCTURED )
		window_cell(nwpv, array1, gridID1, array2, gridsize2, sbox[index].cellidx);
 	      else
		window(nwpv, array1, gridID1, array2, sbox[index].lat1, sbox[index].lat2, sbox[index].lon11, 
		       sbox[index].lon12, sbox[index].lon21, sbox[index].lon22);

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID2, varID);
		  for ( int i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

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

  vlistDestroy(vlistID2);

  if ( vars   ) Free(vars);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  if ( sbox )
    {
      if ( operatorID == SELLONLATBOX  && gridtype == GRID_UNSTRUCTURED )
	{
	  for ( index = 0; index < ngrids; index++ )
	    if ( sbox[index].cellidx ) Free(sbox[index].cellidx);
	}
      Free(sbox);
    }

  cdoFinish();

  return 0;
}
