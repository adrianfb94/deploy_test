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

*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


static
int gentpngrid(int gridID1)
{
  int gridtype, gridID2;
  int nlon1, nlat1;
  int nlon2, nlat2;
  int prec;
  int ilat, ilon, ilonr, k, kr;
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;
  double *xbounds1 = NULL, *ybounds1 = NULL;
  double *xbounds2 = NULL, *ybounds2 = NULL;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = nlon1;
  nlat2 = nlat1+2;

  gridtype = gridInqType(gridID1);
  prec     = gridInqDatatype(gridID1);

  gridID2 = gridCreate(gridtype, nlon2*nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefDatatype(gridID2, prec);

  grid_copy_attributes(gridID1, gridID2);

  char xunits[CDI_MAX_NAME]; xunits[0] = 0;
  char yunits[CDI_MAX_NAME]; yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
	
  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xvals1 = (double*) Malloc(nlon1*nlat1*sizeof(double));
	  yvals1 = (double*) Malloc(nlon1*nlat1*sizeof(double));
	  xvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));

	  gridInqXvals(gridID1, xvals1);
	  gridInqYvals(gridID1, yvals1);

	  for ( ilat = 0; ilat < nlat1; ilat++ )
	    {
	      for ( ilon = 0; ilon < nlon1; ilon++ )
		{
		  xvals2[(ilat+2)*nlon1 + ilon] = xvals1[ilat*nlon1 + ilon];
		  yvals2[(ilat+2)*nlon1 + ilon] = yvals1[ilat*nlon1 + ilon];
		}
	    }

	  for ( ilon = 0; ilon < nlon1; ilon++ )
	    {
	      ilonr = nlon1 - ilon - 1;
	      xvals2[1*nlon1 + ilon] = xvals2[2*nlon1 + ilonr]; /* syncronise line 2 with line 3 */
	      xvals2[0*nlon1 + ilon] = xvals2[3*nlon1 + ilonr]; /* syncronise line 1 with line 4 */
	      yvals2[1*nlon1 + ilon] = yvals2[2*nlon1 + ilonr]; /* syncronise line 2 with line 3 */
	      yvals2[0*nlon1 + ilon] = yvals2[3*nlon1 + ilonr]; /* syncronise line 1 with line 4 */
	    }	
	  
	  gridDefXvals(gridID2, xvals2);
	  gridDefYvals(gridID2, yvals2);
	  
	  Free(xvals1);
	  Free(yvals1);
	  Free(xvals2);
	  Free(yvals2);
	}
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = (double*) Malloc(4*nlon1*nlat1*sizeof(double));
	  ybounds1 = (double*) Malloc(4*nlon1*nlat1*sizeof(double));
	  xbounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));

	  gridInqXbounds(gridID1, xbounds1);
	  gridInqYbounds(gridID1, ybounds1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      gridDefNvertex(gridID2, 4);

	      for ( ilat = 0; ilat < nlat1; ilat++ )
		{
		  for ( ilon = 0; ilon < 4*nlon1; ilon++ )
		    {
		      xbounds2[4*(ilat+2)*nlon1 + ilon] = xbounds1[4*ilat*nlon1 + ilon];
		      ybounds2[4*(ilat+2)*nlon1 + ilon] = ybounds1[4*ilat*nlon1 + ilon];
		    }
		}

	      for ( ilon = 0; ilon < nlon1; ilon++ )
		{
		  ilonr = nlon1 - ilon - 1;
		  for ( k = 0; k < 4; ++k )
		    {
		      kr = 3 - k;
		      xbounds2[4*1*nlon1 + 4*ilon + k] = xbounds2[4*2*nlon1 + 4*ilonr + kr];
		      xbounds2[4*0*nlon1 + 4*ilon + k] = xbounds2[4*3*nlon1 + 4*ilonr + kr];
		      ybounds2[4*1*nlon1 + 4*ilon + k] = ybounds2[4*2*nlon1 + 4*ilonr + kr];
		      ybounds2[4*0*nlon1 + 4*ilon + k] = ybounds2[4*3*nlon1 + 4*ilonr + kr];
		    }
		}	
	      /*
	      for ( ilon = 0; ilon < 4*nlon1; ilon++ )
		{
		  ilonr = 4*nlon1 - ilon - 1;
		    {
		      xbounds2[4*1*nlon1 + ilon ] = xbounds2[4*2*nlon1 + ilonr ];
		      xbounds2[4*0*nlon1 + ilon ] = xbounds2[4*3*nlon1 + ilonr ];
		      ybounds2[4*1*nlon1 + ilon ] = ybounds2[4*2*nlon1 + ilonr ];
		      ybounds2[4*0*nlon1 + ilon ] = ybounds2[4*3*nlon1 + ilonr ]; 
		    }
		}
	      */	
	    }

	  gridDefXbounds(gridID2, xbounds2);
	  gridDefYbounds(gridID2, ybounds2);

	  Free(xbounds1);
	  Free(ybounds1);
	  Free(xbounds2);
	  Free(ybounds2);
	}
    }

  return gridID2;
}

static
int gengrid(int gridID1, int lhalo, int rhalo)
{
  int gridtype, gridID2;
  int nlon1, nlat1;
  int nlon2, nlat2;
  int nmin, nmax;
  int i;
  int prec;
  int ilat, ilon;
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;
  double *xbounds1 = NULL, *ybounds1 = NULL;
  double *xbounds2 = NULL, *ybounds2 = NULL;
  double *pxvals2 = NULL, *pyvals2 = NULL;
  double *pxbounds2 = NULL, *pybounds2 = NULL;
  double cpi2 = M_PI*2;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = nlon1 + lhalo + rhalo;
  nlat2 = nlat1;

  nmin = 0;
  nmax = nlon1;
  if ( lhalo < 0 ) nmin  = -lhalo;
  if ( rhalo < 0 ) nmax +=  rhalo;
  /*
  printf("nlon1=%d, nlon2=%d, lhalo=%d, rhalo=%d nmin=%d, nmax=%d\n",
	 nlon1, nlon2, lhalo, rhalo, nmin, nmax);
  */
  gridtype = gridInqType(gridID1);
  prec     = gridInqDatatype(gridID1);

  gridID2 = gridCreate(gridtype, nlon2*nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefDatatype(gridID2, prec);

  grid_copy_attributes(gridID1, gridID2);

  char xunits[CDI_MAX_NAME]; xunits[0] = 0;
  char yunits[CDI_MAX_NAME]; yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  if ( memcmp(xunits, "degree", 6) == 0 ) cpi2 *= RAD2DEG;

  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xvals1 = (double*) Malloc(nlon1*nlat1*sizeof(double));
	  yvals1 = (double*) Malloc(nlon1*nlat1*sizeof(double));
	  xvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = (double*) Malloc(nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xvals1 = (double*) Malloc(nlon1*sizeof(double));
	  yvals1 = (double*) Malloc(nlat1*sizeof(double));
	  xvals2 = (double*) Malloc(nlon2*sizeof(double));
	  yvals2 = (double*) Malloc(nlat2*sizeof(double));
	}

      pxvals2 = xvals2;
      pyvals2 = yvals2;

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  for ( ilat = 0; ilat < nlat2; ilat++ )
	    {
	      for ( ilon = nlon1-lhalo; ilon < nlon1; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}

	      for ( ilon = nmin; ilon < nmax; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < rhalo; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}
	    }
	}
      else
	{
	  for ( i = nlon1-lhalo; i < nlon1; i++ ) *pxvals2++ = xvals1[i] - cpi2;
	  for ( i = nmin; i < nmax; i++ ) *pxvals2++ = xvals1[i];
	  for ( i = 0; i < rhalo; i++ ) *pxvals2++ = xvals1[i] + cpi2;

	  for ( i = 0; i < nlat1; i++ ) yvals2[i] = yvals1[i];
	}
      /*
      for ( i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
      for ( i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
      */
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      Free(xvals1);
      Free(yvals1);
      Free(xvals2);
      Free(yvals2);
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = (double*) Malloc(4*nlon1*nlat1*sizeof(double));
	  ybounds1 = (double*) Malloc(4*nlon1*nlat1*sizeof(double));
	  xbounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = (double*) Malloc(4*nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xbounds1 = (double*) Malloc(2*nlon1*sizeof(double));
	  ybounds1 = (double*) Malloc(2*nlat1*sizeof(double));
	  xbounds2 = (double*) Malloc(2*nlon2*sizeof(double));
	  ybounds2 = (double*) Malloc(2*nlat2*sizeof(double));
	}

      pxbounds2 = xbounds2;
      pybounds2 = ybounds2;

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  gridDefNvertex(gridID2, 4);
	  for ( ilat = 0; ilat < nlat1; ilat++ )
	    {
	      for ( ilon = 4*(nlon1-lhalo); ilon < 4*nlon1; ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon1 + ilon];
		}

	      for ( ilon = 4*nmin; ilon < 4*nmax; ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < 4*rhalo; ilon++ )
		{
		  *pxbounds2++ = xbounds1[4*ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[4*ilat*nlon1 + ilon];
		}
	    }
	}
      else
	{
	  gridDefNvertex(gridID2, 2);
	  for ( i = 2*(nlon1-lhalo); i < 2*nlon1; i++ ) *pxbounds2++ = xbounds1[i] - cpi2;
	  for ( i = 2*nmin; i < 2*nmax; i++ ) *pxbounds2++ = xbounds1[i];
	  for ( i = 0; i < 2*rhalo; i++ ) *pxbounds2++ = xbounds1[i] + cpi2;

	  for ( i = 0; i < 2*nlat2; i++ ) ybounds2[i] = ybounds1[i];
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


static
int genindexgrid(int gridID1, int *lhalo, int *rhalo)
{
  operatorCheckArgc(2);

  *lhalo = parameter2int(operatorArgv()[0]);
  *rhalo = parameter2int(operatorArgv()[1]);

  int nlon1 = gridInqXsize(gridID1);

  if ( *lhalo > nlon1 )
    {
      *lhalo = nlon1;
      cdoWarning("left halo out of range. Set to %d.", *lhalo);
    }

  if ( *lhalo < 0 && -(*lhalo) > nlon1/2 )
    {
      *lhalo = -nlon1/2;
      cdoWarning("left halo out of range. Set to %d.", *lhalo);
    }

  if ( *rhalo > nlon1 )
    {
      *rhalo = nlon1;
      cdoWarning("right halo out of range. Set to %d.", rhalo);
    }

  if ( *rhalo < 0 && -(*rhalo) > nlon1/2 )
    {
      *rhalo = -nlon1/2;
      cdoWarning("right halo out of range. Set to %d.", rhalo);
    }

  int gridID2 = gengrid(gridID1, *lhalo, *rhalo);

  return gridID2;
}


static
void halo(double *array1, int gridID1, double *array2, int lhalo, int rhalo)
{
  int nlon1 = gridInqXsize(gridID1);
  int nlat  = gridInqYsize(gridID1);

  int nmin = 0;
  int nmax = nlon1;
  if ( lhalo < 0 ) nmin  = -lhalo;
  if ( rhalo < 0 ) nmax +=  rhalo;

  for ( int ilat = 0; ilat < nlat; ilat++ )
    {
      for ( int ilon = nlon1-lhalo; ilon < nlon1; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];

      for ( int ilon = nmin; ilon < nmax; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];

      for ( int ilon = 0; ilon < rhalo; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];
    }
}


static
void tpnhalo(double *array1, int gridID1, double *array2)
{
  int nlon = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  for ( int ilat = 0; ilat < nlat; ilat++ )
    for ( int ilon = 0; ilon < nlon; ilon++ )
      array2[(ilat+2)*nlon + ilon] = array1[ilat*nlon + ilon];

  for ( int ilon = 0; ilon < nlon; ilon++ )
    {
      int ilonr = nlon - ilon - 1;
      array2[1*nlon + ilon] = array2[2*nlon + ilonr]; /* syncronise line 2 with line 3 */
      array2[0*nlon + ilon] = array2[3*nlon + ilonr]; /* syncronise line 1 with line 4 */
    }
}


void *Sethalo(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridsize, gridsize2;
  int gridID1 = -1, gridID2;
  int index, gridtype;
  int nmiss;
  int i;
  int lhalo = 0, rhalo = 0;
  double missval;

  cdoInitialize(argument);

  // clang-format off
  int SETHALO = cdoOperatorAdd("sethalo", 0, 0, NULL);
                cdoOperatorAdd("tpnhalo", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) break;
      if ( gridtype == GRID_CURVILINEAR ) break;
      if ( gridtype == GRID_GENERIC &&
	   gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0 ) break;
    }

  if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if ( index == ngrids ) cdoAbort("No regular grid found!");
  if ( ndiffgrids > 0 )  cdoAbort("Too many different grids!");

  if ( operatorID == SETHALO )
    {
      operatorInputArg("left and right halo");
      gridID2 = genindexgrid(gridID1, &lhalo, &rhalo);
    }
  else
    {
      gridID2 = gentpngrid(gridID1);
    }

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      if ( gridID1 == vlistGrid(vlistID1, index) )
	{
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	  break;
	}
    }

  int nvars = vlistNvars(vlistID1);
  int *vars  = (int*) Malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  gridsize = gridInqSize(gridID1);
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  gridsize2 = gridInqSize(gridID2);
  double *array2 = (double*) Malloc(gridsize2*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      pstreamReadRecord(streamID1, array1, &nmiss);

	      if ( operatorID == SETHALO )
		halo(array1, gridID1, array2, lhalo, rhalo);
	      else
		tpnhalo(array1, gridID1, array2);

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID1, varID);
		  for ( i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      pstreamDefRecord(streamID2, varID, levelID);
	      pstreamWriteRecord(streamID2, array2, nmiss);
	    }
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( vars   ) Free(vars);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  cdoFinish();

  return 0;
}
