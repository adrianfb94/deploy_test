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
void gen_index(int gridID1, int gridID2, int *index)
{
  int i, j, k, i1, i2;

  int gridtype1 = gridInqType(gridID1);
  int gridtype2 = gridInqType(gridID2);

  int gridsize2 = gridInqSize(gridID2);

  if ( gridtype1 != gridtype2 )
    cdoAbort("Input streams have different grid types!");

  if ( index == NULL ) cdoAbort("Internal problem, index not allocated!");

  for ( i = 0; i < gridsize2; i++ ) index[i] = -1;

  if ( gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN )
    {
      int nlon1 = gridInqXsize(gridID1);
      int nlat1 = gridInqYsize(gridID1);

      int nlon2 = gridInqXsize(gridID2);
      int nlat2 = gridInqYsize(gridID2);

      if ( ! (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL)) )
	cdoAbort("Grid 1 has no values!");

      if ( ! (gridInqXvals(gridID2, NULL) && gridInqYvals(gridID2, NULL)) )
	cdoAbort("Grid 2 has no values!");

      double *xvals1 = (double*) Malloc(nlon1*sizeof(double));
      double *yvals1 = (double*) Malloc(nlat1*sizeof(double));
      double *xvals2 = (double*) Malloc(nlon2*sizeof(double));
      double *yvals2 = (double*) Malloc(nlat2*sizeof(double));

      int *xindex = (int*) Malloc(nlon2*sizeof(int));
      int *yindex = (int*) Malloc(nlat2*sizeof(int));

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      /* Convert lat/lon units if required */
      {
	char units[CDI_MAX_NAME];
	gridInqXunits(gridID1, units);
	grid_to_degree(units, nlon1, xvals1, "grid1 center lon");
	gridInqYunits(gridID1, units);
	grid_to_degree(units, nlat1, yvals1, "grid1 center lat");
      }

      gridInqXvals(gridID2, xvals2);
      gridInqYvals(gridID2, yvals2);

      /* Convert lat/lon units if required */
      {
	char units[CDI_MAX_NAME];
	gridInqXunits(gridID2, units);
	grid_to_degree(units, nlon2, xvals2, "grid2 center lon");
	gridInqYunits(gridID2, units);
	grid_to_degree(units, nlat2, yvals2, "grid2 center lat");
      }

      for ( i2 = 0; i2 < nlat2; i2++ )
	{
	  for ( i1 = 0; i1 < nlat1; i1++ )
	    if ( fabs(yvals2[i2]-yvals1[i1]) < 0.001 ) break;

	  if ( i1 == nlat1 )
	    yindex[i2] = -1;
	  else
            yindex[i2] = i1;
	}

      for ( i2 = 0; i2 < nlon2; i2++ )
	{
	  for ( i1 = 0; i1 < nlon1; i1++ )
	    if ( fabs(xvals2[i2]-xvals1[i1]) < 0.001 ) break;

	  if ( i1 == nlon1 )
	    {
	      if ( xvals2[i2] < 0 )
		{
		  for ( i1 = 0; i1 < nlon1; i1++ )
		    if ( fabs(xvals2[i2]+360-xvals1[i1]) < 0.001 ) break;
		}
	      else if ( xvals2[i2] > 180 )
		{
		  for ( i1 = 0; i1 < nlon1; i1++ )
		    if ( fabs(xvals2[i2]-360-xvals1[i1]) < 0.001 ) break;
		}
	    }

	  if ( i1 == nlon1 )
	    xindex[i2] = -1;
	  else
            xindex[i2] = i1;
	}
      /*
      for ( i2 = 0; i2 < nlon2; i2++ )
	printf("x %d %d\n", i2, xindex[i2]);

      for ( i2 = 0; i2 < nlat2; i2++ )
	printf("y %d %d\n", i2, yindex[i2]);
      */
      k = 0;
      for ( j = 0; j < nlat2; j++ )
	for ( i = 0; i < nlon2; i++ )
	  {
	    if ( xindex[i] == -1 || yindex[j] == -1 )
	      index[k++] = -1;
	    else
	      index[k++] = yindex[j]*nlon1 + xindex[i];
	  }


      Free(xindex);
      Free(yindex);
      Free(xvals1);
      Free(yvals1);
      Free(xvals2);
      Free(yvals2);
    }
  else
    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype1));
}


void *Mergegrid(void *argument)
{
  int varID, levelID;
  int nrecs = 0;
  int nmiss1, nmiss2;
  int index;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID2 = pstreamInqVlist(streamID2);

  vlistCompare(vlistID1, vlistID2, CMP_NAME | CMP_NLEVEL);

  int ndiffgrids = 0;
  for ( index = 1; index < vlistNgrids(vlistID1); index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index) )
      ndiffgrids++;

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids in %s!", cdoStreamName(0)->args);

  ndiffgrids = 0;
  for ( index = 1; index < vlistNgrids(vlistID2); index++ )
    if ( vlistGrid(vlistID2, 0) != vlistGrid(vlistID2, index))
      ndiffgrids++;

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids in %s!", cdoStreamName(1)->args);

  int gridID1 = vlistGrid(vlistID1, 0);
  int gridID2 = vlistGrid(vlistID2, 0);

  int gridsize1 = gridInqSize(gridID1);
  int gridsize2 = gridInqSize(gridID2);

  double *array1 = (double*) Malloc(gridsize1*sizeof(double));
  double *array2 = (double*) Malloc(gridsize2*sizeof(double));
  int *gindex = (int*) Malloc(gridsize2*sizeof(int));

  gen_index(gridID1, gridID2, gindex);

  int vlistID3 = vlistDuplicate(vlistID1);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistDefTaxis(vlistID3, taxisID3);
  pstreamDefVlist(streamID3, vlistID3);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      int nrecs2 = pstreamInqTimestep(streamID2, tsID);
      if ( nrecs2 == 0 )
	cdoAbort("Input streams have different number of timesteps!");

      if ( nrecs != nrecs2 )
	cdoAbort("Input streams have different number of records!");

      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID2, &varID, &levelID);
	  pstreamReadRecord(streamID2, array2, &nmiss2);

	  double missval2 = vlistInqVarMissval(vlistID2, varID);

	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss1);

	  double missval1 = vlistInqVarMissval(vlistID1, varID);

	  for ( int i = 0; i < gridsize2; i++ )
	    {
	      if ( gindex[i] >= 0 && !DBL_IS_EQUAL(array2[i], missval2) )
		{
		  array1[gindex[i]] = array2[i];
		}
	    }

	  if ( nmiss1 )
	    {
	      nmiss1 = 0;
	      for ( int i = 0; i < gridsize1; i++ )
		if ( DBL_IS_EQUAL(array1[i], missval1) ) nmiss1++;
	    }

	  pstreamDefRecord(streamID3, varID, levelID);
	  pstreamWriteRecord(streamID3, array1, nmiss1);
	}

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
 
  if ( gindex ) Free(gindex);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  cdoFinish();

  return 0;
}
