/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.m1pg.de>
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

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Arithlat(void *argument)
{
  int gridtype;
  int gridID0 = -1;
  int nrecs;
  int varID, levelID;
  int nmiss;
  long i;
  char units[CDI_MAX_NAME];
  double *scale = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("mulcoslat", func_mul, 0, NULL);
  cdoOperatorAdd("divcoslat", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  long gridsize = vlistGridsizeMax(vlistID1);

  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array, &nmiss);
	  
	  int gridID = vlistInqVarGrid(vlistID1, varID);

	  if ( gridID != gridID0 )
	    {
	      gridID0 = gridID;

	      gridtype = gridInqType(gridID);
              int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
	      if ( gridtype == GRID_LONLAT      ||
		   gridtype == GRID_GAUSSIAN    ||
		   projtype == CDI_PROJ_LCC )
		{
		  gridID = gridToCurvilinear(gridID, 0);
		}
	      else if ( gridtype == GRID_CURVILINEAR ||
			gridtype == GRID_UNSTRUCTURED )
		{
		  /* No conversion necessary */
		}
	      else if ( gridtype == GRID_GME )
		{
		  gridID = gridToUnstructured(gridID, 0);
		}
	      else
		{
		  if ( gridtype == GRID_GAUSSIAN_REDUCED )
		    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
			     gridNamePtr(gridtype));
		  else
		    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
		}

	      gridsize = gridInqSize(gridID);

	      scale = (double*) Realloc(scale, gridsize*sizeof(double));
	      gridInqYvals(gridID, scale);

	      /* Convert lat/lon units if required */
	      
	      gridInqXunits(gridID, units);

	      grid_to_radian(units, gridsize, scale, "grid latitudes");

	      if ( operfunc == func_mul )
		for ( i = 0; i < gridsize; ++i ) scale[i] = cos(scale[i]);
	      else
		for ( i = 0; i < gridsize; ++i ) scale[i] = 1./cos(scale[i]);

	      if ( cdoVerbose ) for ( i = 0; i < 10; ++i ) cdoPrint("coslat  %3d  %g", i+1, scale[i]);
	    }

	  for ( i = 0; i < gridsize; ++i ) array[i] *= scale[i];

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array ) Free(array);
  if ( scale ) Free(scale);

  cdoFinish();

  return 0;
}
