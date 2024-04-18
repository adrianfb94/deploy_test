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


void *Fldrms(void *argument)
{
  int lastgrid = -1;
  int index;
  int nrecs;
  int varID, levelID;
  int nmiss;
  double sglval;

  cdoInitialize(argument);

  bool needWeights = true;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  double slon = 0;
  double slat = 0;
  int gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  index = 0;
  int gridID1 = vlistGrid(vlistID1, index);
  int gridID2 = vlistGrid(vlistID2, index);

  if ( gridInqSize(gridID1) != gridInqSize(gridID2) )
    cdoAbort("Fields have different grid size!");
  
  if ( needWeights &&
       gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID3, index, gridID3);

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  field_type field1, field2, field3;
  field_init(&field1);
  field_init(&field2);
  field_init(&field3);

  int lim = vlistGridsizeMax(vlistID1);
  field1.ptr    = (double*) Malloc(lim*sizeof(double));
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) Malloc(lim*sizeof(double));

  field2.ptr    = (double*) Malloc(lim*sizeof(double));
  field2.weight = NULL;

  field3.ptr  = &sglval;
  field3.grid = gridID3;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      nrecs = pstreamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = (size_t) nmiss;
	  pstreamInqRecord(streamID2, &varID, &levelID);
	  pstreamReadRecord(streamID2, field2.ptr, &nmiss);
          field2.nmiss = (size_t) nmiss;

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);

          if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      field1.weight[0] = 1;
	      if ( field1.size > 1 )
		{
		  int wstatus = gridWeights(field1.grid, field1.weight);
		  if ( wstatus != 0 && tsID == 0 && levelID == 0 )
		    {
		      char varname[CDI_MAX_NAME];
		      vlistInqVarName(vlistID1, varID, varname);
		      cdoWarning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", varname);
		    }
		}
	    }

	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);
	  field3.missval = vlistInqVarMissval(vlistID1, varID);

	  fldrms(field1, field2, &field3);

	  pstreamDefRecord(streamID3, varID,  levelID);
	  pstreamWriteRecord(streamID3, &sglval, (int)field3.nmiss);
	}
      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( field1.ptr )    Free(field1.ptr);
  if ( field1.weight ) Free(field1.weight);
  if ( field2.ptr )    Free(field2.ptr);

  cdoFinish();

  return 0;
}
