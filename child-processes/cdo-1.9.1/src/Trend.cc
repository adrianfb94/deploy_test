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

      Trend      trend           Trend
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Trend(void *argument)
{
  int vdate = 0, vtime = 0;
  int varID, levelID;
  int nmiss;
  int nrecs;
  double temp1, temp2;
  enum {nwork = 5};
  field_type **work[5];

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID2);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));

  for ( int w = 0; w < nwork; w++ )
    work[w] = field_calloc(vlistID1, FIELD_PTR);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      double zj = (double)tsID;
      
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);

	  double missval  = vlistInqVarMissval(vlistID1, varID);
	  int gridID   = vlistInqVarGrid(vlistID1, varID);
	  size_t gridsize = gridInqSize(gridID);

	  for ( size_t i = 0; i < gridsize; i++ )
	    if ( !DBL_IS_EQUAL(field1.ptr[i], missval) )
	      {
		work[0][varID][levelID].ptr[i] += zj;
		work[1][varID][levelID].ptr[i] += zj * zj;
		work[2][varID][levelID].ptr[i] += zj * field1.ptr[i];
		work[3][varID][levelID].ptr[i] += field1.ptr[i];
		work[4][varID][levelID].ptr[i]++;
	      }      
	}

      tsID++;
    }
	  

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  pstreamDefTimestep(streamID2, 0);
  pstreamDefTimestep(streamID3, 0);

  for ( int recID = 0; recID < nrecords; recID++ )
    {
      varID   = recVarID[recID];
      levelID = recLevelID[recID];

      double missval  = vlistInqVarMissval(vlistID1, varID);
      int gridID   = vlistInqVarGrid(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);

      double missval1  = missval;
      double missval2  = missval;

      for ( size_t i = 0; i < gridsize; i++ )
	{
	  temp1 = SUBMN(work[2][varID][levelID].ptr[i],
		      DIVMN( MULMN(work[0][varID][levelID].ptr[i], work[3][varID][levelID].ptr[i]), work[4][varID][levelID].ptr[i]));
	  temp2 = SUBMN(work[1][varID][levelID].ptr[i],
		      DIVMN( MULMN(work[0][varID][levelID].ptr[i], work[0][varID][levelID].ptr[i]), work[4][varID][levelID].ptr[i]));

	  field2.ptr[i] = DIVMN(temp1, temp2);
	  field1.ptr[i] = SUBMN( DIVMN(work[3][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]),
			      MULMN( DIVMN(work[0][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]), field2.ptr[i]));
	}

      nmiss = 0;
      for ( size_t i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(field1.ptr[i], missval) ) nmiss++;

      pstreamDefRecord(streamID2, varID, levelID);
      pstreamWriteRecord(streamID2, field1.ptr, nmiss);

      nmiss = 0;
      for ( size_t i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(field2.ptr[i], missval) ) nmiss++;

      pstreamDefRecord(streamID3, varID, levelID);
      pstreamWriteRecord(streamID3, field2.ptr, nmiss);
    }


  for ( int w = 0; w < nwork; w++ ) field_free(work[w], vlistID1);

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
