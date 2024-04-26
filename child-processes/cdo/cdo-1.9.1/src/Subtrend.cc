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

     Subtrend   subtrend        Subtract trend
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Subtrend(void *argument)
{
  int gridID, varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));
  int streamID3 = pstreamOpenRead(cdoStreamName(2));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = pstreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_DIM);
  vlistCompare(vlistID1, vlistID3, CMP_DIM);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = pstreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field4;
  field_init(&field1);
  field_init(&field4);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field4.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type **vars2 = field_malloc(vlistID1, FIELD_PTR);
  field_type **vars3 = field_malloc(vlistID1, FIELD_PTR);


  int tsID = 0;
  int nrecs = pstreamInqTimestep(streamID2, tsID);

  for ( int recID = 0; recID < nrecs; recID++ )
    {
      pstreamInqRecord(streamID2, &varID, &levelID);
      pstreamReadRecord(streamID2, vars2[varID][levelID].ptr, &nmiss);
    }

  tsID = 0;
  nrecs = pstreamInqTimestep(streamID3, tsID);

  for ( int recID = 0; recID < nrecs; recID++ )
    {
      pstreamInqRecord(streamID3, &varID, &levelID);
      pstreamReadRecord(streamID3, vars3[varID][levelID].ptr, &nmiss);
    }


  tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID4, taxisID1);
      pstreamDefTimestep(streamID4, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  double missval = vlistInqVarMissval(vlistID1, varID);
	  double missval1 = missval;
	  double missval2 = missval;
	  for ( int i = 0; i < gridsize; i++ )
	    field4.ptr[i] = SUBMN(field1.ptr[i], ADDMN(vars2[varID][levelID].ptr[i], MULMN(vars3[varID][levelID].ptr[i], tsID)));
    
	  nmiss = 0;
	  for ( int i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(field4.ptr[i], missval) ) nmiss++;

	  pstreamDefRecord(streamID4, varID, levelID);
	  pstreamWriteRecord(streamID4, field4.ptr, nmiss);
	}

      tsID++;
    }

  field_free(vars2, vlistID1);
  field_free(vars3, vlistID1);

  if ( field1.ptr ) Free(field1.ptr);
  if ( field4.ptr ) Free(field4.ptr);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
