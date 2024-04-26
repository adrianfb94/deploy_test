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
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024


void *Duplicate(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nalloc = 0;
  int nmiss;
  int *vdate = NULL, *vtime = NULL;
  int ndup = 2;
  field_type ***vars = NULL;

  cdoInitialize(argument);

  if      ( operatorArgc()  > 1 ) cdoAbort("Too many arguments!");
  else if ( operatorArgc() == 1 ) ndup = parameter2int(operatorArgv()[0]);

  if ( cdoVerbose ) cdoPrint("ndup = %d", ndup);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ntsteps  = vlistNtsteps(vlistID1);
  int nvars    = vlistNvars(vlistID1);

  if ( ntsteps == 1 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	if ( vlistInqVarTimetype(vlistID1, varID) != TIME_CONSTANT ) break;

      if ( varID == nvars ) ntsteps = 0;
    }

  if ( ntsteps == 0 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }
 
  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vdate = (int*) Realloc(vdate, nalloc*sizeof(int));
	  vtime = (int*) Realloc(vtime, nalloc*sizeof(int));
	  vars  = (field_type***) Realloc(vars, nalloc*sizeof(field_type**));
	}

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  int gridID   = vlistInqVarGrid(vlistID1, varID);
	  int gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
	  pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  int nts = tsID;

  for ( int idup = 0; idup < ndup; idup++ )
    {
      for ( tsID = 0; tsID < nts; tsID++ )
	{
	  taxisDefVdate(taxisID2, vdate[tsID]);
	  taxisDefVtime(taxisID2, vtime[tsID]);
	  pstreamDefTimestep(streamID2, idup*nts+tsID);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( vars[tsID][varID][levelID].ptr )
		    {
		      nmiss = vars[tsID][varID][levelID].nmiss;
		      pstreamDefRecord(streamID2, varID, levelID);
		      pstreamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}
    }

  for ( tsID = 0; tsID < nts; tsID++ ) field_free(vars[tsID], vlistID1);

  if ( vars  ) Free(vars);
  if ( vdate ) Free(vdate);
  if ( vtime ) Free(vtime);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
