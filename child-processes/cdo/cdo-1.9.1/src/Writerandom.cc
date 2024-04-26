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

      Writerandom writerandom
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Writerandom(void *argument)
{
  int gridsize;
  int nrecs;
  int varID, levelID;
  int rindex;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      double **recdata = (double**) Malloc(nrecs*sizeof(double*));
      int *recvarID   = (int*) Malloc(nrecs*sizeof(int));
      int *reclevelID = (int*) Malloc(nrecs*sizeof(int));
      int *recnmiss   = (int*) Malloc(nrecs*sizeof(int));
      int *recindex   = (int*) Malloc(nrecs*sizeof(int));

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  recvarID[recID] = varID;
	  reclevelID[recID] = levelID;
	  recdata[recID] = (double*) Malloc(gridsize*sizeof(double));
	  pstreamReadRecord(streamID1, recdata[recID], &recnmiss[recID]);
	}

      for ( int recID = 0; recID < nrecs; recID++ ) recindex[recID] = -1;

      for ( rindex = nrecs-1; rindex >= 0; rindex-- )
	{
	  int index = (int) (rindex*((double)rand())/((double)RAND_MAX));
	  /*	printf("rindex %d %d\n", rindex, index); */
	  int ipos = -1;
	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      if ( recindex[recID] == -1 ) ipos++;
	      if ( recindex[recID] == -1 && ipos == index )
		{
		  recindex[recID] = rindex;
		  break;
		}
	    }
	}

      /*
      for ( int recID = 0; recID < nrecs; recID++ )
	printf("recID %d %d\n", recID, recindex[recID]);
      */
      for ( int recID = 0; recID < nrecs; recID++ )
	if ( recindex[recID] == -1 )
	  cdoAbort("Internal problem! Random initialize.");

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  rindex   = recindex[recID];
	  varID    = recvarID[rindex];
	  levelID  = reclevelID[rindex];
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, recdata[rindex], recnmiss[rindex]);
	}

      for ( int recID = 0; recID < nrecs; recID++ ) Free(recdata[recID]);

      Free(recdata);
      Free(recvarID);
      Free(reclevelID);
      Free(recnmiss);
      Free(recindex);

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
