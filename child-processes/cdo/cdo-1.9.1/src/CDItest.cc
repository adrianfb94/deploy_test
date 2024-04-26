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


void *CDItest(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int max_copy = 3;
  double s_utime, s_stime;
  double e_utime, e_stime;

  cdoInitialize(argument);

  bool lcopy = false;
  //bool lcopy = UNCHANGED_RECORD;

  int NCOPY = cdoOperatorAdd("ncopy",   0, 0, NULL);
  UNUSED(NCOPY);

  int operatorID = cdoOperatorID();
  UNUSED(operatorID);

  //  operatorInputArg("Number of copies");
  if ( operatorArgc() == 1 ) max_copy = parameter2int(operatorArgv()[0]);

  processStartTime(&s_utime, &s_stime);

  int n = 0;
  while ( TRUE )
    {
      int streamID1 = pstreamOpenRead(cdoStreamName(0));

      int vlistID1 = pstreamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      int vlistID2 = vlistDuplicate(vlistID1);
      int taxisID2 = taxisDuplicate(taxisID1);
      vlistDefTaxis(vlistID2, taxisID2);

      int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
      pstreamDefVlist(streamID2, vlistID2);

      int gridsize = vlistGridsizeMax(vlistID1);
      double *array = (double*) Malloc(gridsize*sizeof(double));

      int tsID1 = 0;
      int tsID2 = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  pstreamDefTimestep(streamID2, tsID2);
	       
	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);
	      pstreamDefRecord(streamID2,  varID,  levelID);
	  
	      if ( lcopy )
		{
		  pstreamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  pstreamReadRecord(streamID1, array, &nmiss);
		  pstreamWriteRecord(streamID2, array, nmiss);
		}
	    }
	  tsID1++;
	  tsID2++;
	}

      pstreamClose(streamID1);
      pstreamClose(streamID2);

      vlistDestroy(vlistID2);
      taxisDestroy(taxisID2);

      if ( array ) Free(array);

      n++;

      cdoProcessTime(&e_utime, &e_stime);

      double c_usertime = e_utime - s_utime;
      double c_systime  = e_stime - s_stime;
      double c_cputime  = c_usertime + c_systime;

      s_utime = e_utime;
      s_stime = e_stime;

      cdoPrint("Copy number %d: %.2fs %.2fs %.2fs", n, c_usertime, c_systime, c_cputime);

      if ( n == max_copy ) break;
    }

  cdoFinish();

  return 0;
}
