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


void *Deltime(void *argument)
{
  int nrecs;
  int varID, levelID;
  int vdate /*, vtime */;
  int copytimestep;
  int gridsize;
  int nmiss;
  int year, month, day;
  int dday, dmon;
  double *array = NULL;
  const char *cmons[]={"", "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"};

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int DELDAY   = cdoOperatorAdd("delday",   0, 0, NULL);
  int DEL29FEB = cdoOperatorAdd("del29feb", 0, 0, NULL);
  // clang-format on

  UNUSED(DELDAY);

  int operatorID = cdoOperatorID();

  if ( operatorID == DEL29FEB )
    {
      dday = 29;
      dmon = 2;
    }
  else
    {
      int im;
      int nsel;
      char *sarg;
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Too few arguments!");
      if ( nsel > 1 ) cdoAbort("Too many arguments!");
      sarg = operatorArgv()[0];
      dday = atoi(sarg);
      dmon = 0;
      while ( isdigit(*sarg) ) sarg++;
      if ( isalpha(*sarg) )
	{
	  char smon[32];
	  strncpy(smon, sarg, sizeof(smon)-1);
	  smon[sizeof(smon)-1] = 0;
	  strtolower(smon);
	  for ( im = 0; im < 12; ++im )
	    if ( memcmp(smon, cmons[im+1], 3) == 0 ) break;

	  if ( im < 12 ) dmon = im + 1;
	}
    }

  if ( cdoVerbose ) cdoPrint("delete day %d%s", dday, cmons[dmon]);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisDefCalendar(taxisID2, CALENDAR_365DAYS);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) Malloc(gridsize*sizeof(double));
    }
      
  int nfound = 0;
  int tsID  = 0;
  int tsID2 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      // vtime = taxisInqVtime(taxisID1);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( day == dday && (month == dmon || dmon == 0) )
	{
	  nfound++;
	  copytimestep = FALSE;
	  if ( cdoVerbose )
	    cdoPrint("Delete %4.4d-%2.2d-%2.2d at timestep %d", year, month, day, tsID+1);
	}
      else
	copytimestep = TRUE;

      if ( copytimestep )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  pstreamDefTimestep(streamID2, tsID2++);

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);
	      pstreamDefRecord(streamID2, varID, levelID);
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
	}
       
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( nfound == 0 )
    cdoWarning("Day %d%s not found!", dday, cmons[dmon]);

  if ( ! lcopy )
    if ( array ) Free(array);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
