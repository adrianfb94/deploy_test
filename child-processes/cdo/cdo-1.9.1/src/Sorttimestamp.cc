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

     Sorttimestamp    sorttimestamp         Sort all timesteps
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024

typedef struct
{
  int      index;
  double   datetime;
}
timeinfo_t;

static
int cmpdatetime(const void *s1, const void *s2)
{
  int cmp = 0;
  const timeinfo_t *x = (timeinfo_t *)s1;
  const timeinfo_t *y = (timeinfo_t *)s2;
  /*
  printf("%g %g  %d %d\n", x->datetime, y->datetime, x, y);
  */
  if      ( x->datetime < y->datetime ) cmp = -1;
  else if ( x->datetime > y->datetime ) cmp =  1;

  return cmp;
}


void *Sorttimestamp(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID;
  int tsID, lasttsID = -1;
  int nalloc = 0;
  int vlistID2 = -1, taxisID2 = -1;
  int nmiss;
  int nvars = 0, nlevel;
  int *vdate = NULL, *vtime = NULL;
  field_type ***vars = NULL;

  cdoInitialize(argument);

  int nfiles = cdoStreamCnt() - 1;

  int xtsID = 0;
  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      int streamID1 = pstreamOpenRead(cdoStreamName(fileID));

      int vlistID1 = pstreamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      if ( fileID == 0 )
	{
	  vlistID2 = vlistDuplicate(vlistID1);
	  taxisID2 = taxisDuplicate(taxisID1);
	  if ( taxisHasBounds(taxisID2) )
	    {
	      cdoWarning("Time bounds unsupported by this operator, removed!");
	      taxisDeleteBounds(taxisID2);
	    }
	}
      else
	{
	  vlistCompare(vlistID2, vlistID1, CMP_ALL);
	}

      nvars = vlistNvars(vlistID1);

      int tsID = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
	{
	  if ( xtsID >= nalloc )
	    {
	      nalloc += NALLOC_INC;
	      vdate = (int*) Realloc(vdate, nalloc*sizeof(int));
	      vtime = (int*) Realloc(vtime, nalloc*sizeof(int));
	      vars  = (field_type ***) Realloc(vars, nalloc*sizeof(field_type **));
	    }

	  vdate[xtsID] = taxisInqVdate(taxisID1);
	  vtime[xtsID] = taxisInqVtime(taxisID1);

	  vars[xtsID] = field_malloc(vlistID1, FIELD_NONE);

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      vars[xtsID][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
	      pstreamReadRecord(streamID1, vars[xtsID][varID][levelID].ptr, &nmiss);
	      vars[xtsID][varID][levelID].nmiss = nmiss;
	    }

	  tsID++;
	  xtsID++;
	}

      pstreamClose(streamID1);
    }

  int nts = xtsID;

  timeinfo_t *timeinfo= (timeinfo_t*) Malloc(nts*sizeof(timeinfo_t));

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      int calendar = taxisInqCalendar(taxisID2);
      int julday = date_to_julday(calendar, vdate[tsID]);
      int secofday = time_to_sec(vtime[tsID]);
      double vdatetime = julday + secofday / 86400.;
      timeinfo[tsID].index    = tsID;
      timeinfo[tsID].datetime = vdatetime;
    }
  

  qsort(timeinfo, nts, sizeof(timeinfo_t), cmpdatetime);  	      


  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID2 = 0;
  for ( tsID = 0; tsID < nts; tsID++ )
    {
      xtsID = timeinfo[tsID].index;

      if ( tsID > 0 )
	{
	  if ( IS_EQUAL(timeinfo[tsID].datetime, timeinfo[lasttsID].datetime) )
	    {
	      if ( cdoVerbose )
		{
		  char vdatestr[32], vtimestr[32];
		  date2str(vdate[xtsID], vdatestr, sizeof(vdatestr));
		  time2str(vtime[xtsID], vtimestr, sizeof(vtimestr));
		  cdoPrint("Timestep %4d %s %s already exists, skipped!", xtsID, vdatestr, vtimestr);
		}
	      continue;
	    }
	}

      lasttsID = tsID;

      taxisDefVdate(taxisID2, vdate[xtsID]);
      taxisDefVtime(taxisID2, vtime[xtsID]);
      pstreamDefTimestep(streamID2, tsID2++);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( vars[xtsID][varID][levelID].ptr )
		{
		  nmiss = vars[xtsID][varID][levelID].nmiss;
		  pstreamDefRecord(streamID2, varID, levelID);
		  pstreamWriteRecord(streamID2, vars[xtsID][varID][levelID].ptr, nmiss);
		}
	    }
	}

      field_free(vars[xtsID], vlistID2);      
    }

  if ( vars  ) Free(vars);
  if ( vdate ) Free(vdate);
  if ( vtime ) Free(vtime);

  pstreamClose(streamID2);

  cdoFinish();

  return 0;
}
