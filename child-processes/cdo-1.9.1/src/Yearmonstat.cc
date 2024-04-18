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

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "calendar.h"
#include "pstream.h"


void *Yearmonstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID;
  int i;
  int dpm;
  int year0 = 0, month0 = 0;
  int year, month, day;
  int nmiss;
  int nlevel;
  long nsets;
  double dsets;
  char vdatestr[32], vtimestr[32];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("yearmonmean",  func_mean, 0, NULL);
  cdoOperatorAdd("yearmonavg",   func_avg,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int calendar = taxisInqCalendar(taxisID1);
  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, calendar);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_type **samp1 = field_malloc(vlistID1, FIELD_NONE);

  int tsID    = 0;
  int otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      dsets = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  vdate = dtlist_get_vdate(dtlist, nsets);
	  vtime = dtlist_get_vtime(dtlist, nsets);
	  cdiDecodeDate(vdate, &year, &month, &day);

	  if ( nsets == 0 ) year0 = year;

	  if ( year != year0 ) break;

	  if ( nsets > 0 && month == month0 )
	    {
	      date2str(vdate0, vdatestr, sizeof(vdatestr));
	      time2str(vtime0, vtimestr, sizeof(vtimestr));
	      cdoWarning("   last timestep: %s %s", vdatestr, vtimestr);
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      cdoWarning("current timestep: %s %s", vdatestr, vtimestr);
	      cdoAbort("Month does not change!");
	    }

	  dpm = days_per_month(calendar, year, month);

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  pstreamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
		  vars1[varID][levelID].nmiss = (size_t) nmiss;

		  farcmul(&vars1[varID][levelID], dpm);

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = dpm;
		    }
		}
	      else
		{
		  pstreamReadRecord(streamID1, field.ptr, &nmiss);
                  field.nmiss   = (size_t) nmiss;
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  farcmul(&field, dpm);

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = dsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] += dpm;
		    }

		  farfun(&vars1[varID][levelID], field, operfunc);
		}
	    }

	  month0 = month;
	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  dsets += dpm;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( samp1[varID][levelID].ptr == NULL )
		farcdiv(&vars1[varID][levelID], dsets);
	      else
		fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	    }
	}

      if ( cdoVerbose )
	{
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  time2str(vtime0, vtimestr, sizeof(vtimestr));
	  cdoPrint("%s %s  nsets = %d", vdatestr, vtimestr, nsets);
	}

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      pstreamDefTimestep(streamID2, otsID);

      for ( int recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, vars1[varID][levelID].ptr,  (int)vars1[varID][levelID].nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);

  if ( field.ptr ) Free(field.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  dtlist_delete(dtlist);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
