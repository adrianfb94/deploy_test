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

      Seasstat   seasrange       Seasonal range
      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasvar1        Seasonal variance [Normalize by (n-1)]
      Seasstat   seasstd         Seasonal standard deviation
      Seasstat   seasstd1        Seasonal standard deviation [Normalize by (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Seasstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int vdate0 = 0, vtime0 = 0;
  int vdate1 = 0, vtime1 = 0;
  int nrecs;
  int varID, levelID;
  int year, month, day, seas0 = 0;
  int nmiss;
  int oldmon = 0;
  int nseason = 0;
  const char *seas_name[4];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("seasrange", func_range, 0, NULL);
  cdoOperatorAdd("seasmin",   func_min,   0, NULL);
  cdoOperatorAdd("seasmax",   func_max,   0, NULL);
  cdoOperatorAdd("seassum",   func_sum,   0, NULL);
  cdoOperatorAdd("seasmean",  func_mean,  0, NULL);
  cdoOperatorAdd("seasavg",   func_avg,   0, NULL);
  cdoOperatorAdd("seasvar",   func_var,   0, NULL);
  cdoOperatorAdd("seasvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("seasstd",   func_std,   0, NULL);
  cdoOperatorAdd("seasstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  int season_start = get_season_start();
  get_season_name(seas_name);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  if ( taxisInqType(taxisID2) == TAXIS_FORECAST ) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<recinfo_type> recinfo(maxrecs);

  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  int gridsizemax = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsizemax*sizeof(double));

  field_type **samp1 = field_malloc(vlistID1, FIELD_NONE);
  field_type **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_type **vars2 = NULL;
  if ( lvarstd || lrange ) vars2 = field_malloc(vlistID1, FIELD_PTR);

  int tsID  = 0;
  int otsID = 0;
  while ( TRUE )
    {
      long nsets = 0;
      bool newseas = false;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  int vdate = dtlist_get_vdate(dtlist, nsets);
	  int vtime = dtlist_get_vtime(dtlist, nsets);

	  cdiDecodeDate(vdate, &year, &month, &day);

	  int newmon = month;

	  if ( season_start == START_DEC && newmon == 12 ) newmon = 0;

          int seas = month_to_season(month);

	  if ( nsets == 0 )
	    {
	      nseason++;
	      vdate0 = vdate;
	      vtime0 = vtime;
	      seas0  = seas;
	      oldmon = newmon;
	    }

	  if ( newmon < oldmon ) newseas = true;

	  if ( (seas != seas0) || newseas ) break;

	  oldmon = newmon;

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
                  recinfo[recID].varID   = varID;
                  recinfo[recID].levelID = levelID;
                  recinfo[recID].lconst  = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
		}

              field_type *psamp1 = &samp1[varID][levelID];
              field_type *pvars1 = &vars1[varID][levelID];
              field_type *pvars2 = vars2 ? &vars2[varID][levelID] : NULL;

	      int gridsize = pvars1->size;

	      if ( nsets == 0 )
		{
		  pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
		  pvars1->nmiss = (size_t)nmiss;
                  if ( lrange )
                    {
                      pvars2->nmiss = pvars1->nmiss;
		      for ( int i = 0; i < gridsize; i++ )
                        pvars2->ptr[i] = pvars1->ptr[i];
                    }

		  if ( nmiss > 0 || psamp1->ptr )
		    {
		      if ( psamp1->ptr == NULL )
			psamp1->ptr = (double*) Malloc(gridsize*sizeof(double));

		      for ( int i = 0; i < gridsize; i++ )
                        psamp1->ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
		    }
		}
	      else
		{
		  pstreamReadRecord(streamID1, field.ptr, &nmiss);
                  field.nmiss   = (size_t)nmiss;
		  field.grid    = pvars1->grid;
		  field.missval = pvars1->missval;

		  if ( field.nmiss > 0 || psamp1->ptr )
		    {
		      if ( psamp1->ptr == NULL )
			{
			  psamp1->ptr = (double*) Malloc(gridsize*sizeof(double));
			  for ( int i = 0; i < gridsize; i++ )
			    psamp1->ptr[i] = nsets;
			}

		      for ( int i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], pvars1->missval) )
			  psamp1->ptr[i]++;
		    }

		  if ( lvarstd )
		    {
		      farsumq(pvars2, field);
		      farsum(pvars1, field);
		    }
                  else if ( lrange )
                    {
                      farmin(pvars2, field);
                      farmax(pvars1, field);
                    }
		  else
		    {
		      farfun(pvars1, field, operfunc);
		    }
		}
	    }

	  if ( nsets == 0 && lvarstd )
            for ( int recID = 0; recID < maxrecs; recID++ )
              {
                if ( recinfo[recID].lconst ) continue;

                int varID   = recinfo[recID].varID;
                int levelID = recinfo[recID].levelID;
                field_type *pvars1 = &vars1[varID][levelID];
                field_type *pvars2 = &vars2[varID][levelID];

                farmoq(pvars2, *pvars1);
	      }

	  vdate1 = vdate;
	  vtime1 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      for ( int recID = 0; recID < maxrecs; recID++ )
        {
          if ( recinfo[recID].lconst ) continue;

          int varID   = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          field_type *psamp1 = &samp1[varID][levelID];
          field_type *pvars1 = &vars1[varID][levelID];
          field_type *pvars2 = vars2 ? &vars2[varID][levelID] : NULL;

          if ( lmean )
            {
              if ( psamp1->ptr ) fardiv(pvars1, *psamp1);
              else               farcdiv(pvars1, (double)nsets);
            }
          else if ( lvarstd )
            {
              if ( psamp1->ptr )
                {
                  if ( lstd ) farstd(pvars1, *pvars2, *psamp1, divisor);
                  else        farvar(pvars1, *pvars2, *psamp1, divisor);
                }
              else
                {
                  if ( lstd ) farcstd(pvars1, *pvars2, nsets, divisor);
                  else        farcvar(pvars1, *pvars2, nsets, divisor);
                }
            }
          else if ( lrange )
            {
              farsub(pvars1, *pvars2);
            }
        }

      if ( cdoVerbose )
	{
	  char vdatestr0[32], vtimestr0[32];
	  char vdatestr1[32], vtimestr1[32];
	  date2str(vdate0, vdatestr0, sizeof(vdatestr0));
	  time2str(vtime0, vtimestr0, sizeof(vtimestr0));
	  date2str(vdate1, vdatestr1, sizeof(vdatestr1));
	  time2str(vtime1, vtimestr1, sizeof(vtimestr1));
	  cdoPrint("season: %3d %3s  start: %s %s  end: %s %s ntimesteps: %d", 
		   nseason, seas_name[seas0], vdatestr0, vtimestr0, vdatestr1, vtimestr1, nsets);
	}

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      pstreamDefTimestep(streamID2, otsID);

      if ( nsets < 3 )
	{
	  char vdatestr[32];
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  cdoWarning("Season %3d (%s) has only %d input time step%s!", 
		     otsID+1, vdatestr, nsets, nsets == 1 ? "" : "s");
	}

      for ( int recID = 0; recID < maxrecs; recID++ )
	{
          if ( otsID && recinfo[recID].lconst ) continue;

          int varID   = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          field_type *pvars1 = &vars1[varID][levelID];

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, pvars1->ptr, (int)pvars1->nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if ( lvarstd ) field_free(vars2, vlistID1);

  dtlist_delete(dtlist);

  if ( field.ptr ) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
