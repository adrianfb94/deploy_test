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

      Timselstat    timselrange        Time selection range
      Timselstat    timselmin          Time selection minimum
      Timselstat    timselmax          Time selection maximum
      Timselstat    timselsum          Time selection sum
      Timselstat    timselmean         Time selection mean
      Timselstat    timselavg          Time selection average
      Timselstat    timselvar          Time selection variance
      Timselstat    timselvar1         Time selection variance [Normalize by (n-1)]
      Timselstat    timselstd          Time selection standard deviation
      Timselstat    timselstd1         Time selection standard deviation [Normalize by (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Timselstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int nrecs = 0;
  int varID, levelID;
  int tsID;
  int nsets;
  int nmiss;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("timselrange", func_range, 0, NULL);
  cdoOperatorAdd("timselmin",   func_min,   0, NULL);
  cdoOperatorAdd("timselmax",   func_max,   0, NULL);
  cdoOperatorAdd("timselsum",   func_sum,   0, NULL);
  cdoOperatorAdd("timselmean",  func_mean,  0, NULL);
  cdoOperatorAdd("timselavg",   func_avg,   0, NULL);
  cdoOperatorAdd("timselvar",   func_var,   0, NULL);
  cdoOperatorAdd("timselvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("timselstd",   func_std,   0, NULL);
  cdoOperatorAdd("timselstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;

  operatorInputArg("nsets <noffset <nskip>>");

  int nargc  = operatorArgc();
  int ndates = parameter2int(operatorArgv()[0]);
  int noffset = (nargc > 1) ? parameter2int(operatorArgv()[1]) : 0;
  int nskip   = (nargc > 2) ? parameter2int(operatorArgv()[2]) : 0;
  // clang-format on

  if ( cdoVerbose ) cdoPrint("nsets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
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

  for ( tsID = 0; tsID < noffset; tsID++ )
    {
      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
              recinfo[recID].varID   = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst  = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
	    }
	}
    }

  int otsID = 0;
  if ( tsID < noffset )
    {
      cdoWarning("noffset is larger than number of timesteps!");
      goto LABEL_END;
    }

  while ( TRUE )
    {
      for ( nsets = 0; nsets < ndates; nsets++ )
	{
	  nrecs = pstreamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;

	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);

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
                      pvars2->nmiss = (size_t)nmiss;
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

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      pstreamDefTimestep(streamID2, otsID);

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

      for ( int i = 0; i < nskip; i++ )
	{
	  nrecs = pstreamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      if ( nrecs == 0 ) break;
    }

 LABEL_END:


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
