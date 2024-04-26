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

      Yseasstat  yseasrange      Multi-year seasonal range
      Yseasstat  yseasmin        Multi-year seasonal minimum
      Yseasstat  yseasmax        Multi-year seasonal maximum
      Yseasstat  yseassum        Multi-year seasonal sum
      Yseasstat  yseasmean       Multi-year seasonal mean
      Yseasstat  yseasavg        Multi-year seasonal average
      Yseasstat  yseasvar        Multi-year seasonal variance
      Yseasstat  yseasvar1       Multi-year seasonal variance [Normalize by (n-1)]
      Yseasstat  yseasstd        Multi-year seasonal standard deviation
      Yseasstat  yseasstd1       Multi-year seasonal standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  NSEAS       4

typedef struct {
  int vdate;
  int vtime;
}
date_time_t;


void set_date(int vdate_new, int vtime_new, date_time_t *datetime)
{
  int year, month, day;
  cdiDecodeDate(vdate_new, &year, &month, &day);
  if ( month == 12 ) vdate_new = cdiEncodeDate(year-1, month, day);

  if ( vdate_new > datetime->vdate )
    {
      datetime->vdate = vdate_new;
      datetime->vtime = vtime_new;
    }
}


void *Yseasstat(void *argument)
{
  int varID;
  int year, month, day;
  int nrecs;
  int levelID;
  int seas_nsets[NSEAS];
  int nmiss;
  date_time_t datetime[NSEAS];
  field_type **vars1[NSEAS], **vars2[NSEAS], **samp1[NSEAS];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("yseasrange", func_range, 0, NULL);
  cdoOperatorAdd("yseasmin",   func_min,   0, NULL);
  cdoOperatorAdd("yseasmax",   func_max,   0, NULL);
  cdoOperatorAdd("yseassum",   func_sum,   0, NULL);
  cdoOperatorAdd("yseasmean",  func_mean,  0, NULL);
  cdoOperatorAdd("yseasavg",   func_avg,   0, NULL);
  cdoOperatorAdd("yseasvar",   func_var,   0, NULL);
  cdoOperatorAdd("yseasvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("yseasstd",   func_std,   0, NULL);
  cdoOperatorAdd("yseasstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  for ( int seas = 0; seas < NSEAS; seas++ )
    {
      vars1[seas]  = NULL;
      vars2[seas]  = NULL;
      samp1[seas]  = NULL;
      seas_nsets[seas]  = 0;
      datetime[seas].vdate = 0;
      datetime[seas].vtime = 0;
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<recinfo_type> recinfo(maxrecs);

  int gridsizemax = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsizemax*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);
      cdiDecodeDate(vdate, &year, &month, &day);

      int seas = month_to_season(month);

      set_date(vdate, vtime, &datetime[seas]);

      if ( vars1[seas] == NULL )
	{
	  vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[seas] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd || lrange )
	    vars2[seas] = field_malloc(vlistID1, FIELD_PTR);
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
              recinfo[recID].varID   = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst  = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
	    }

          field_type *psamp1 = &samp1[seas][varID][levelID];
          field_type *pvars1 = &vars1[seas][varID][levelID];
          field_type *pvars2 = vars2[seas] ? &vars2[seas][varID][levelID] : NULL;
          int nsets = seas_nsets[seas];

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

      if ( seas_nsets[seas] == 0 && lvarstd )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
	    if ( recinfo[recID].lconst ) continue;

            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[seas][varID][levelID];
            field_type *pvars2 = &vars2[seas][varID][levelID];

            farmoq(pvars2, *pvars1);
	  }

      seas_nsets[seas]++;
      tsID++;
    }

  for ( int seas = 0; seas < NSEAS; seas++ )
    if ( seas_nsets[seas] )
      {
        int nsets = seas_nsets[seas];
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
	    if ( recinfo[recID].lconst ) continue;

            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *psamp1 = &samp1[seas][varID][levelID];
            field_type *pvars1 = &vars1[seas][varID][levelID];
            field_type *pvars2 = vars2[seas] ? &vars2[seas][varID][levelID] : NULL;

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

	taxisDefVdate(taxisID2, datetime[seas].vdate);
	taxisDefVtime(taxisID2, datetime[seas].vtime);
	pstreamDefTimestep(streamID2, otsID);

	for ( int recID = 0; recID < maxrecs; recID++ )
	  {
            if ( otsID && recinfo[recID].lconst ) continue;

            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[seas][varID][levelID];

	    pstreamDefRecord(streamID2, varID, levelID);
	    pstreamWriteRecord(streamID2, pvars1->ptr, (int)pvars1->nmiss);
	  }

	otsID++;
      }

  for ( int seas = 0; seas < NSEAS; seas++ )
    {
      if ( vars1[seas] != NULL )
	{
	  field_free(vars1[seas], vlistID1);
	  field_free(samp1[seas], vlistID1);
	  if ( lvarstd ) Free(vars2[seas]);
	}
    }

  if ( field.ptr ) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
