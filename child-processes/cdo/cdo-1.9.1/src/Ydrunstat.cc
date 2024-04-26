/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Ydrunstat    ydrunmin          Multi-year daily running minimum
      Ydrunstat    ydrunmax          Multi-year daily running maximum
      Ydrunstat    ydrunsum          Multi-year daily running sum
      Ydrunstat    ydrunmean         Multi-year daily running mean
      Ydrunstat    ydrunavg          Multi-year daily running average
      Ydrunstat    ydrunvar          Multi-year daily running variance
      Ydrunstat    ydrunvar1         Multi-year daily running variance [Normalize by (n-1)]
      Ydrunstat    ydrunstd          Multi-year daily running standard deviation
      Ydrunstat    ydrunstd1         Multi-year daily running standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "calendar.h"
#include "pstream.h"


#define NDAY 373


typedef struct {
  int       vdate[NDAY];
  int       vtime[NDAY];  
  field_type **vars1[NDAY]; 
  field_type **vars2[NDAY];
  int       nsets[NDAY];
  int       vlist;
}
YDAY_STATS;


static YDAY_STATS *ydstatCreate(int vlistID);
static void ydstatDestroy(YDAY_STATS *stats);
static void ydstatUpdate(YDAY_STATS *stats, int vdate, int vtime, 
                         field_type **vars1, field_type **vars2, int nsets, int operfunc);
static void ydstatFinalize(YDAY_STATS *stats, int operfunc);


void *Ydrunstat(void *argument)
{
  int varID;
  int nrecs;
  int levelID;
  int tsID;
  int inp, its;
  int nmiss;
    
  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("ydrunmin",   func_min,   0, NULL);
  cdoOperatorAdd("ydrunmax",   func_max,   0, NULL);
  cdoOperatorAdd("ydrunsum",   func_sum,   0, NULL);
  cdoOperatorAdd("ydrunmean",  func_mean,  0, NULL);
  cdoOperatorAdd("ydrunavg",   func_avg,   0, NULL);
  cdoOperatorAdd("ydrunvar",   func_var,   0, NULL);
  cdoOperatorAdd("ydrunvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("ydrunstd",   func_std,   0, NULL);
  cdoOperatorAdd("ydrunstd1",  func_std1,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  operatorInputArg("number of timesteps");
  int ndates = parameter2int(operatorArgv()[0]);

  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  
  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int calendar = taxisInqCalendar(taxisID1);
  int dpy      = calendar_dpy(calendar);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<recinfo_type> recinfo(maxrecs);

  cdo_datetime_t *datetime = (cdo_datetime_t*) Malloc((ndates+1)*sizeof(cdo_datetime_t));
  
  YDAY_STATS *stats = ydstatCreate(vlistID1);
  field_type ***vars1 = (field_type ***) Malloc((ndates+1)*sizeof(field_type **));
  field_type ***vars2 = NULL;
  if ( lvarstd )
    vars2 = (field_type ***) Malloc((ndates+1)*sizeof(field_type **));
  
  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( lvarstd )
	vars2[its] = field_malloc(vlistID1, FIELD_PTR);
    }
  
  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
	cdoAbort("File has less then %d timesteps!", ndates);

      datetime[tsID].date = taxisInqVdate(taxisID1);
      datetime[tsID].time = taxisInqVtime(taxisID1);
	
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
              recinfo[recID].varID   = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst  = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
	    }
	  
          field_type *pvars1 = &vars1[tsID][varID][levelID];
          field_type *pvars2 = (vars2 && vars2[tsID]) ? &vars2[tsID][varID][levelID] : NULL;

	  pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
	  pvars1->nmiss = nmiss;

	  if ( lvarstd )
	    {
	      farmoq(pvars2, *pvars1);
	      for ( int inp = 0; inp < tsID; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], *pvars1);
		  farsum(&vars1[inp][varID][levelID], *pvars1);
		}
	    }
	  else
	    {
	      for ( int inp = 0; inp < tsID; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], *pvars1, operfunc);
		}
	    }
	}
    }
  
  while ( TRUE )
    {
      datetime_avg(dpy, ndates, datetime);
      
      int vdate = datetime[ndates].date;
      int vtime = datetime[ndates].time;
      
      if ( lvarstd )   
        ydstatUpdate(stats, vdate, vtime, vars1[0], vars2[0], ndates, operfunc);
      else
        ydstatUpdate(stats, vdate, vtime, vars1[0], NULL, ndates, operfunc);
        
      datetime[ndates] = datetime[0];
      vars1[ndates] = vars1[0];
      if ( lvarstd )
        vars2[ndates] = vars2[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  datetime[inp] = datetime[inp+1];
	  vars1[inp] = vars1[inp+1];
	  if ( lvarstd )
	    vars2[inp] = vars2[inp+1];
	}

      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      datetime[ndates-1].date = taxisInqVdate(taxisID1);
      datetime[ndates-1].time = taxisInqVtime(taxisID1);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  
          field_type *pvars1 = &vars1[ndates-1][varID][levelID];
          field_type *pvars2 = (vars2 && vars2[ndates-1]) ? &vars2[ndates-1][varID][levelID] : NULL;

	  pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
	  pvars1->nmiss = nmiss;

	  if ( lvarstd )
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], *pvars1);
		  farsum(&vars1[inp][varID][levelID], *pvars1);
		}
	      farmoq(pvars2, *pvars1);
	    }
	  else
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], *pvars1, operfunc);
		}
	    }
	}

      tsID++;
    }

  /*
  // set the year to the minimum of years found on output timestep
  int outyear = 1e9;
  int year, month, day;
  for ( int dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	cdiDecodeDate(stats->vdate[dayoy], &year, &month, &day);
	if ( year < outyear ) outyear = year;
      }
  for ( int dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	cdiDecodeDate(stats->vdate[dayoy], &year, &month, &day);
        // printf("vdates[%d] = %d  nsets = %d\n", dayoy, stats->vdate[dayoy], stats->nsets[dayoy]);
	if ( year > outyear ) stats->vdate[dayoy] = cdiEncodeDate(outyear, month, day);
      }
  */
  ydstatFinalize(stats, operfunc);

  int otsID = 0;

  for ( int dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	taxisDefVdate(taxisID2, stats->vdate[dayoy]);
	taxisDefVtime(taxisID2, stats->vtime[dayoy]);
	pstreamDefTimestep(streamID2, otsID);

        for ( int recID = 0; recID < maxrecs; recID++ )
          {
	    if ( otsID && recinfo[recID].lconst ) continue;

            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &stats->vars1[dayoy][varID][levelID];

	    pstreamDefRecord(streamID2, varID, levelID);
	    pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);
	  }

	otsID++;
      }
  
  for ( its = 0; its < ndates; its++ )
    {
      field_free(vars1[its], vlistID1);
      if ( lvarstd ) field_free(vars2[its], vlistID1);
    }
  
  ydstatDestroy(stats);
  Free(vars1);
  if ( lvarstd ) Free(vars2);

  if ( datetime ) Free(datetime);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

static
YDAY_STATS *ydstatCreate(int vlistID)
{  
  YDAY_STATS *stats = (YDAY_STATS*) Malloc(sizeof(YDAY_STATS));
  
  for ( int dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      stats->vdate[dayoy] = 0;
      stats->vtime[dayoy] = 0;
      stats->vars1[dayoy] = NULL;
      stats->vars2[dayoy] = NULL;
      stats->nsets[dayoy] = 0;
    }

  stats->vlist = vlistID;
  
  return stats;
}

static
void ydstatDestroy(YDAY_STATS *stats)
{
  int varID, levelID, nlevels;
  
  if ( stats != NULL )
    {
      int nvars = vlistNvars(stats->vlist);
      
      for ( int dayoy = 0; dayoy < NDAY; dayoy++ )
        {
          if ( stats->vars1[dayoy] != NULL )
            {
              for ( varID = 0; varID < nvars; varID++ )
                {
              	  nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
              	  for ( levelID = 0; levelID < nlevels; levelID++ )
              	    Free(stats->vars1[dayoy][varID][levelID].ptr);
              	  Free(stats->vars1[dayoy][varID]);
                }
              Free(stats->vars1[dayoy]);
            }
          if ( stats->vars2[dayoy] != NULL )
            {
              for ( varID = 0; varID < nvars; varID++ )
                {
              	  nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
              	  for ( levelID = 0; levelID < nlevels; levelID++ )
              	    Free(stats->vars2[dayoy][varID][levelID].ptr);
              	  Free(stats->vars2[dayoy][varID]);
                }
              Free(stats->vars2[dayoy]);
            }
        }
      Free(stats);    
    }
}

static
void ydstatUpdate(YDAY_STATS *stats, int vdate, int vtime, 
		  field_type **vars1, field_type **vars2, int nsets, int operfunc)
{
  int varID, levelID, nlevels;
  int gridsize;
  int year, month, day, dayoy;

  bool lvarstd = vars2 != NULL;

  int nvars = vlistNvars(stats->vlist);

  cdiDecodeDate(vdate, &year, &month, &day);

  if ( month >= 1 && month <= 12 )
    dayoy = (month - 1) * 31 + day;
  else
    dayoy = 0;

  if ( dayoy < 0 || dayoy >= NDAY )
    cdoAbort("day %d out of range!", dayoy);

  stats->vdate[dayoy] = vdate;
  stats->vtime[dayoy] = vtime;

  if ( stats->vars1[dayoy] == NULL )
    {
      stats->vars1[dayoy] = field_malloc(stats->vlist, FIELD_PTR);
      if ( lvarstd )
	stats->vars2[dayoy] = field_malloc(stats->vlist, FIELD_PTR);
    }

  for ( varID = 0; varID  < nvars; varID++ )
    {
      if ( vlistInqVarTimetype(stats->vlist, varID) == TIME_CONSTANT ) continue;
        
      gridsize = gridInqSize(vlistInqVarGrid(stats->vlist, varID));
      nlevels  = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
          
      for ( levelID = 0; levelID < nlevels; levelID++ )
        {
	  if ( stats->nsets[dayoy] == 0 )
	    {
	      memcpy(stats->vars1[dayoy][varID][levelID].ptr, vars1[varID][levelID].ptr, gridsize * sizeof(double));
	      stats->vars1[dayoy][varID][levelID].nmiss = vars1[varID][levelID].nmiss;
	       
	      if ( lvarstd )
	        {
	          memcpy(stats->vars2[dayoy][varID][levelID].ptr, vars2[varID][levelID].ptr, gridsize * sizeof(double));
	          stats->vars2[dayoy][varID][levelID].nmiss = vars2[varID][levelID].nmiss;
	        }
	    }
	  else
	    {
	      if ( lvarstd )
	        {
		  farsum(&stats->vars1[dayoy][varID][levelID], vars1[varID][levelID]);
		  farsum(&stats->vars2[dayoy][varID][levelID], vars2[varID][levelID]);
		}
	      else
		{
	          farfun(&stats->vars1[dayoy][varID][levelID], vars1[varID][levelID], operfunc);
		}
	    }
        }
    }

  stats->nsets[dayoy] += nsets;
}

static
void ydstatFinalize(YDAY_STATS *stats, int operfunc)
{
  int varID, levelID, nlevels;
  int dayoy;
  int divisor = operfunc == func_std1 || operfunc == func_var1;

  int nvars = vlistNvars(stats->vlist);
  
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
      	switch ( operfunc )
      	  {
	    case func_avg:
	    case func_mean:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTimetype(stats->vlist, varID) == TIME_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcdiv(&stats->vars1[dayoy][varID][levelID], (double) stats->nsets[dayoy]);
	        }
	      break;
	      
	    case func_std:
	    case func_std1:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTimetype(stats->vlist, varID) == TIME_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcstd(&stats->vars1[dayoy][varID][levelID], stats->vars2[dayoy][varID][levelID],
                            stats->nsets[dayoy], divisor);
	        }
	      break;
	      
	    case func_var:
	    case func_var1:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTimetype(stats->vlist, varID) == TIME_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcvar(&stats->vars1[dayoy][varID][levelID], stats->vars2[dayoy][varID][levelID],
			    stats->nsets[dayoy], divisor);
	        }
	      break;
      	  }
      }
}
