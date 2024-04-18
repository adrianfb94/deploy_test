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

      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_MON    12


void *Ymonarith(void *argument)
{
  enum {MONTHLY, SEASONAL};
  int nrecs, nlev;
  int varID, levelID;
  int offset;
  int nmiss;
  int vdate, year, mon, day;
  int **varnmiss2[MAX_MON];
  double **vardata2[MAX_MON];
  const char *seas_name[4];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("ymonadd",  func_add, MONTHLY, NULL);
  cdoOperatorAdd("ymonsub",  func_sub, MONTHLY, NULL);
  cdoOperatorAdd("ymonmul",  func_mul, MONTHLY, NULL);
  cdoOperatorAdd("ymondiv",  func_div, MONTHLY, NULL);
  cdoOperatorAdd("yseasadd", func_add, SEASONAL, NULL);
  cdoOperatorAdd("yseassub", func_sub, SEASONAL, NULL);
  cdoOperatorAdd("yseasmul", func_mul, SEASONAL, NULL);
  cdoOperatorAdd("yseasdiv", func_div, SEASONAL, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int opertype = cdoOperatorF2(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int nvars  = vlistNvars(vlistID2);

  if ( opertype == SEASONAL ) get_season_name(seas_name);

  for ( mon = 0; mon < MAX_MON ; mon++ ) vardata2[mon] = NULL;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID2, tsID)) )
    {
      vdate = taxisInqVdate(taxisID2);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if ( mon < 1 || mon > MAX_MON ) cdoAbort("Month %d out of range!", mon);
      mon--;

      if ( opertype == SEASONAL ) mon = month_to_season(mon+1);

      if ( vardata2[mon] != NULL )
	{
	  if ( opertype == SEASONAL )
	    cdoAbort("Season %s already allocatd!", seas_name[mon]);
	  else
	    cdoAbort("Month %d already allocatd!", mon);
	}

      vardata2[mon]  = (double **) Malloc(nvars*sizeof(double *));
      varnmiss2[mon] = (int **) Malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[mon][varID]  = (double*) Malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[mon][varID] = (int*) Malloc(nlev*sizeof(int));
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID2, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  pstreamReadRecord(streamID2, vardata2[mon][varID]+offset, &nmiss);
	  varnmiss2[mon][varID][levelID] = nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if ( mon < 1 || mon > MAX_MON ) cdoAbort("Month %d out of range!", mon);
      mon--;

      if ( opertype == SEASONAL ) mon = month_to_season(mon+1);

      if ( vardata2[mon] == NULL )
	{
	  if ( opertype == SEASONAL )
	    cdoAbort("Season %s not found!", seas_name[mon]);
	  else
	    cdoAbort("Month %d not found!", mon);
	}

      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = (size_t) nmiss;
	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
          
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  memcpy(field2.ptr, vardata2[mon][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss = varnmiss2[mon][varID][levelID];
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  field2.missval = vlistInqVarMissval(vlistID2, varID);

	  farfun(&field1, field2, operfunc);

	  pstreamDefRecord(streamID3, varID, levelID);
	  pstreamWriteRecord(streamID3, field1.ptr, (int)field1.nmiss);
	}
      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  for ( mon = 0; mon < MAX_MON ; mon++ ) 
    if ( vardata2[mon] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    Free(vardata2[mon][varID]);
	    Free(varnmiss2[mon][varID]);
	  }

        Free(vardata2[mon]);
        Free(varnmiss2[mon]);
      }

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);

  cdoFinish();

  return 0;
}
