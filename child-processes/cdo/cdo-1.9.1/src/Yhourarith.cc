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

      Yhourarith  yhouradd         Add multi-year hourly time series
      Yhourarith  yhoursub         Subtract multi-year hourly time series
      Yhourarith  yhourmul         Multiply multi-year hourly time series
      Yhourarith  yhourdiv         Divide multi-year hourly time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_HOUR  9301  /* 31*12*25 + 1 */

static
int hour_of_year(int vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
      
  if ( month >= 1 && month <= 12 && day >= 1 && day <=31 && hour >= 0 && hour < 24 )
    houroy = ((month-1)*31 + day - 1)*25 + hour + 1;
  else
    houroy = 0;

  if ( houroy < 0 || houroy >= MAX_HOUR )
    {
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoAbort("Hour of year %d out of range (%s %s)!", houroy, vdatestr, vtimestr);
    }

  return houroy;
}


void *Yhourarith(void *argument)
{
  int nrecs, nlev;
  int varID, levelID;
  int offset;
  int vdate, vtime;
  int nmiss;
  int houroy;
  int **varnmiss2[MAX_HOUR];
  double **vardata2[MAX_HOUR];

  cdoInitialize(argument);

  cdoOperatorAdd("yhouradd", func_add, 0, NULL);
  cdoOperatorAdd("yhoursub", func_sub, 0, NULL);
  cdoOperatorAdd("yhourmul", func_mul, 0, NULL);
  cdoOperatorAdd("yhourdiv", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

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

  for ( houroy = 0; houroy < MAX_HOUR ; ++houroy ) vardata2[houroy] = NULL;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID2, tsID)) )
    {
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);

      houroy = hour_of_year(vdate, vtime);
      if ( vardata2[houroy] != NULL ) cdoAbort("Hour of year %d already allocatd!", houroy);

      vardata2[houroy]  = (double **) Malloc(nvars*sizeof(double *));
      varnmiss2[houroy] = (int **) Malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[houroy][varID]  = (double*) Malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[houroy][varID] = (int*) Malloc(nlev*sizeof(int));
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID2, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  pstreamReadRecord(streamID2, vardata2[houroy][varID]+offset, &nmiss);
	  varnmiss2[houroy][varID][levelID] = nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      houroy = hour_of_year(vdate, vtime);
      if ( vardata2[houroy] == NULL ) cdoAbort("Hour of year %d not found!", houroy);

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
	  memcpy(field2.ptr, vardata2[houroy][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss   = varnmiss2[houroy][varID][levelID];
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

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    if ( vardata2[houroy] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    Free(vardata2[houroy][varID]);
	    Free(varnmiss2[houroy][varID]);
	  }

        Free(vardata2[houroy]);
        Free(varnmiss2[houroy]);
      }

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);

  cdoFinish();

  return 0;
}
