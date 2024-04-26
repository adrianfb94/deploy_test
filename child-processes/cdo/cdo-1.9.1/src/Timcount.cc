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

      Timcount    timcount          Time counts
      Hourcount   hourcount         Hourly counts
      Daycount    daycount          Daily counts
      Moncount    moncount          Monthly counts
      Yearcount   yearcount         Yearly counts
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Timcount(void *argument)
{
  char indate1[DATE_LEN+1], indate2[DATE_LEN+1];
  int vdate0 = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID;
  int nmiss;
  int nwpv; // number of words per value; real:1  complex:2

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("timcount",  0, 31, NULL);
  cdoOperatorAdd("yearcount", 0, 10, NULL);
  cdoOperatorAdd("moncount",  0,  8, NULL);
  cdoOperatorAdd("daycount",  0,  6, NULL);
  cdoOperatorAdd("hourcount", 0,  4, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
      vlistDefVarUnits(vlistID2, varID, "No.");

  if ( cdoOperatorF2(operatorID) == 16 ) vlistDefNtsteps(vlistID2, 1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nrecords = vlistNrecs(vlistID1);
  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type **vars1 = field_malloc(vlistID1, FIELD_PTR);

  int tsID    = 0;
  int otsID   = 0;
  while ( TRUE )
    {
      int nsets = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
	{
	  int vdate = taxisInqVdate(taxisID1);
	  int vtime = taxisInqVtime(taxisID1);

	  if ( nsets == 0 ) SET_DATE(indate2, vdate, vtime);
	  SET_DATE(indate1, vdate, vtime);

	  if ( DATE_IS_NEQ(indate1, indate2, cmplen) ) break;

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      nwpv     = vars1[varID][levelID].nwpv;
	      gridsize = gridInqSize(vars1[varID][levelID].grid);

	      if ( nsets == 0 )
		{
		  for ( int i = 0; i < nwpv*gridsize; i++ )
		    vars1[varID][levelID].ptr[i] = vars1[varID][levelID].missval;
		  vars1[varID][levelID].nmiss = gridsize;
		}

              pstreamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t)nmiss;
              field.grid    = vars1[varID][levelID].grid;
	      field.missval = vars1[varID][levelID].missval;

              farcount(&vars1[varID][levelID], field);
	    }

	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      taxisDefVdate(taxisID2, vdate0);
      taxisDefVtime(taxisID2, vtime0);
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

  if ( field.ptr ) Free(field.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
