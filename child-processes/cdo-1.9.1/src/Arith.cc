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

      Arith      add             Add two fields
      Arith      sub             Subtract two fields
      Arith      mul             Multiply two fields
      Arith      div             Divide two fields
      Arith      min             Minimum of two fields
      Arith      max             Maximum of two fields
      Arith      atan2           Arc tangent of two fields
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Arith(void *argument)
{
  enum {FILL_NONE, FILL_TS, FILL_VAR, FILL_VARTS, FILL_FILE};
  int filltype = FILL_NONE;
  int nmiss;
  int nrecs, nvars = 0;
  int nlevels2 = 1;
  int varID, levelID;
  int levelID2;
  int *varnmiss2 = NULL;
  int **varnmiss = NULL;
  double *vardata2 = NULL;
  double **vardata = NULL;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("add",     func_add,     0, NULL);
  cdoOperatorAdd("sub",     func_sub,     0, NULL);
  cdoOperatorAdd("mul",     func_mul,     0, NULL);
  cdoOperatorAdd("div",     func_div,     0, NULL);
  cdoOperatorAdd("min",     func_min,     0, NULL);
  cdoOperatorAdd("max",     func_max,     0, NULL);
  cdoOperatorAdd("atan2",   func_atan2,   0, NULL);
  cdoOperatorAdd("setmiss", func_setmiss, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  operatorCheckArgc(0);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int streamIDx1 = streamID1;
  int streamIDx2 = streamID2;

  field_type field1, field2;
  field_type *fieldx1 = &field1;
  field_type *fieldx2 = &field2;

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistIDx1 = vlistID1;
  int vlistIDx2 = vlistID2;

  if ( cdoVerbose ) vlistPrint(vlistID1);
  if ( cdoVerbose ) vlistPrint(vlistID2);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisIDx1 = taxisID1;

  int ntsteps1 = vlistNtsteps(vlistID1);
  int ntsteps2 = vlistNtsteps(vlistID2);
  if ( ntsteps1 == 0 ) ntsteps1 = 1;
  if ( ntsteps2 == 0 ) ntsteps2 = 1;

  bool lfill1, lfill2;
  if ( vlistNvars(vlistID1) == 1 && vlistNvars(vlistID2) == 1 )
    {
      lfill2 = vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1;
      lfill1 = vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1;
    }
  else
    {
      lfill2 = vlistNvars(vlistID1) != 1 && vlistNvars(vlistID2) == 1;
      lfill1 = vlistNvars(vlistID1) == 1 && vlistNvars(vlistID2) != 1;
    }

  if ( lfill2 )
    {
      nlevels2 = vlistCompareX(vlistID1, vlistID2, CMP_DIM);

      if ( ntsteps1 != 1 && ntsteps2 == 1 )
	{
	  filltype = FILL_VAR;
	  cdoPrint("Filling up stream2 >%s< by copying the first variable.", cdoStreamName(1)->args);
	}
      else
	{
	  filltype = FILL_VARTS;
	  cdoPrint("Filling up stream2 >%s< by copying the first variable of each timestep.", cdoStreamName(1)->args);
	}
    }
  else if ( lfill1 )
    {
      nlevels2 = vlistCompareX(vlistID2, vlistID1, CMP_DIM);

      if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_VAR;
	  cdoPrint("Filling up stream1 >%s< by copying the first variable.", cdoStreamName(0)->args);
	}
      else
	{
	  filltype = FILL_VARTS;
	  cdoPrint("Filling up stream1 >%s< by copying the first variable of each timestep.", cdoStreamName(0)->args);
	}
      streamIDx1 = streamID2;
      streamIDx2 = streamID1;
      vlistIDx1 = vlistID2;
      vlistIDx2 = vlistID1;
      taxisIDx1 = taxisID2;
      fieldx1 = &field2;
      fieldx2 = &field1;
    }

  if ( filltype == FILL_NONE ) vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int gridsize = vlistGridsizeMax(vlistIDx1);

  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));
  if ( filltype == FILL_VAR || filltype == FILL_VARTS )
    {
      vardata2 = (double*) Malloc(gridsize*nlevels2*sizeof(double));
      varnmiss2 = (int*) Malloc(nlevels2*sizeof(int));
    }

  if ( cdoVerbose ) cdoPrint("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if ( filltype == FILL_NONE )
    {
      if ( ntsteps1 != 1 && ntsteps2 == 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream2 >%s< by copying the first timestep.", cdoStreamName(1)->args);
	}
      else if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoStreamName(0)->args);
	  streamIDx1 = streamID2;
          streamIDx2 = streamID1;
	  vlistIDx1 = vlistID2;
	  vlistIDx2 = vlistID1;
	  taxisIDx1 = taxisID2;
	  fieldx1 = &field2;
	  fieldx2 = &field1;
	}

      if ( filltype == FILL_TS )
	{
	  nvars  = vlistNvars(vlistIDx2);
	  vardata  = (double **) Malloc(nvars*sizeof(double *));
	  varnmiss = (int **) Malloc(nvars*sizeof(int *));
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      int nlev     = zaxisInqSize(vlistInqVarZaxis(vlistIDx2, varID));
	      vardata[varID]  = (double*) Malloc(nlev*gridsize*sizeof(double));
	      varnmiss[varID] = (int*) Malloc(nlev*sizeof(int));
	    }
	}
    }

  int vlistID3 = vlistDuplicate(vlistIDx1);
  if ( filltype == FILL_TS && vlistIDx1 != vlistID1 )
    {
      nvars  = vlistNvars(vlistID1);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarMissval(vlistID3, varID, vlistInqVarMissval(vlistID1, varID));
    }

  int taxisID3 = taxisDuplicate(taxisIDx1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int tsID = 0;
  int tsID2 = 0;
  while ( (nrecs = pstreamInqTimestep(streamIDx1, tsID)) )
    {
      if ( tsID == 0 || filltype == FILL_NONE || filltype == FILL_FILE || filltype == FILL_VARTS )
	{
	  int nrecs2 = pstreamInqTimestep(streamIDx2, tsID2);
	  if ( nrecs2 == 0 )
	    {
	      if ( filltype == FILL_NONE && streamIDx2 == streamID2 )
		{
		  filltype = FILL_FILE;
		  cdoPrint("Filling up stream2 >%s< by copying all timesteps.", cdoStreamName(1)->args);
		}

	      if ( filltype == FILL_FILE )
		{
		  tsID2 = 0;
		  pstreamClose(streamID2);
		  streamID2 = pstreamOpenRead(cdoStreamName(1));
		  streamIDx2 = streamID2;

		  vlistID2 = pstreamInqVlist(streamID2);
		  vlistIDx2 = vlistID2;

		  vlistCompare(vlistID1, vlistID2, CMP_DIM);

		  nrecs2 = pstreamInqTimestep(streamIDx2, tsID2);
		  if ( nrecs2 == 0 )
		    cdoAbort("Empty input stream %s!", cdoStreamName(1)->args);
		}
	      else
		cdoAbort("Input streams have different number of timesteps!");
	    }
	}

      taxisCopyTimestep(taxisID3, taxisIDx1);

      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamIDx1, &varID, &levelID);
	  pstreamReadRecord(streamIDx1, fieldx1->ptr, &nmiss);
          fieldx1->nmiss = (size_t) nmiss;
          int varID2 = varID;
          
	  if ( tsID == 0 || filltype == FILL_NONE || filltype == FILL_FILE || filltype == FILL_VARTS )
	    {
	      bool lstatus = nlevels2 > 1 ? varID == 0 : recID == 0;

	      if ( lstatus || (filltype != FILL_VAR && filltype != FILL_VARTS) )
		{
		  pstreamInqRecord(streamIDx2, &varID2, &levelID2);
		  pstreamReadRecord(streamIDx2, fieldx2->ptr, &nmiss);
                  fieldx2->nmiss = (size_t) nmiss;
                  if ( varID   != varID2 ) cdoAbort("Internal error, varIDs of input streams differ!");
                  if ( levelID != levelID2 ) cdoAbort("Internal error, levelIDs of input streams differ!");
		}

	      if ( filltype == FILL_TS )
		{
		  int gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
		  int offset   = gridsize*levelID;
		  memcpy(vardata[varID]+offset, fieldx2->ptr, gridsize*sizeof(double));
		  varnmiss[varID][levelID] = fieldx2->nmiss;
		}
	      else if ( lstatus && (filltype == FILL_VAR || filltype == FILL_VARTS) )
		{
		  int gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
		  int offset   = gridsize*levelID2;
		  memcpy(vardata2+offset, fieldx2->ptr, gridsize*sizeof(double));
		  varnmiss2[levelID2] = fieldx2->nmiss;
		}
	    }
	  else if ( filltype == FILL_TS )
	    {
	      int gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID2));
	      int offset   = gridsize*levelID;
	      memcpy(fieldx2->ptr, vardata[varID]+offset, gridsize*sizeof(double));
	      fieldx2->nmiss = varnmiss[varID][levelID];
	    }

	  fieldx1->grid    = vlistInqVarGrid(vlistIDx1, varID);
	  fieldx1->missval = vlistInqVarMissval(vlistIDx1, varID);

	  if ( filltype == FILL_VAR || filltype == FILL_VARTS )
	    {
	      levelID2 = (nlevels2 > 1) ? levelID : 0;
	      int gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
	      int offset   = gridsize*levelID2;
	      memcpy(fieldx2->ptr, vardata2+offset, gridsize*sizeof(double));
	      fieldx2->nmiss   = varnmiss2[levelID2];
	      fieldx2->grid    = vlistInqVarGrid(vlistIDx2, 0);
	      fieldx2->missval = vlistInqVarMissval(vlistIDx2, 0);
	    }
	  else
	    {
	      fieldx2->grid    = vlistInqVarGrid(vlistIDx2, varID2);
	      fieldx2->missval = vlistInqVarMissval(vlistIDx2, varID2);
	    }

	  farfun(&field1, field2, operfunc);

	  pstreamDefRecord(streamID3, varID, levelID);
	  pstreamWriteRecord(streamID3, field1.ptr, (int)field1.nmiss);
	}

      tsID++;
      tsID2++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID3);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  Free(vardata[varID]);
	  Free(varnmiss[varID]);
	}

      Free(vardata);
      Free(varnmiss);
    }

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);
  if ( vardata2   ) Free(vardata2);
  if ( varnmiss2  ) Free(varnmiss2);

  cdoFinish();

  return 0;
}
