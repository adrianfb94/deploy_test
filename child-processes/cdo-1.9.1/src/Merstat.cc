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

      Merstat    merrange        Meridional range
      Merstat    mermin          Meridional minimum
      Merstat    mermax          Meridional maximum
      Merstat    mersum          Meridional sum
      Merstat    mermean         Meridional mean
      Merstat    meravg          Meridional average
      Merstat    merstd          Meridional standard deviation
      Merstat    merstd          Meridional standard deviation [Normalize by (n-1)]
      Merstat    mervar          Meridional variance
      Merstat    mervar          Meridional variance [Normalize by (n-1)]
      Merstat    merpctl         Meridional percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


void *Merstat(void *argument)
{
  int gridID1, gridID2 = -1, lastgrid = -1;
  int wstatus = FALSE;
  int index;
  int nmiss;
  int nrecs;
  int varID, levelID;
  char varname[CDI_MAX_NAME];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("merrange", func_range, 0, NULL);
  cdoOperatorAdd("mermin",   func_min,   0, NULL);
  cdoOperatorAdd("mermax",   func_max,   0, NULL);
  cdoOperatorAdd("mersum",   func_sum,   0, NULL);
  cdoOperatorAdd("mermean",  func_meanw, 1, NULL);
  cdoOperatorAdd("meravg",   func_avgw,  1, NULL);
  cdoOperatorAdd("mervar",   func_varw,  1, NULL);
  cdoOperatorAdd("mervar1",  func_var1w, 1, NULL);
  cdoOperatorAdd("merstd",   func_stdw,  1, NULL);
  cdoOperatorAdd("merstd1",  func_std1w, 1, NULL);
  cdoOperatorAdd("merpctl",  func_pctl,  0, NULL);
  // clang-format on
  
  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  bool needWeights = cdoOperatorF2(operatorID) != 0;

  double pn = 0;
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2double(operatorArgv()[0]);
      percentile_check_number(pn);
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  index = 0;
  gridID1 = vlistGrid(vlistID1, index);

  if ( gridInqType(gridID1) == GRID_LONLAT   ||
       gridInqType(gridID1) == GRID_GAUSSIAN ||
       gridInqType(gridID1) == GRID_GENERIC )
    {
      gridID2 = gridToMeridional(gridID1);
    }
  else
    {
      cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));
    }

  vlistChangeGridIndex(vlistID2, index, gridID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  int nlonmax = gridInqXsize(gridID1); /* max nlon ? */
  int lim = vlistGridsizeMax(vlistID1);

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);

  field1.ptr    = (double*) Malloc(lim*sizeof(double));
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) Malloc(lim*sizeof(double));

  field2.ptr  = (double*) Malloc(nlonmax*sizeof(double));
  field2.grid = gridID2;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = (size_t) nmiss;
	  field1.grid = vlistInqVarGrid(vlistID1, varID);
	  if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      wstatus = gridWeights(field1.grid, field1.weight);
	    }
	  if ( wstatus != 0 && tsID == 0 && levelID == 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      cdoWarning("Using constant grid cell area weights for variable %s!", varname);
	    }
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operfunc == func_pctl )
	    merpctl(field1, & field2, pn);
	  else  
	    merfun(field1, &field2, operfunc);

	  pstreamDefRecord(streamID2, varID,  levelID);
	  pstreamWriteRecord(streamID2, field2.ptr, (int)field2.nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( field1.ptr )    Free(field1.ptr);
  if ( field1.weight ) Free(field1.weight);
  if ( field2.ptr )    Free(field2.ptr);

  cdoFinish();

  return 0;
}
