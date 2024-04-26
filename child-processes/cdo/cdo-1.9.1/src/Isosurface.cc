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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


double intlin(double x, double y1, double x1, double y2, double x2);

static
void isosurface(double isoval, long nlev1, double *lev1, field_type *field3D, field_type *field2D)
{
  bool lmiss1, lmiss2;

  long gridsize = gridInqSize(field3D->grid);
  long nmiss    = field3D->nmiss;
  double missval  = field3D->missval;
  double *data3D   = field3D->ptr;
  double *data2D   = field2D->ptr;

  for ( long i = 0; i < gridsize; ++i )
    {
      data2D[i] = missval;

      for ( long k = 0; k < (nlev1-1); ++k )
	{
	  double val1 = data3D[k*gridsize+i];
	  double val2 = data3D[(k+1)*gridsize+i];

	  if ( nmiss > 0 )
	    {
	      lmiss1 = DBL_IS_EQUAL(val1, missval);
	      lmiss2 = DBL_IS_EQUAL(val2, missval);
	      if ( lmiss1 && lmiss2 ) continue;
	      if ( lmiss1 && IS_EQUAL(isoval, val2) ) data2D[i] = lev1[k+1];
	      if ( lmiss2 && IS_EQUAL(isoval, val1) ) data2D[i] = lev1[k]  ;
	      if ( lmiss1 || lmiss2 ) continue;
	    }

	  if ( (isoval >= val1 && isoval <= val2) || (isoval >= val2 && isoval <= val1) )
	    {
	      if ( IS_EQUAL(val1, val2) )
		data2D[i] = lev1[k];
	      else
		data2D[i] = intlin(isoval, lev1[k], val1, lev1[k+1], val2);

	      break;
	    }
	}
    }

  nmiss = 0;
  for ( long i = 0; i < gridsize; ++i )
    if ( DBL_IS_EQUAL(data2D[i], missval) ) nmiss++;

  field2D->missval = missval;
  field2D->nmiss   = nmiss;
}


void *Isosurface(void *argument)
{
  int nlevel = 0;
  int nrecs;
  int gridID;
  int i, offset;
  int varID, levelID;
  int nmiss;
  int zaxisID, zaxisID1 = -1;
  double missval;
  double *single;

  cdoInitialize(argument);

  operatorInputArg("isoval");

  operatorCheckArgc(1);

  double isoval = parameter2double(operatorArgv()[0]);

  if ( cdoVerbose ) cdoPrint("Isoval: %g", isoval);


  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nzaxis = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF )
	if ( nlevel > 1 )
	  {
	    zaxisID1 = zaxisID;
	    break;
	  }
    }

  if ( i == nzaxis ) cdoAbort("No processable variable found!");

  int nlev1 = nlevel;
  double *lev1  = (double*) Malloc(nlev1*sizeof(double));
  cdoZaxisInqLevels(zaxisID1, lev1);

  int zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID2, i, zaxisIDsfc);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  int nvars = vlistNvars(vlistID1);

  bool *liso = (bool*)     Malloc(nvars*sizeof(bool));
  int *vars  = (int*)     Malloc(nvars*sizeof(int));
  field_type *vars1 = (field_type*) Malloc(nvars*sizeof(field_type));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      liso[varID] = (zaxisID == zaxisID1);

      field_init(&vars1[varID]);
      vars1[varID].grid    = gridID;
      vars1[varID].zaxis   = zaxisID;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  vars[varID] = FALSE;
	  vars1[varID].nmiss = 0;
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vars1[varID].ptr + offset;
	  
	  pstreamReadRecord(streamID1, single, &nmiss);
	  vars1[varID].nmiss += nmiss;
	  vars[varID] = TRUE;
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      if ( liso )
		{
		  isosurface(isoval, nlev1, lev1, &vars1[varID], &field);

		  pstreamDefRecord(streamID2, varID, 0);
		  pstreamWriteRecord(streamID2, field.ptr, field.nmiss);
		}
	      else
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
		  missval  = vlistInqVarMissval(vlistID2, varID);

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      offset   = gridsize*levelID;
		      single   = vars1[varID].ptr + offset;

		      nmiss = 0;
		      for ( i = 0; i < gridsize; ++i )
			if ( DBL_IS_EQUAL(single[i], missval) ) nmiss++;

		      pstreamDefRecord(streamID2, varID, levelID);
		      pstreamWriteRecord(streamID2, single, nmiss);
		    }
		}
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ ) Free(vars1[varID].ptr);
  Free(vars1);

  Free(vars);
  Free(liso);
  if (lev1) Free(lev1);

  Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
