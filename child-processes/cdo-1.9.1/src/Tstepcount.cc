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

      Tstepcount  tstepcount  Count number of timesteps
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024


static
double tstepcount(long nts, double missval1, double *array1, double refval)
{
  long j;
  long n = 0;

  if ( DBL_IS_EQUAL(refval, missval1) ) return missval1;

  for ( j = 0; j < nts; j++ )
    {
      n++;
      if ( DBL_IS_EQUAL(array1[j], refval) ) break;  
    }

  if ( j == nts )
    return missval1;
  else
    return (double) n;
}


void *Tstepcount(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID;
  int nalloc = 0;
  int nmiss;
  int nlevel;
  int vdate = 0, vtime = 0;
  double missval;
  double refval = 0;
  field_type ***vars = NULL;
  typedef struct
  {
    double *array1;
  } memory_t;

  cdoInitialize(argument);

  if ( operatorArgc() == 1 ) refval = parameter2double(operatorArgv()[0]);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistDefVarUnits(vlistID2, varID, "steps");
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vars  = (field_type ***) Realloc(vars, nalloc*sizeof(field_type **));
	}

      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
	  pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  int nts = tsID;

  memory_t *mem = (memory_t*) Malloc(ompNumThreads*sizeof(memory_t));
  for ( int i = 0; i < ompNumThreads; i++ )
    {
      mem[i].array1 = (double*) Malloc(nts*sizeof(double));
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(gridsize,mem,vars,varID,levelID,nts,missval,refval) schedule(dynamic,1)
#endif
	  for ( int i = 0; i < gridsize; i++ )
	    {
	      int ompthID = cdo_omp_get_thread_num();

	      for ( int tsID = 0; tsID < nts; tsID++ )
		mem[ompthID].array1[tsID] = vars[tsID][varID][levelID].ptr[i];

	      double count = tstepcount(nts, missval, mem[ompthID].array1, refval);
	      vars[0][varID][levelID].ptr[i] = count;
	    }
	}
    }

  for ( int i = 0; i < ompNumThreads; i++ )
    Free(mem[i].array1);
  Free(mem);

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  pstreamDefTimestep(streamID2, 0);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID2, varID);
      missval  = vlistInqVarMissval(vlistID2, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  pstreamDefRecord(streamID2, varID, levelID);

	  nmiss = 0;
	  for ( int i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(vars[0][varID][levelID].ptr[i], missval) ) nmiss++;

	  pstreamWriteRecord(streamID2, vars[0][varID][levelID].ptr, nmiss);
	}
    }

  for ( tsID = 0; tsID < nts; tsID++ ) field_free(vars[tsID], vlistID1);

  if ( vars ) Free(vars);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
