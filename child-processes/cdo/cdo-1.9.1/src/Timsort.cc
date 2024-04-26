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

     Timsort    timsort         Sort over the time
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024

static
int cmpdarray(const void *s1, const void *s2)
{
  int cmp = 0;
  const double *x = (double *)s1;
  const double *y = (double *)s2;

  if      ( *x < *y ) cmp = -1;
  else if ( *x > *y ) cmp =  1;

  return cmp;
}


void *Timsort(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID;
  int nalloc = 0;
  int nmiss;
  int nlevel;
  int *vdate = NULL, *vtime = NULL;
  field_type ***vars = NULL;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID1);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vdate = (int*) Realloc(vdate, nalloc*sizeof(int));
	  vtime = (int*) Realloc(vtime, nalloc*sizeof(int));
	  vars  = (field_type ***) Realloc(vars, nalloc*sizeof(field_type **));
	}

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

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

  double **sarray = (double **) Malloc(ompNumThreads*sizeof(double *));
  for ( int i = 0; i < ompNumThreads; i++ )
    sarray[i] = (double*) Malloc(nts*sizeof(double));

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;

      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(gridsize,nts,sarray,vars,varID,levelID)
#endif
	  for ( int i = 0; i < gridsize; i++ )
	    {
	      int ompthID = cdo_omp_get_thread_num();

	      for ( int tsID = 0; tsID < nts; tsID++ )
		sarray[ompthID][tsID] = vars[tsID][varID][levelID].ptr[i];

	      qsort(sarray[ompthID], nts, sizeof(double), cmpdarray);  	      

	      for ( int tsID = 0; tsID < nts; tsID++ )
		vars[tsID][varID][levelID].ptr[i] = sarray[ompthID][tsID];
	    }
	}
    }

  for ( int i = 0; i < ompNumThreads; i++ )
    if ( sarray[i] ) Free(sarray[i]);

  if ( sarray ) Free(sarray);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      taxisDefVdate(taxisID2, vdate[tsID]);
      taxisDefVtime(taxisID2, vtime[tsID]);
      pstreamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( vars[tsID][varID][levelID].ptr )
		{
		  nmiss = vars[tsID][varID][levelID].nmiss;
		  pstreamDefRecord(streamID2, varID, levelID);
		  pstreamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		}
	    }
	}

      field_free(vars[tsID], vlistID1);      
    }

  if ( vars  ) Free(vars);
  if ( vdate ) Free(vdate);
  if ( vtime ) Free(vtime);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
