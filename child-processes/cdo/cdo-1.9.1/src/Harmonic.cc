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

      Harmonic   harmonic        Harmonic
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Harmonic(void *argument)
{
  int gridsize;
  int nrecs;
  int varID, levelID;
  int nmiss;
  int offset;
  int nlevel;
  int vdate = 0, vtime = 0;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  double missval;

  cdoInitialize(argument);

  operatorInputArg("wave number and wave length of first harmonic in number of timesteps");

  operatorCheckArgc(2);

  int n_out = parameter2int(operatorArgv()[0]);
  int n     = parameter2int(operatorArgv()[1]);

  if ( n_out > 9 ) cdoAbort("Maximum number of wave numbers is 9!");

  if ( n < 1 || n < 2 * n_out )
    cdoAbort("The wave length must be positive and smaller than "
	     "2 times the number of requested harmonics (=%d)!", n_out);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  int *streamIDs = (int*) Malloc(n_out*sizeof(int));

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  for ( int j = 0; j < n_out; ++j )
    {
      sprintf(filename+nchars, "%1d", j+1);
      if ( filesuffix[0] )
	sprintf(filename+nchars+1, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      int streamID2 = pstreamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamIDs[j] = streamID2;

      pstreamDefVlist(streamID2, vlistID2);
    }

  int nvars = vlistNvars(vlistID1);

  double ***out  = (double ***) Malloc(n_out*sizeof(double **));
  double ***work = (double ***) Malloc(2*n_out*sizeof(double **));

  for ( int j = 0; j < n_out; ++j )
    {
      out[j] = (double **) Malloc(nvars*sizeof(double *));
      for ( varID = 0; varID < nvars; ++varID )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  out[j][varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
	}
    }

  for ( int j = 0; j < n_out*2; ++j )
    {
      work[j] = (double **) Malloc(nvars*sizeof(double *));
      for ( varID = 0; varID < nvars; ++varID )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  work[j][varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
	  memset(work[j][varID], 0, gridsize*nlevel*sizeof(double));
	}
    }


  gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID == 0 )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array, &nmiss);

	  if ( nmiss > 0 ) cdoAbort("Missing values are not allowed!");

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;

	  for ( int j = 0; j < n_out; ++j )
	    {
	      double sine   = sin(2 * M_PI * (((j + 1) * (tsID+1)) % n) / n);
	      double cosine = cos(2 * M_PI * (((j + 1) * (tsID+1)) % n) / n);
	      for ( int i = 0; i < gridsize; i++ )
		{
		  work[j][varID][i+offset]         += array[i] * sine;
		  work[n_out + j][varID][i+offset] += array[i] * cosine;
		}
	    }
	}

      tsID++;
    }

  int nts = tsID;

  if ( array ) Free(array);

  pstreamClose(streamID1);

  if ( nts%n )
    {
      cdoAbort("The length of first harmonic (=%d)"
	       " does not divide the number of timesteps (=%d)!", n, nts);
    }

  for ( int j = 0; j < n_out && 2*(j+1) < n; j++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( int i = 0; i < gridsize; i++ )
		out[j][varID][i+offset] = sqrt(work[j][varID][i+offset] * work[j][varID][i+offset] +
					work[n_out+j][varID][i+offset] * work[n_out+j][varID][i+offset]) * 2 / nts;
	    }
	}
    }

  if ( 2*n_out == n )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( int i = 0; i < gridsize; i++ )
		out[n_out - 1][varID][i+offset] = work[2 * n_out - 1][varID][i+offset] / nts;
	    }
	}
    }

  int nout = n_out;

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  for ( int j = 0; j < nout; j++ )
    {
      int streamID2 = streamIDs[j];
      pstreamDefTimestep(streamID2, 0);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      pstreamDefRecord(streamID2, varID, levelID);
	      pstreamWriteRecord(streamID2, out[j][varID]+offset, 0);
	    }
	}
    }

  for ( int j = 0; j < n_out && 2 * (j + 1) < n; j++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  missval  = vlistInqVarMissval(vlistID2, varID);
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( int i = 0; i < gridsize; i++ )
		{
		  out[j][varID][i+offset] = work[j][varID][i+offset] || work[n_out+j][varID][i+offset]
		              ? atan2 (work[j][varID][i+offset], work[n_out+j][varID][i+offset]) *
		                n / (j + 1) / 2 / M_PI : missval;
	  
		  if ( out[j][varID][i+offset] < 0 )
		    out[j][varID][i+offset] += n / (j + 1.);
		}
	    }
	}
    }

  nout = n_out;
  if ( 2*n_out == n ) nout -= 1;

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  for ( int j = 0; j < nout; j++ )
    {
      int streamID2 = streamIDs[j];
      pstreamDefTimestep(streamID2, 1);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  missval  = vlistInqVarMissval(vlistID2, varID);
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      pstreamDefRecord(streamID2, varID, levelID);
	      nmiss = 0;
	      for ( int i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(out[j][varID][i+offset], missval) ) nmiss++;
	      pstreamWriteRecord(streamID2, out[j][varID]+offset, nmiss);
	    }
	}
    }

  for ( int j = 0; j < n_out; j++ )
    {
      int streamID2 = streamIDs[j];
      pstreamClose(streamID2);
    }

  Free(streamIDs);

  for ( int j = 0; j < n_out; ++j )
    {
      for ( varID = 0; varID < nvars; ++varID )
        Free(out[j][varID]);

      Free(out[j]);
    }

  Free(out);

  for ( int j = 0; j < n_out*2; ++j )
    {
      for ( varID = 0; varID < nvars; ++varID )
        Free(work[j][varID]);

      Free(work[j]);
    }

  Free(work);

  cdoFinish();

  return 0;
}
