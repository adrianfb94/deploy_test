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

        Timstat2        timcor      correlates two data files on the same grid
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

// correlation in time
static
size_t correlation_t(size_t gridsize, double missval1, double missval2, size_t *nofvals, 
                     double *work0, double *work1, double *work2, double *work3, double *work4)
{
  size_t nmiss = 0;
  double temp0, temp1, temp2, temp3, temp4, temp5, temp6;
  double cor;

  for ( size_t i = 0; i < gridsize; ++i )
    {	  
      size_t nvals = nofvals[i];

      if ( nvals > 0 )
	{
	  temp0 = MULMN(work0[i], work1[i]);
	  temp1 = SUBMN(work4[i], DIVMN(temp0, nvals));
	  temp2 = MULMN(work0[i], work0[i]);
	  temp3 = MULMN(work1[i], work1[i]);
	  temp4 = SUBMN(work2[i], DIVMN(temp2, nvals));
	  temp5 = SUBMN(work3[i], DIVMN(temp3, nvals));
	  temp6 = MULMN(temp4, temp5);

	  cor = DIVMN(temp1, SQRTMN(temp6));

          if      ( cor < -1 )  cor = -1;
          else if ( cor >  1 )  cor =  1;

	  if ( DBL_IS_EQUAL(cor, missval1) ) nmiss++;
	}
      else
	{
	  nmiss++;
	  cor = missval1;
	}

      work0[i] = cor;
    }

  return nmiss;
}

// covariance in time
static
size_t covariance_t(size_t gridsize, double missval1, double missval2, size_t *nofvals, 
                    double *work0, double *work1, double *work2)
{
  size_t nmiss = 0;
  double covar;

  for ( size_t i = 0; i < gridsize; ++i )
    {	  
      size_t nvals = nofvals[i];

      if ( nvals > 0 )
	{
          double dnvals = nvals;
	  double temp = DIVMN( MULMN(work0[i], work1[i]), dnvals*dnvals);
	  covar = SUBMN( DIVMN(work2[i], dnvals), temp);

	  if ( DBL_IS_EQUAL(covar, missval1) ) nmiss++;
	}
      else
	{
	  nmiss++;
	  covar = missval1;
	}

      work0[i] = covar;
    }

  return nmiss;
}


void *Timstat2(void *argument)
{
  int vdate = 0, vtime = 0;
  int nrecs2, nlevs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("timcor",   func_cor,   0, NULL);
  cdoOperatorAdd("timcovar", func_covar, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  int nwork = 0;
  if      ( operfunc == func_cor   ) nwork = 5;
  else if ( operfunc == func_covar ) nwork = 3;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
 
  int nvars  = vlistNvars(vlistID1);
  int nrecs  = vlistNrecs(vlistID1);
  int nrecs3 = nrecs;
  int *recVarID   = (int*) Malloc(nrecs*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecs*sizeof(int));

  int taxisID1 = vlistInqTaxis(vlistID1);
  //int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);
 
  size_t gridsize = vlistGridsizeMax(vlistID1);

  double *array1  = (double*) Malloc(gridsize*sizeof(double));
  double *array2  = (double*) Malloc(gridsize*sizeof(double));
  				 
  double ****work = (double ****) Malloc(nvars*sizeof(double ***));
  size_t ***nofvals = (size_t ***) Malloc(nvars*sizeof(size_t **));

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID1, 0);  
      gridsize = gridInqSize(gridID);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      work[varID]    = (double ***) Malloc(nlevs*sizeof(double **));
      nofvals[varID] = (size_t **) Malloc(nlevs*sizeof(size_t *));  

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  nofvals[varID][levelID] = (size_t*) Malloc(gridsize*sizeof(size_t));
	  memset(nofvals[varID][levelID], 0, gridsize*sizeof(size_t));
      
	  work[varID][levelID] = (double **) Malloc(nwork*sizeof(double *));
	  for ( int i = 0; i < nwork; i++ )
	    {
	      work[varID][levelID][i] = (double*) Malloc(gridsize*sizeof(double));
	      memset(work[varID][levelID][i], 0, gridsize*sizeof(double));
	    }
	}
    }
 
  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      nrecs2 = pstreamInqTimestep(streamID2, tsID);
      if ( nrecs != nrecs2 )
        cdoWarning("Input streams have different number of records!");

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamInqRecord(streamID2, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;	     	     
	    }	 

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  double missval1 = vlistInqVarMissval(vlistID1, varID);
	  double missval2 = vlistInqVarMissval(vlistID2, varID);

	  pstreamReadRecord(streamID1, array1, &nmiss);
	  pstreamReadRecord(streamID2, array2, &nmiss);

	  if ( operfunc == func_cor )
	    {
	      for ( size_t i = 0; i < gridsize; ++i )
		{
		  if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
		       ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
		    {
		      work[varID][levelID][0][i] += array1[i];
		      work[varID][levelID][1][i] += array2[i];
		      work[varID][levelID][2][i] += array1[i]*array1[i];
		      work[varID][levelID][3][i] += array2[i]*array2[i];
		      work[varID][levelID][4][i] += array1[i]*array2[i];
		      nofvals[varID][levelID][i]++;
		    }
		}	 
	    }
	  else if ( operfunc == func_covar )
	    {
	      for ( size_t i = 0; i < gridsize; ++i )
		{
		  if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
		       ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
		    {
		      work[varID][levelID][0][i] += array1[i];
		      work[varID][levelID][1][i] += array2[i];
		      work[varID][levelID][2][i] += array1[i]*array2[i];
		      nofvals[varID][levelID][i]++;
		    }
		}	 
	    }
	}

      tsID++;
    }

  tsID = 0;
  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  pstreamDefTimestep(streamID3, tsID);

  for ( int recID = 0; recID < nrecs3; recID++ )
    {
      varID    = recVarID[recID];
      levelID  = recLevelID[recID];
   
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

      double missval1 = vlistInqVarMissval(vlistID1, varID);
      double missval2 = vlistInqVarMissval(vlistID2, varID);

      if ( operfunc == func_cor )
	{
	  nmiss = correlation_t(gridsize, missval1, missval2, nofvals[varID][levelID],
				work[varID][levelID][0], work[varID][levelID][1],
				work[varID][levelID][2], work[varID][levelID][3], 
				work[varID][levelID][4]);
	}
      else if ( operfunc == func_covar )
	{
	  nmiss = covariance_t(gridsize, missval1, missval2, nofvals[varID][levelID],
			       work[varID][levelID][0], work[varID][levelID][1],
			       work[varID][levelID][2]);
	}

      pstreamDefRecord(streamID3, varID, levelID);
      pstreamWriteRecord(streamID3, work[varID][levelID][0], nmiss);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  Free(nofvals[varID][levelID]);
	  for ( int i = 0; i < nwork; i++ )
	    Free(work[varID][levelID][i]);
	  Free(work[varID][levelID]);
	}
    
      Free(nofvals[varID]);
      Free(work[varID]);
    }
    
  Free(nofvals);
  Free(work);

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array1 )     Free(array1);
  if ( array2 )     Free(array2);
  if ( recVarID )   Free(recVarID);
  if ( recLevelID ) Free(recLevelID);
    
  cdoFinish();   
 
  return 0;
}
