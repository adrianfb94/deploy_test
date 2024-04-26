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

        Timstat3        varquot2test
        Timstat3        meandiff2test
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"


#define  NIN     2
#define  NOUT    1
#define  NFWORK  4
#define  NIWORK  2

void *Timstat3(void *argument)
{
  int streamID[NIN];
  int vlistID[NIN], vlistID2 = -1;
  int vdate = 0, vtime = 0;
  int nlevs;
  int is;
  int varID, levelID, gridID;
  int nmiss;
  double missval, missval1, missval2;
  double fractil_1, fractil_2, statistic;
  int ***iwork[NIWORK];
  field_type **fwork[NFWORK];
  field_type in[NIN], out[NOUT];
  int reached_eof[NIN];
  int n_in = NIN;

  cdoInitialize(argument);

  // clang-format off
  int VARQUOT2TEST  = cdoOperatorAdd("varquot2test",  0, 0, NULL);
  int MEANDIFF2TEST = cdoOperatorAdd("meandiff2test", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg("constant and risk (e.g. 0.05)");
  operatorCheckArgc(2);
  double rconst = parameter2double(operatorArgv()[0]);
  double risk   = parameter2double(operatorArgv()[1]);

  if ( operatorID == VARQUOT2TEST )
    {
      if ( rconst <= 0 )
	cdoAbort("Constant must be positive!");
      
      if ( risk <= 0 || risk >= 1 )
	cdoAbort("Risk must be greater than 0 and lower than 1!");
    }

  for ( is = 0; is < NIN; ++is )
    {
      streamID[is] = pstreamOpenRead(cdoStreamName(is));

      vlistID[is] = pstreamInqVlist(streamID[is]);
      if ( is > 0 )
	{
	  vlistID2 = pstreamInqVlist(streamID[is]);
	  vlistCompare(vlistID[0], vlistID2, CMP_ALL);
	}
    }

  int vlistID3 = vlistDuplicate(vlistID[0]);

  int gridsize = vlistGridsizeMax(vlistID[0]);
  int nvars = vlistNvars(vlistID[0]);
  int nrecs = vlistNrecs(vlistID[0]);
  int nrecs3 = nrecs;
  int *recVarID   = (int*) Malloc(nrecs*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecs*sizeof(int));

  int taxisID1 = vlistInqTaxis(vlistID[0]);
  int taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  for ( int i = 0; i < NIN; ++i ) reached_eof[i] = 0;

  for ( int i = 0; i < NIN; ++i )
    {
      field_init(&in[i]);
      in[i].ptr = (double*) Malloc(gridsize*sizeof(double));
    }
				 
  for ( int i = 0; i < NOUT; ++i )
    {
      field_init(&out[i]);
      out[i].ptr = (double*) Malloc(gridsize*sizeof(double));
    }
				 
  for ( int iw = 0; iw < NFWORK; ++iw )
    fwork[iw] = (field_type **) Malloc(nvars*sizeof(field_type *));
  for ( int iw = 0; iw < NIWORK; ++iw )
    iwork[iw] = (int ***) Malloc(nvars*sizeof(int **));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID   = vlistInqVarGrid(vlistID[0], varID);      
      gridsize = vlistGridsizeMax(vlistID[0]);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID[0], varID));
      missval  = missval1 = vlistInqVarMissval(vlistID[0], varID);
      // missval2 = vlistInqVarMissval(vlistID[1], varID); 

      for ( int iw = 0; iw < NFWORK; ++iw )
	fwork[iw][varID] = (field_type*) Malloc(nlevs*sizeof(field_type));
      for ( int iw = 0; iw < NIWORK; ++iw )
	iwork[iw][varID] = (int **) Malloc(nlevs*sizeof(int *));  

      for ( levelID = 0; levelID < nlevs; ++levelID )
	{
	  for ( int iw = 0; iw < NFWORK; ++iw )
	    {
	      field_init(&fwork[iw][varID][levelID]);
	      fwork[iw][varID][levelID].grid    = gridID;
	      fwork[iw][varID][levelID].nmiss   = 0;
	      fwork[iw][varID][levelID].missval = missval;
	      fwork[iw][varID][levelID].ptr     = (double*) Malloc(gridsize*sizeof(double));
	      memset(fwork[iw][varID][levelID].ptr, 0, gridsize*sizeof(double));
	    }

	  for ( int iw = 0; iw < NIWORK; ++iw )
	    {
	      iwork[iw][varID][levelID] = (int*) Malloc(gridsize*sizeof(int));
	      memset(iwork[iw][varID][levelID], 0, gridsize*sizeof(int));
	    }
	}
    }
 
  int tsID = 0;
  while ( TRUE )
    {
      for ( is = 0; is < NIN; ++is )
	{
	  if ( reached_eof[is] ) continue;

	  nrecs = pstreamInqTimestep(streamID[is], tsID);
	  if ( nrecs == 0 )
	    {
	      reached_eof[is] = 1;
	      continue;
	    }

	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID[is], &varID, &levelID);

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID[is], varID));

	      in[is].missval = vlistInqVarMissval(vlistID[is], varID);

	      if ( tsID == 0 && is == 0 )
		{
		  recVarID[recID] = varID;
		  recLevelID[recID] = levelID;	     	     
		}	 

	      pstreamReadRecord(streamID[is], in[is].ptr, &nmiss);
	      in[is].nmiss = (size_t) nmiss;
              
	      for ( int i = 0; i < gridsize; ++i )
		{
		  /*
		  if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
		  ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
		  */
		    {
		      fwork[NIN*is + 0][varID][levelID].ptr[i] += in[is].ptr[i];
		      fwork[NIN*is + 1][varID][levelID].ptr[i] += in[is].ptr[i] * in[is].ptr[i];
		      iwork[is][varID][levelID][i]++;
		    }
		}	 
	    }
	}

      for ( is = 0; is < NIN; ++is )
	if ( ! reached_eof[is] ) break;

      if ( is == NIN ) break;

      tsID++;
    }


  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  pstreamDefTimestep(streamID3, 0);

  for ( int recID = 0; recID < nrecs3; recID++ )
    {
      varID    = recVarID[recID];
      levelID  = recLevelID[recID];
  
      missval1 = fwork[0][varID][levelID].missval;
      missval2 = missval1;

      if ( operatorID == VARQUOT2TEST )
	{
	  for ( int i = 0; i < gridsize; ++i )
	    {
	      double fnvals0 = iwork[0][varID][levelID][i];
	      double fnvals1 = iwork[1][varID][levelID][i];

	      double temp0 = DIVMN( MULMN(fwork[0][varID][levelID].ptr[i], fwork[0][varID][levelID].ptr[i]), fnvals0);
	      double temp1 = DIVMN( MULMN(fwork[2][varID][levelID].ptr[i], fwork[2][varID][levelID].ptr[i]), fnvals1);
	      double temp2 = SUBMN(fwork[1][varID][levelID].ptr[i], temp0);
	      double temp3 = SUBMN(fwork[3][varID][levelID].ptr[i], temp1);
	      statistic = DIVMN(temp2, ADDMN(temp2, MULMN(rconst, temp3)));
	      
	      if ( fnvals0 <= 1 || fnvals1 <= 1 )
		fractil_1 = fractil_2 = missval1;
	      else
		beta_distr_constants((fnvals0 - 1) / 2,
				     (fnvals1 - 1) / 2, 1 - risk,
				     &fractil_1, &fractil_2, __func__);
	      out[0].ptr[i] = DBL_IS_EQUAL(statistic, missval1) ? missval1 : 
	       	              statistic <= fractil_1 || statistic >= fractil_2 ? 1 : 0;
	    }
	}									       
      else if ( operatorID == MEANDIFF2TEST )
	{
	  double fractil;
	  double mean_factor[NIN], var_factor[NIN];
	  
	  mean_factor[0] = 1;
	  mean_factor[1] = -1;
	  var_factor[0] = var_factor[1] = 1;

	  for ( int i = 0; i < gridsize; ++i )
	    {
	      double temp0 = 0;
	      double deg_of_freedom = -n_in;
	      for ( int j = 0; j < n_in; j++ )
		{
		  double fnvals = iwork[j][varID][levelID][i];
		  double tmp   = DIVMN( MULMN(fwork[2*j][varID][levelID].ptr[i], fwork[2*j][varID][levelID].ptr[i]), fnvals);
		  temp0 = ADDMN(temp0, DIVMN( SUBMN(fwork[2*j+1][varID][levelID].ptr[i], tmp), var_factor[j]));
		  deg_of_freedom = ADDMN(deg_of_freedom, fnvals);
		}

	      if ( !DBL_IS_EQUAL(temp0, missval1) && temp0 < 0 ) /* This is possible because */
		temp0 = 0;	                                 /* of rounding errors       */

	      double stddev_estimator = SQRTMN( DIVMN(temp0, deg_of_freedom));
	      double mean_estimator = -rconst;
	      for ( int j = 0; j < n_in; j++ )
		{
		  double fnvals = iwork[j][varID][levelID][i];
		  mean_estimator = ADDMN(mean_estimator,
				       MULMN(mean_factor[j],
					   DIVMN(fwork[2*j][varID][levelID].ptr[i], fnvals)));
		}

	      double temp1 = 0;
	      for ( int j = 0; j < n_in; j++ )
		{
		  double fnvals = iwork[j][varID][levelID][i];
		  temp1 = ADDMN(temp1, DIVMN( MULMN( MULMN(mean_factor[j], mean_factor[j]),
					     var_factor[j]), fnvals));
		}

	      double norm = SQRTMN(temp1);
	      
	      double temp2 = DIVMN( DIVMN(mean_estimator, norm), stddev_estimator);
	      fractil = deg_of_freedom < 1 ? missval1 :
		student_t_inv (deg_of_freedom, 1 - risk/2, __func__);

	      out[0].ptr[i] = DBL_IS_EQUAL(temp2, missval1)|| DBL_IS_EQUAL(fractil, missval1) ? 
		              missval1 : fabs(temp2) >= fractil;
	    }
	}

      nmiss = 0;
      for ( int i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(out[0].ptr[i], missval1) ) nmiss++;

      pstreamDefRecord(streamID3, varID, levelID);
      pstreamWriteRecord(streamID3, out[0].ptr, nmiss);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID[0], varID));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  for ( int iw = 0; iw < NFWORK; ++iw )
	    Free(fwork[iw][varID][levelID].ptr);
	  for ( int iw = 0; iw < NIWORK; ++iw )
	    Free(iwork[iw][varID][levelID]);
	}
    
      for ( int iw = 0; iw < NFWORK; ++iw ) Free(fwork[iw][varID]);
      for ( int iw = 0; iw < NIWORK; ++iw ) Free(iwork[iw][varID]);
    }

  for ( int iw = 0; iw < NFWORK; iw++ ) Free(fwork[iw]);
  for ( int iw = 0; iw < NIWORK; iw++ ) Free(iwork[iw]);


  pstreamClose(streamID3);
  for ( is = 0; is < NIN; ++is )
    pstreamClose(streamID[is]);

  for ( int i = 0; i < NIN; ++i ) Free(in[i].ptr);
  for ( int i = 0; i < NOUT; ++i ) Free(out[i].ptr);

  Free(recVarID);
  Free(recLevelID);
    
  cdoFinish();   
 
  return 0;
}
