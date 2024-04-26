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

     EOF3d        eof3d             3D-EOF in spatial or time space
     EOF3d        eof3dspatial      3D-EOF in spatial space
     EOF3d        eof3dtime         3D-EOF in time space
*/
/*
 * TODO: 
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "statistic.h"


enum T_EIGEN_MODE get_eigenmode(void);
enum T_WEIGHT_MODE get_weightmode(void);


// NO MISSING VALUE SUPPORT ADDED SO FAR

void *EOF3d(void * argument)
{
  enum {EOF3D_, EOF3D_TIME, EOF3D_SPATIAL};

  size_t temp_size = 0, npack = 0;
  int varID, levelID;
  bool missval_warning = false;
  int nmiss, ngrids;
  int n = 0;
  size_t nlevs = 0;
  int timer_cov = 0, timer_eig = 0;

  int calendar = CALENDAR_STANDARD;

  double sum_w;
  double **cov = NULL;                                /* TODO: covariance matrix / eigenvectors after solving */
  double *eigv;

  if ( cdoTimer )
    {
      timer_cov  = timer_new("Timeof cov");
      timer_eig  = timer_new("Timeof eig");
    }

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("eof3d",        EOF3D_,        0, NULL);
  cdoOperatorAdd("eof3dtime",    EOF3D_TIME,    0, NULL);
  cdoOperatorAdd("eof3dspatial", EOF3D_SPATIAL, 0, NULL);
  // clang-format on

  int operatorID  = cdoOperatorID();
  int operfunc    = cdoOperatorF1(operatorID);

  operatorInputArg("Number of eigen functions to write out");
  int n_eig = parameter2int(operatorArgv()[0]);

  enum T_EIGEN_MODE eigen_mode = get_eigenmode();
  enum T_WEIGHT_MODE weight_mode = get_weightmode();

  /*  eigenvalues */

  if ( operfunc == EOF3D_SPATIAL )
    cdoAbort("Operator not Implemented - use eof3d or eof3dtime instead");

  int streamID1  = pstreamOpenRead(cdoStreamName(0));
  int vlistID1   = pstreamInqVlist(streamID1);

  /* COUNT NUMBER OF TIMESTEPS if EOF3D_ or EOF3D_TIME */
  int nts = vlistNtsteps(vlistID1);
  if ( nts == -1 )
    {
      nts = 0;
      while ( pstreamInqTimestep(streamID1, nts) ) nts++;

      if ( cdoVerbose ) cdoPrint("Counted %i timeSteps", nts);

      pstreamClose(streamID1);

      streamID1 = pstreamOpenRead(cdoStreamName(0));
      vlistID1  = pstreamInqVlist(streamID1);
    }
  else
    if ( cdoVerbose ) cdoPrint("Found %i timeSteps", nts);

  int taxisID1  = vlistInqTaxis(vlistID1);

  /* reset the requested number of eigen-function to the maximum if neccessary */
  if ( n_eig > nts )
    {
      cdoWarning("Solving in time-space:");
      cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
      cdoWarning("Setting n_eig to %i.", nts);
      n_eig = nts;
    }

  n = nts;

  if ( cdoVerbose )  cdoPrint("counted %i timesteps",n);

  int nvars      = vlistNvars(vlistID1);
  int nrecs;

  int gridID1    = vlistInqVarGrid(vlistID1, 0);
  size_t gridsizemax  = vlistGridsizeMax(vlistID1);

  /* allocation of temporary fields and output structures */
  double *in     = (double *) Malloc(gridsizemax*sizeof(double));
  int **datacounts = (int **) Malloc(nvars*sizeof(int*));
  double ***datafields   = (double ***) Malloc(nvars*sizeof(double **));
  double ***eigenvectors = (double ***) Malloc(nvars*sizeof(double **));
  double ***eigenvalues  = (double ***) Malloc(nvars*sizeof(double **));

  size_t maxlevs = 0;
  for ( varID = 0; varID < nvars; ++varID )
    {
      size_t gridsize = vlistGridsizeMax(vlistID1);
      nlevs     = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size = gridsize * nlevs;
      double missval = vlistInqVarMissval(vlistID1, varID);

      if ( nlevs > maxlevs ) maxlevs = nlevs;

      datacounts[varID]   = (int*) Malloc(nlevs*sizeof(int));
      datafields[varID]   = (double **) Malloc(nts*sizeof(double *));

      for ( int tsID = 0; tsID < nts; tsID++ )
	{
	  datafields[varID][tsID] = (double *) Malloc(temp_size*sizeof(double));
	  for ( size_t i = 0; i < temp_size; ++i ) datafields[varID][tsID][i] = 0;
	}
      datacounts[varID] = (int *) Malloc(temp_size*sizeof(int));	      
      for( size_t i = 0; i < temp_size; i++) datacounts[varID][i] = 0;
      
      eigenvectors[varID] = (double **) Malloc(n_eig*sizeof(double *));
      eigenvalues[varID]  = (double **) Malloc(nts*sizeof(double *));

      for ( int i = 0; i < n; i++ )
	{
	  if ( i < n_eig )
	    {
	      eigenvectors[varID][i] = (double *) Malloc(temp_size*sizeof(double));
	      for ( size_t i2 = 0; i2 < temp_size; ++i2 )
		eigenvectors[varID][i][i2] = missval;
	    }
	  
	  eigenvalues[varID][i]    = (double *) Malloc(1*sizeof(double));
	  eigenvalues[varID][i][0] = missval;
	}
    }

  if ( cdoVerbose)
    cdoPrint("Allocated eigenvalue/eigenvector with nts=%i, n=%i, gridsize=%zu for processing in %s",
	     nts, n, gridsizemax, "time_space");
  
  double *weight = (double *) Malloc(maxlevs*gridsizemax*sizeof(double));
  for ( size_t i = 0; i < maxlevs*gridsizemax; ++i ) weight[i] = 1.;

  if ( weight_mode == WEIGHT_ON )
    {
      int wstatus = gridWeights(gridID1, weight);
      if ( wstatus != 0 )
	{
	  weight_mode = WEIGHT_OFF;
	  cdoWarning("Using constant grid cell area weights!");
	}
      else
        {
          for ( size_t k = 1; k < maxlevs; ++k )
            for ( size_t i = 0; i < gridsizemax; ++i )
              weight[k*gridsizemax+i] =  weight[i];
        }
    }

  int tsID = 0;

  /* read the data and create covariance matrices for each var & level */
  while ( TRUE )
    {
      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          double missval = vlistInqVarMissval(vlistID1, varID);

          pstreamReadRecord(streamID1, in, &nmiss);

	  size_t offset = gridsize * levelID;
	  for ( size_t i = 0; i < gridsize; ++i )
	    {
	      if ( ! DBL_IS_EQUAL(in[i], missval ) )
		{
		  datafields[varID][tsID][offset + i] = in[i];
		  datacounts[varID][offset + i]++;
		}
	      else
		{
                  if ( datacounts[varID][offset + i] != 0 ) cdoAbort("Missing values unsupported!");
		  if ( missval_warning == false )
		    {
		      // cdoWarning("Missing Value Support not checked for this Operator!");
		      // cdoWarning("Does not work with changing locations of missing values in time.");
		      missval_warning = true;
		    }
		  datafields[varID][tsID][i+offset] = 0;
		}
	    }
        }
      tsID++;
    }

  if ( cdoVerbose ) 
    cdoPrint("Read data for %i variables",nvars);
  
  size_t *pack = (size_t*) Malloc(temp_size*sizeof(size_t)); //TODO

  for ( varID = 0; varID < nvars; varID++ )
    {
      size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size = gridsize * nlevs;
      double missval = vlistInqVarMissval(vlistID1, varID);

      if ( cdoVerbose )
        {
          char vname[64];
          vlistInqVarName(vlistID1,varID,&vname[0]);
          cdoPrint("============================================================================");
          cdoPrint("Calculating covariance matrix and SVD for var%i (%s)",varID,vname);
        }

      npack = 0;    // TODO already set to 0

      if ( cdoTimer ) timer_start(timer_cov);
      
      for ( size_t i = 0; i < temp_size ; i++ )
	{
	  if ( datacounts[varID][i] > 1 )
	    {
	      pack[npack] = i;
	      npack++;
	    }
	}

      sum_w = 1;
      if ( weight_mode == WEIGHT_ON )
	{
	  sum_w = 0;
	  for ( size_t i = 0; i < npack; i++ )  sum_w += weight[pack[i]];
	}

      if ( npack < 1 )
        {
          char vname[64];
          vlistInqVarName(vlistID1,varID,&vname[0]);
          cdoWarning("Refusing to calculate EOF from a single time step for var%i (%s)", varID+1, &vname[0]);
          continue;
        }

	  
      cov = (double **) Malloc(nts*sizeof(double*));
      for ( int j1 = 0; j1 < nts; j1++)
	cov[j1] = (double *) Malloc(nts*sizeof(double));
      eigv = (double *) Malloc(n*sizeof(double));

      if ( cdoVerbose )
        {
          cdoPrint("varID %i allocated eigv and cov with nts=%i and n=%i", varID, nts, n);
          cdoPrint("   npack=%zu, nts=%i temp_size=%zu", npack, nts, temp_size);
        }


#if defined(_OPENMP)
#pragma omp parallel for default(shared) schedule(static,2000)
#endif 
      for ( int j1 = 0; j1 < nts; j1++ )
        {
          double *df1p = datafields[varID][j1];
          for ( int j2 = j1; j2 < nts; j2++ )
            {
              double *df2p = datafields[varID][j2];
              double sum = 0;
              for ( size_t i = 0; i < npack; i++ )
                sum += weight[pack[i]%gridsizemax]*df1p[pack[i]]*df2p[pack[i]];
              cov[j2][j1] = cov[j1][j2] = sum / sum_w / nts;
            }
        }
      
      if ( cdoVerbose ) cdoPrint("calculated cov-matrix");

      /* SOLVE THE EIGEN PROBLEM */
      if ( cdoTimer ) timer_stop(timer_cov);


      if ( cdoTimer ) timer_start(timer_eig);

      if ( cdoVerbose ) 
	cdoPrint("Processed correlation matrix for var %2i | npack: %zu", varID, n);

      if ( eigen_mode == JACOBI ) 
	parallel_eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);
      else 
	eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);
      /* NOW: cov contains the eigenvectors, eigv the eigenvalues */

      if ( cdoVerbose ) 
	cdoPrint("Processed SVD decomposition for var %i from %zu x %zu matrix",varID,n,n);

      for( int eofID = 0; eofID < n; eofID++ )
	eigenvalues[varID][eofID][0] = eigv[eofID];
      
      if ( cdoTimer ) timer_stop(timer_eig);

      for ( int eofID = 0; eofID < n_eig; eofID++ )
	{
	  double *eigenvec = eigenvectors[varID][eofID];

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(varID,nts,eofID,npack,pack,cov,datafields,eigenvec)
#endif 
	  for ( size_t i = 0; i < npack; i++ )
	    {
	      double sum = 0;
	      for ( int j = 0; j < nts; j++ )
		sum += datafields[varID][j][pack[i]] * cov[eofID][j];

	      eigenvec[pack[i]] = sum;
	    }

	  // NORMALIZING
	  double sum = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none)  shared(eigenvec,weight,pack,npack,gridsizemax) reduction(+:sum)
#endif 
	  for ( size_t i = 0; i < npack; i++ )
	    sum +=  weight[pack[i]%gridsizemax] *
	            eigenvec[pack[i]] * eigenvec[pack[i]];

	  if ( sum > 0 )
	    {
	      sum = sqrt(sum);
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(sum,npack,eigenvec,pack)
#endif
	      for ( size_t i = 0; i < npack; i++ )
		eigenvec[pack[i]] /= sum;
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(eigenvec,pack,missval,npack)
#endif
	      for ( size_t i = 0; i < npack; i++ )
		eigenvec[pack[i]] = missval;
	    }
	}     /* for ( eofID = 0; eofID < n_eig; eofID++ )     */

      if ( eigv ) Free(eigv);
      for ( int i = 0; i < n; i++ )
	if ( cov[i] ) Free(cov[i]);
    }         /* for ( varID = 0; varID < nvars; varID++ )    */

  /* write files with eigenvalues (ID3) and eigenvectors (ID2) */

  /*  eigenvalues */
  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  int vlistID2 = vlistDuplicate(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  double xvals = 0, yvals = 0;
  gridDefXvals(gridID2, &xvals);
  gridDefYvals(gridID2, &yvals);

  ngrids = vlistNgrids(vlistID2);
  for ( int i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID2, i, gridID2);

  int zaxisID2 = zaxisCreate(ZAXIS_GENERIC, 1);
  double zvals = 0;
  zaxisDefLevels(zaxisID2, &zvals);
  zaxisDefName(zaxisID2, "zaxis_Reduced");
  zaxisDefLongname(zaxisID2, "Reduced zaxis from EOF3D - only one eigen value per 3D eigen vector");

  int nzaxis = vlistNzaxis(vlistID2);
  for ( int i = 0; i < nzaxis; i++ )
    vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

  /*  eigenvectors */
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID3 = vlistDuplicate(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID3);

  int vdate = 10101;
  int vtime = 0;
  juldate_t juldate = juldate_encode(calendar, vdate, vtime);
  for ( tsID = 0; tsID < n; tsID++ )
    {
      juldate = juldate_add_seconds(60, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      taxisDefVdate(taxisID2, vdate);
      taxisDefVtime(taxisID2, vtime);
      pstreamDefTimestep(streamID2, tsID);

      if ( tsID < n_eig )
        {
          taxisDefVdate(taxisID3, vdate);
          taxisDefVtime(taxisID3, vtime);
          pstreamDefTimestep(streamID3, tsID);
        }

      for ( varID = 0; varID < nvars; varID++ )
        {
          double missval = vlistInqVarMissval(vlistID1, varID);
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < (int)nlevs; levelID++ )
            {
	      size_t offset = levelID * gridsizemax;
              if ( tsID < n_eig )
                {
                  nmiss = 0;
                  for ( size_t i = 0; i < gridsizemax; i++ )
                    if ( DBL_IS_EQUAL(eigenvectors[varID][tsID][offset + i], missval) ) nmiss++;

                  pstreamDefRecord(streamID3, varID, levelID);
                  pstreamWriteRecord(streamID3, &eigenvectors[varID][tsID][offset], nmiss);
                }
	    }

	  nmiss = (DBL_IS_EQUAL(eigenvalues[varID][tsID][0], missval)) ? 1 : 0;

	  pstreamDefRecord(streamID2, varID, 0);
	  pstreamWriteRecord(streamID2, eigenvalues[varID][tsID],nmiss);
        } // for ( varID = 0; ... )
    } // for ( tsID = 0; ... )

  for ( varID = 0; varID < nvars; varID++)
    {
      for ( int i = 0; i < nts; i++)
	{
	  Free(datafields[varID][i]);
	  if ( i < n_eig )
	    Free(eigenvectors[varID][i]);
	  Free(eigenvalues[varID][i]);
	}

      Free(datafields[varID]);
      Free(datacounts[varID]);
      Free(eigenvectors[varID]);
      Free(eigenvalues[varID]);
    }

  Free(datafields);
  Free(datacounts);
  Free(eigenvectors);
  Free(eigenvalues);
  Free(in);

  Free(pack);
  Free(weight);

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
  
  cdoFinish();
 
  return 0;
}
