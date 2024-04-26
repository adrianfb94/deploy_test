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

     Timeof        eof             EOF in spatial or time space
     Timeof        eofspatial      EOF in spatial space
     Timeof        eoftime         EOF in time space
*/
/*
 * TODO: 
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <limits.h>  // LONG_MAX
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "statistic.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

static
void scale_eigvec_grid(double *restrict out, int tsID, size_t npack, const size_t *restrict pack, const double *restrict weight, double **covar, double sum_w)
{
  for ( size_t i = 0; i < npack; ++i )
    out[pack[i]] = covar[tsID][i] / sqrt(weight[pack[i]]/sum_w);
}

static
void scale_eigvec_time(double *restrict out, int tsID, int nts, size_t npack, const size_t *restrict pack, const double *restrict weight,
		       double **covar, double **data, double missval, double sum_w)
{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(npack, nts, tsID, pack, data, covar, out)
#endif
  for ( size_t i = 0; i < npack; ++i )
    {
      double sum = 0;
      for ( int j = 0; j < nts; ++j )
	sum += data[j][i] * covar[tsID][j];
      
      out[pack[i]] = sum;
    }
  /*
  for ( size_t j = 0; j < nts; ++j )
    {
      for ( size_t i = 0; i < npack; ++i )
	out[pack[i]] += data[j][i] * covar[tsID][j];
    }
  */

  // Normalizing
  double sum = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) reduction(+:sum)	\
  shared(out,weight,pack,npack)
#endif
  for ( size_t i = 0; i < npack; ++i )
    {
      // do not need to account for weights as eigenvectors are non-weighted                               
      sum += weight[pack[i]] * out[pack[i]] * out[pack[i]];
    }

  if ( sum > 0 )
    {
      sum = sqrt(sum/sum_w);
#if defined(_OPENMP)
#pragma omp parallel for default(none)  shared(npack,pack,sum,out)
#endif
      for ( size_t i = 0; i < npack; ++i ) out[pack[i]] /= sum;
    }
  else
    {
#if defined(_OPENMP)
#pragma omp parallel for default(none)  shared(npack,pack,out,missval)
#endif
      for ( size_t i = 0; i < npack; ++i ) out[pack[i]] = missval;
    }
}


enum T_EIGEN_MODE get_eigenmode(void)
{
  enum T_EIGEN_MODE eigen_mode = JACOBI;

  char *envstr = getenv("CDO_SVD_MODE");  
  if ( envstr )
    {
      if ( !strncmp(envstr, "danielson_lanczos", 17) ) 
	eigen_mode = DANIELSON_LANCZOS;
      else if ( !strncmp(envstr, "jacobi", 6) )
	eigen_mode = JACOBI;
      else
	{
	  cdoWarning("Unknown environmental setting %s for CDO_SVD_MODE. Available options are", envstr);
	  cdoWarning("  - 'jacobi' for a one-sided parallelized jacobi algorithm");
	  cdoWarning("  - 'danielson_lanzcos' for the D/L algorithm");
	}
    }

  if ( cdoVerbose ) 
    cdoPrint("Using CDO_SVD_MODE '%s' from %s",
	     eigen_mode==JACOBI?"jacobi":"danielson_lanczos",
	     envstr?"Environment":" default");  

#if defined(_OPENMP)
  if ( omp_get_max_threads() > 1 && eigen_mode == DANIELSON_LANCZOS )  {
    cdoWarning("Requested parallel computation with %i Threads ",omp_get_max_threads());
    cdoWarning("  but environmental setting CDO_SVD_MODE causes sequential ");
    cdoWarning("  Singular value decomposition");
  }
#endif 

  return eigen_mode;
}


enum T_WEIGHT_MODE get_weightmode(void)
{  
  enum T_WEIGHT_MODE weight_mode = WEIGHT_OFF;

  char *envstr = getenv("CDO_WEIGHT_MODE");
  if ( envstr )
    {
      if ( !strncmp(envstr, "off", 3) ) 
	weight_mode = WEIGHT_OFF;
      else if ( !strncmp(envstr, "on", 2) )
	weight_mode = WEIGHT_ON;
      else
	cdoWarning("Unknown environmental setting %s for CDO_WEIGHT_MODE. Available options are: on/off", envstr);
    }

  if ( cdoVerbose ) 
    cdoPrint("Using CDO_WEIGHT_MODE '%s' from %s",
	     weight_mode==WEIGHT_OFF?"off":"on",
	     envstr?"Environment":" default");  

  return weight_mode;
}


void *EOFs(void * argument)
{
  enum {EOF_, EOF_TIME, EOF_SPATIAL};

  int nlevs = 0 ;
  int nmiss;
  int varID, levelID;
  int nts = 0;
  size_t n = 0;
  int grid_space = 0, time_space = 0;
  int timer_cov = 0, timer_eig = 0;

  int calendar = CALENDAR_STANDARD;

  typedef struct {
    bool init;
    bool first_call;
    double *eig_val;
    double *covar_array;
    double **covar;
    double **data;
  }
  eofdata_t;

  if ( cdoTimer )
    {
      timer_cov  = timer_new("Timeof cov");
      timer_eig  = timer_new("Timeof eig");
    }
  
  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("eof",        EOF_,        0, NULL);
  cdoOperatorAdd("eoftime",    EOF_TIME,    0, NULL);
  cdoOperatorAdd("eofspatial", EOF_SPATIAL, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  operatorInputArg("Number of eigen functions to write out");
  int n_eig = parameter2int(operatorArgv()[0]);
  
  enum T_EIGEN_MODE eigen_mode = get_eigenmode();
  enum T_WEIGHT_MODE weight_mode = get_weightmode();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int vlistID1  = pstreamInqVlist(streamID1);
  int taxisID1  = vlistInqTaxis(vlistID1);
  int gridID1   = vlistInqVarGrid(vlistID1, 0);
  size_t gridsize  = vlistGridsizeMax(vlistID1);
  int nvars     = vlistNvars(vlistID1);
  int nrecs;

  int ngrids = vlistNgrids(vlistID1);
  for ( int index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      cdoAbort("Too many different grids!");

  /* eigenvalues */

  /* COUNT NUMBER OF TIMESTEPS if EOF_ or EOF_TIME */
  if ( operfunc == EOF_ || operfunc == EOF_TIME )
    {
      if ( cdoVerbose ) cdoPrint("Counting timesteps in ifile");

      nts = vlistNtsteps(vlistID1);
      if ( nts == -1 )
	{
          nts = 0;
	  while ( pstreamInqTimestep(streamID1, nts) ) nts++;

	  if ( cdoVerbose ) cdoPrint("Counted %i timeSteps", nts);

          pstreamClose(streamID1);

          streamID1 = pstreamOpenRead(cdoStreamName(0));
          vlistID1  = pstreamInqVlist(streamID1);
          taxisID1  = vlistInqTaxis(vlistID1);
	}
      else
        if ( cdoVerbose ) cdoPrint("Found %i timeSteps", nts);

      if ( (size_t)nts < gridsize || operfunc == EOF_TIME )
	{
	  time_space = 1;
	  grid_space = 0;
	}
      else
        {
          time_space = 0;
          grid_space = 1;
        }
    }
  else if ( operfunc == EOF_SPATIAL )
    {
      time_space = 0;
      grid_space = 1;
    }

  /* reset the requested number of eigen-function to the maximum if neccessary */
  if ( time_space )
    {
      if ( n_eig > nts )
        {
          cdoWarning("Solving in time-space:");
          cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
          cdoWarning("Setting n_eig to %d.", nts);
          cdoWarning("If You want to force a solution in grid-space use operator eofspatial");
          n_eig = nts;
        }
      n = nts;
    }
  else if ( grid_space )
    {
      if ( ((double)gridsize)*gridsize > (double)LONG_MAX ) cdoAbort("Grid space too large!");

      if ( (size_t)n_eig > gridsize )
        {
          cdoWarning("Solving in spatial space");
          cdoWarning("Number of eigen-functions to write out is bigger than grid size");
          cdoWarning("Setting n_eig to %zu", gridsize);
          cdoWarning("If You want to force a solution in time-space use operator eoftime");
          n_eig = gridsize;
        }
      n = gridsize;
    }

  if ( cdoVerbose ) 
    cdoPrint("Calculating %d eigenvectors and %d eigenvalues in %s",
	     n_eig, n, grid_space==1?"grid_space" : "time_space");

  double *weight = (double *) Malloc(gridsize*sizeof(double));
  for ( size_t i = 0; i < gridsize; ++i ) weight[i] = 1.;

  if ( weight_mode == WEIGHT_ON )
    {
      int wstatus = gridWeights(gridID1, weight);
      if ( wstatus != 0  )
	{
	  weight_mode = WEIGHT_OFF;
	  cdoWarning("Using constant grid cell area weights!");
	}
    }

  /* allocation of temporary fields and output structures */
  size_t npack = ULONG_MAX;
  size_t *pack         = (size_t *) Malloc(gridsize*sizeof(size_t));
  double *in           = (double *) Malloc(gridsize*sizeof(double));
  eofdata_t **eofdata  = (eofdata_t **) Malloc(nvars*sizeof(eofdata_t*));

  for ( varID = 0; varID < nvars; ++varID )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      eofdata[varID] = (eofdata_t *) Malloc(nlevs*sizeof(eofdata_t));

      for ( levelID = 0; levelID < nlevs; ++levelID )
        {
	  eofdata[varID][levelID].init = false;
	  eofdata[varID][levelID].first_call = true;
	  eofdata[varID][levelID].eig_val = NULL;
	  eofdata[varID][levelID].covar_array = NULL;
	  eofdata[varID][levelID].covar = NULL;
	  eofdata[varID][levelID].data = NULL;

	  if ( time_space )
	    eofdata[varID][levelID].data = (double **) Malloc(nts*sizeof(double *));
        }
    }

  if ( cdoVerbose )
    cdoPrint("Allocated eigenvalue/eigenvector structures with nts=%d gridsize=%zu", nts, gridsize);

  double *covar_array = NULL;
  double **covar = NULL;
  double sum_w = 1.;

  int tsID = 0;

  /* read the data and create covariance matrices for each var & level */
  while ( TRUE )
    {
      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, in, &nmiss);

	  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          double missval = vlistInqVarMissval(vlistID1, varID);
	  if ( npack == ULONG_MAX )
	    {
	      npack = 0;
	      for ( size_t i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(weight[i], 0.0) && !DBL_IS_EQUAL(weight[i], missval) &&
		       !DBL_IS_EQUAL(in[i], missval) )
                    {
                      pack[npack] = i;
                      npack++;
                    }
                }

	      if ( weight_mode == WEIGHT_ON )
		{
		  sum_w = 0;
		  for ( size_t i = 0; i < npack; i++ ) sum_w += weight[pack[i]];
		}
	    }

	  size_t ipack = 0;
	  for ( size_t i = 0; i < gridsize; ++i )
	    {
	      if ( !DBL_IS_EQUAL(weight[i], 0.0) && !DBL_IS_EQUAL(weight[i], missval) &&
		   !DBL_IS_EQUAL(in[i], missval) )
		{
                  if ( pack[ipack] != i ) cdoAbort("Missing values unsupported!");
                  ipack++;
		}
	    }
          if ( ipack != npack ) cdoAbort("Missing values unsupported!");

	  if ( grid_space )
            {
	      if ( !eofdata[varID][levelID].init )
		{
                  n = npack;
		  double *covar_array = (double *) Malloc(npack*npack*sizeof(double));
		  covar = (double **) Malloc(npack*sizeof(double *));
		  for ( size_t i = 0; i < npack; ++i ) covar[i] = covar_array + npack*i;
		  for ( size_t i = 0; i < npack; ++i )
                    for ( size_t j = 0; j < npack; ++j ) covar[i][j] = 0;

		  eofdata[varID][levelID].covar_array = covar_array;
		  eofdata[varID][levelID].covar       = covar;
		}
	      else
		{
		  covar = eofdata[varID][levelID].covar;
		}
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( size_t ipack = 0; ipack < npack; ++ipack )
                for ( size_t jpack = ipack; jpack < npack; ++jpack )
                  covar[ipack][jpack] += in[pack[ipack]] * in[pack[jpack]];
	    }
          else if ( time_space )
	    {
	      double *data = (double *) Malloc(npack*sizeof(double));
	      eofdata[varID][levelID].data[tsID] = data;

	      for ( size_t ipack = 0; ipack < npack; ipack++ )
		data[ipack] = in[pack[ipack]];
	    }

	  eofdata[varID][levelID].init = true;
        }

      tsID++;
    }

  if ( grid_space ) nts = tsID;

  if ( tsID == 1 ) cdoAbort("File consists of only one timestep!");

  /* write files with eigenvalues (ID3) and eigenvectors (ID2) */

  /* eigenvalues */
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
  for ( int i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID2, i, gridID2);

  /*  eigenvectors */
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID3  = vlistDuplicate(vlistID1);
  int taxisID3  = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID3);

  int vdate = 10101;
  int vtime = 0;
  juldate_t juldate = juldate_encode(calendar, vdate, vtime);

  double *out = in;
  double *eig_val = NULL;

  int nts_out = nts;
  if ( npack < (size_t)nts ) nts_out = npack;

  for ( tsID = 0; tsID < nts_out; tsID++ )
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
	  char vname[256];
	  vlistInqVarName(vlistID1, varID, vname);
	  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          double missval = vlistInqVarMissval(vlistID1, varID);

	  for ( levelID = 0; levelID < nlevs; levelID++ )
	    {
	      double **data = eofdata[varID][levelID].data;

	      if ( eofdata[varID][levelID].first_call )
		{
		  eofdata[varID][levelID].first_call = false;

		  if ( cdoVerbose )
		    cdoPrint("Calculating covar matrices for %i levels of var%i (%s)", nlevs, varID, vname);
		  
		  if ( cdoTimer ) timer_start(timer_cov);

		  if ( cdoVerbose ) cdoPrint("processing level %i",levelID);
		  
		  if ( grid_space )
		    {
		      if ( npack )
			{
			  eig_val = (double *) Malloc(npack*sizeof(double));
			  eofdata[varID][levelID].eig_val = eig_val;
			}

		      covar = eofdata[varID][levelID].covar;

		      for ( size_t ipack = 0; ipack < npack; ++ipack )
			{
                          size_t j;
			  size_t i = pack[ipack];
			  for ( size_t jpack = 0; jpack < npack; ++jpack )
			    {
			      if ( jpack < ipack )
				{
				  covar[ipack][jpack] = covar[jpack][ipack];
				}
			      else
				{
				  j = pack[jpack];
				  covar[ipack][jpack] = 
				    covar[ipack][jpack] *   // covariance
				    sqrt(weight[i]) * sqrt(weight[j]) / sum_w /       // weights
				    nts;   // number of data contributing
				}
			    }
			}
		    }
		  else if ( time_space )
		    {		      
		      if ( cdoVerbose )
			cdoPrint("allocating covar with %i x %i elements | npack=%zu", nts, nts, npack);

		      covar_array = (double *) Malloc(nts*nts*sizeof(double));
		      covar = (double **) Malloc(nts*sizeof(double *));
		      for ( int i = 0; i < nts; ++i ) covar[i] = covar_array + nts*i;

		      eig_val = (double *) Malloc(nts*sizeof(double));
		      eofdata[varID][levelID].eig_val     = eig_val;
		      eofdata[varID][levelID].covar_array = covar_array;
		      eofdata[varID][levelID].covar       = covar;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
		      for ( int j1 = 0; j1 < nts; j1++ )
			{
                          double *df1p = data[j1];
			  for ( int j2 = 0; j2 < j1; j2++ ) covar[j1][j2] = covar[j2][j1];
			  for ( int j2 = j1; j2 < nts; j2++ )
			    {
			      double *df2p = data[j2];
			      double sum = 0;
			      for ( size_t i = 0; i < npack; i++ )
				sum += weight[pack[i]]*df1p[i]*df2p[i];
			      covar[j1][j2] = sum / sum_w / nts;
			    }
			}
		      
		      if ( cdoVerbose )
			cdoPrint("finished calculation of covar-matrix for var %s", vname);
		    }

		  if ( cdoTimer ) timer_stop(timer_cov);

		  /* SOLVE THE EIGEN PROBLEM */
		  if ( cdoTimer ) timer_start(timer_eig);

		  if ( eigen_mode == JACOBI ) 
		    // TODO: use return status (>0 okay, -1 did not converge at all) 
		    parallel_eigen_solution_of_symmetric_matrix(covar, eig_val, n, __func__);
		  else 
		    eigen_solution_of_symmetric_matrix(covar, eig_val, n, __func__);

		  if ( cdoTimer ) timer_stop(timer_eig);
		  /* NOW: covar contains the eigenvectors, eig_val the eigenvalues */

		  for ( size_t i = 0; i < gridsize; ++i ) out[i] = missval;
	  
		  // for ( int i = 0; i < n; i++ ) eig_val[i] *= sum_w;
		} // first_call
	      else
		{
		  covar   = eofdata[varID][levelID].covar;
		  eig_val = eofdata[varID][levelID].eig_val;
		}

              if ( tsID < n_eig )
                {
		  if      ( grid_space ) scale_eigvec_grid(out, tsID, npack, pack, weight, covar, sum_w);
		  else if ( time_space ) scale_eigvec_time(out, tsID, nts, npack, pack, weight, covar, data, missval, sum_w);

                  nmiss = 0;
                  for ( size_t i = 0; i < gridsize; i++ ) if ( DBL_IS_EQUAL(out[i], missval) ) nmiss++;

                  pstreamDefRecord(streamID3, varID, levelID);
                  pstreamWriteRecord(streamID3, out, nmiss);
		} // loop n_eig

              nmiss = (DBL_IS_EQUAL(eig_val[tsID], missval)) ? 1 : 0;

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, &eig_val[tsID], nmiss);
	    } // loop nlevs
	} // loop nvars
    }

  for ( varID = 0; varID < nvars; varID++)
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      
      for( levelID = 0; levelID < nlevs; levelID++ )
        { 	  
	  if ( eofdata[varID][levelID].eig_val ) Free(eofdata[varID][levelID].eig_val);
	  if ( eofdata[varID][levelID].covar_array ) Free(eofdata[varID][levelID].covar_array);
	  if ( eofdata[varID][levelID].covar ) Free(eofdata[varID][levelID].covar);
	  if ( time_space && eofdata[varID][levelID].data )
	    {
	      for ( tsID = 0; tsID < nts; tsID++ )
		if ( eofdata[varID][levelID].data[tsID] ) Free(eofdata[varID][levelID].data[tsID]);
	      Free(eofdata[varID][levelID].data);
	    }
	}
      
      Free(eofdata[varID]);
    }

  Free(eofdata);
  Free(in);
  Free(pack);
  Free(weight);

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  gridDestroy(gridID2);

  //  taxisDestroy(taxisID2);
  //  taxisDestroy(taxisID3);

  cdoFinish();

  return 0;
}
