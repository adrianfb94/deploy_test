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
   Ensstat3       ensrkhist_space  Ensemble ranked histogram averaged over time
   Ensstat3       ensrkhist_time   Ensemble ranked histogram averaged over space
   Ensstat3       ensroccurve      Ensamble Receiver Operating Characteristics
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"

// Defines for rank histogram
enum TDATA_TYPE {TIME, SPACE};
#define time_data TIME
#define space_data SPACE


// Defines for Receiver Operating Characteristics (ROC)
#define DEBUG_ROC 0
enum CONTINGENCY_TYPE {TP, FP, FN, TN};
/* TP - True positive  ( event     forecast and     occured)  HIT               */
/* FP - False positive ( event     forecast and not occured)  FALSE ALARM       */
/* TN - True negative  ( event not forecast and not ocurred)  CORRECT REJECTION */
/* FN - False negative ( event not forecast and     ocurred)  MISSED            */

enum ROC_ENUM_TYPE {TPR, FPR};
/* TPR = True Positive Rate = TP / ( TP + FN ) */
/* FNR = False Negtive Rate = FN / ( FP + TN ) */

double roc_curve_integrate(const double **roc, const int n);

void *Ensstat3(void *argument)
{
  int i,j;
  int nrecs = 0, nrecs0, nmiss;
  int cum;
  int chksum;                  // for check of histogram population 
  int levelID, varID, binID = 0;
  int vlistID;
  int gridID, gridID2;
  int have_miss = 0;
  int streamID = 0, streamID2 = 0;
  int **array2 = NULL;
  int **ctg_tab = NULL, *hist = NULL;         // contingency table and histogram
  double missval;
  double *dat;                  // pointer to ensemble data for ROC
  double *uThresh = NULL, *lThresh = NULL;    // thresholds for histograms
  double **roc = NULL;                 // receiver operating characteristics table
  double val;
  int ival;

  typedef struct
  {
    int streamID;
    int vlistID;
    double *array;
  } ens_file_t;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("ensroc",          func_roc,  0,          NULL);
  cdoOperatorAdd("ensrkhist_space", func_rank, space_data, NULL);
  cdoOperatorAdd("ensrkhist_time",  func_rank, time_data,  NULL);
  // clang-format on
  
  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int datafunc = cdoOperatorF2(operatorID);

  int nbins = 0;
  if ( operfunc == func_roc ) {
    operatorInputArg("Number of eigen functions to write out");
    nbins = parameter2int(operatorArgv()[0]);
  }
  
  int nfiles = cdoStreamCnt() - 1;
  int nens = nfiles-1;

  if ( cdoVerbose )
    cdoPrint("Ensemble over %d files.", nfiles);

  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename) )
    cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t*) Malloc(nfiles*sizeof(ens_file_t));

  /* *************************************************** */
  /* should each thread be allocating memory locally???? */
  /* ("first touch strategy")                            */
  /* --> #pragma omp parallel for ...                    */
  /* *************************************************** */
  field_type *field = (field_type*) Malloc(ompNumThreads*sizeof(field_type));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      field_init(&field[i]);
      field[i].size   = nfiles;
      field[i].ptr    = (double*) Malloc(nfiles*sizeof(double));
      field[i].weight = (double*) Malloc(nfiles*sizeof(double));
      for ( int fileID = 0; fileID < nfiles; fileID++ )
	field[i].weight[fileID] = 1;
    }

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = pstreamOpenRead(cdoStreamName(fileID));

      vlistID = pstreamInqVlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
    }

  /* check for identical contents of all ensemble members */
  for ( int fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  int vlistID1 = ef[0].vlistID;
  int vlistID2 = vlistCreate();
  int nvars = vlistNvars(vlistID1);
  int *varID2 = (int*) Malloc(nvars*sizeof(int));

  double *levs = (double*) Calloc(nfiles, sizeof(double));
  int zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nfiles);
  for ( i=0; i<nfiles; i++ )
    levs[i] = i;
  zaxisDefLevels(zaxisID2,levs);
  zaxisDefName(zaxisID2, "histogram_binID");

  int time_mode = datafunc == TIME? TIME_VARYING : TIME_CONSTANT;

  for ( varID=0; varID<nvars; varID++) {

    /* **************************************************************** */
    /* nfiles includes the observation, so there are nfiles-1 ensembles */
    /* and exactly nfiles bins, in which the observation could fall     */
    /* **************************************************************** */

    if ( datafunc == TIME ) 
      {
	val = 0;
	gridID2 = gridCreate(GRID_LONLAT, 1);
	gridDefXsize(gridID2, 1);
	gridDefYsize(gridID2, 1);
	gridDefXvals(gridID2, &val);
	gridDefYvals(gridID2, &val);
      }
    else // datafunc == SPACE
      gridID2 = vlistInqVarGrid(vlistID1,varID);

    varID2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, time_mode);
  }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  
  for ( varID=0; varID< nvars; varID++ ){
    if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) > 1 ) {
      cdoWarning("More than one level not supported when processing ranked histograms.");
      cdoWarning("Try to use `cdo splitlevel` to split the dataset into levels and apply");
      cdoWarning("the operator seperately to each level.");
      cdoAbort("Exit due to unsupported file structure");
    }
  } 

  if ( operfunc != func_roc )
    {
      streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
      pstreamDefVlist(streamID2, vlistID2);
    }

  int gridsize = vlistGridsizeMax(vlistID1);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double*) Malloc(gridsize*sizeof(double));

  if ( operfunc == func_rank && datafunc == SPACE ) 
    { /*need to memorize data for entire grid before writing          */
      array2 = (int**) Malloc((nfiles+1)*sizeof(int*));
      for ( binID=0; binID<nfiles; binID++ ) 
	array2[binID] = (int*) Calloc( gridsize, sizeof(int));
    }
  else if ( operfunc == func_rank )
    {  /* can process data separately for each timestep and only need */
       /* to cumulate values over the grid                            */
      array2    = (int **) Malloc( (nfiles+1)*sizeof(int *));
      for ( binID=0; binID<nfiles; binID++ )
	array2[binID] = (int*) Calloc( 1, sizeof(int));
    }

  if ( operfunc == func_roc ) {
    ctg_tab = (int**) Malloc((nbins+1)*sizeof(int*));
    hist =    (int*) Malloc(nbins*sizeof(int));
    uThresh = (double*) Malloc(nbins*sizeof(double));
    lThresh = (double*) Malloc(nbins*sizeof(double));
    roc     = (double**) Malloc((nbins+1)*sizeof(double*));
    
    for  ( i=0; i<nbins; i++ ) {
      ctg_tab[i] = (int*) Calloc( 4,sizeof(int));
      roc[i]     = (double*) Calloc( 2,sizeof(double));
      uThresh[i] = ((double)i+1)/nbins;
      lThresh[i] = (double)i/nbins;
    }
    ctg_tab[nbins] = (int*) Calloc(4,sizeof(int));
    roc[nbins]     = (double*) Calloc(2,sizeof(double));
  }
  
  
  int tsID = 0;
  do
    {
      nrecs0 = pstreamInqTimestep(ef[0].streamID, tsID);
      for ( int fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = pstreamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    {
	      if ( nrecs == 0 )
		cdoAbort("Inconsistent ensemble file, too few time steps in %s!", cdoStreamName(fileID)->args);
	      else
		cdoAbort("Inconsistent ensemble file, number of records at time step %d of %s and %s differ!",
			   tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
	    }
	}

      if ( operfunc == func_rank && ( datafunc == TIME || tsID == 0 ) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);
	  if ( nrecs0 > 0 ) pstreamDefTimestep(streamID2, tsID);
	}

      //      fprintf(stderr,"TIMESTEP %i varID %i rec %i\n",tsID,varID,recID);
      
      for ( int recID = 0; recID < nrecs0; recID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(ef, nfiles)      \
                                      private(streamID, nmiss) \
                                  lastprivate(varID, levelID)
#endif
	  for ( int fileID = 0; fileID < nfiles; fileID++ )
	    {
	      streamID = ef[fileID].streamID;
	      pstreamInqRecord(streamID, &varID, &levelID);
	      pstreamReadRecord(streamID, ef[fileID].array, &nmiss);
	    }

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  nmiss = 0;
	  if ( datafunc == TIME && operfunc == func_rank) 
	    for ( binID=0;binID<nfiles;binID++ )
	      array2[binID][0] = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(binID)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
	      int ompthID = cdo_omp_get_thread_num();

	      field[ompthID].missval = missval;
	      field[ompthID].nmiss = 0;
	      have_miss = 0;
	      for ( int fileID = 0; fileID < nfiles; fileID++ )
		{
		  field[ompthID].ptr[fileID] = ef[fileID].array[i];
		  if ( DBL_IS_EQUAL(field[ompthID].ptr[fileID], missval) ) 
		    {
		      have_miss = 1;
		      break;
		    }
		}
	      
	      // need to ignore all data for a gridpoint if a single ensemble
	      // has a missing value at that gridpoint.
	      if ( ! have_miss )  // only process if no missing value in ensemble
		{
		  switch( operfunc ) 
		    {
		    case ( func_rank ): 
		      /* ****************/
		      /* RANK HISTOGRAM */
		      /* ************** */
		      // for ( j=0; j<nfiles; j++ )
		      //   fprintf(stderr,"%5.2g ",field[ompthID].ptr[j]);
#if defined(_OPENMP)
#pragma omp critical
#endif
		      binID = (int) fldfun(field[ompthID], operfunc);
		      // fprintf(stderr,"-->%i\n",binID);
		      
		      if ( datafunc == SPACE && ! have_miss) 
			array2[binID][i]++;
		      else if ( ! have_miss ) 
			array2[binID][0]++;
		      break;

		    case ( func_roc ):
		      /* ********************************** */
		      /* RECEIVER OPERATING CHARACTERISTICS */
		      /* ********************************** */
		      dat = &field[ompthID].ptr[1];
		      ival = field[ompthID].ptr[0] > 0.5? 1 : 0;

		      for ( binID=0; binID<nbins; binID++ )
			hist[binID] = 0;

		      for ( j=0; j<nens; j++ ) 
			for ( binID=0; binID<nbins; binID++ ) 
			  if ( dat[j] >= lThresh[binID] && dat[j] < uThresh[binID] )
			    hist[binID]++;

		      chksum = 0;
		      for ( binID=0; binID<nbins; binID++ )
			chksum += hist[binID];

		      if ( chksum != nens )  exit(1);

		      cum = 0;
		      if ( ival == 1 ) 
			{
			  // all true positives in first bin
			  ctg_tab[0][TP] += nens;
			  
			  cum += hist[0];
			  for ( binID=1; binID<nbins; binID++ ) 
			    {
			      ctg_tab[binID][TP] += nens-cum;
			      ctg_tab[binID][FN] += cum;
			      cum += hist[binID];
			    }
			  ctg_tab[binID][TP] += nens-cum;
			  ctg_tab[binID][FN] += cum;
			}
		      else if ( ival == 0 ) 
			{
			  // all false positives in first bin
			  ctg_tab[0][FP] += nens;
			  cum += hist[0];
			  for ( binID=1; binID<nbins; binID++ ) 
			    {
			      ctg_tab[binID][FP] += nens-cum;
			      ctg_tab[binID][TN] += cum;
			      cum += hist[binID];
			    }
			  ctg_tab[binID][FP] += nens-cum;
			  ctg_tab[binID][TN] += cum;
			}
		      break;
		      
		    }// switch ( operfunc )
		}    // if ( ! have_miss ) 
	    }        // for ( i=0; i<gridsize; i++ )

	  if ( datafunc == TIME && operfunc == func_rank ) 
	    {
	      for ( binID=0; binID<nfiles; binID++ ) {
		val = (double)array2[binID][0];
		//		fprintf(stderr,"%i ",(int)val);
		pstreamDefRecord(streamID2, varID2[varID], binID);
		pstreamWriteRecord(streamID2,&val, nmiss);
	      }
	      //fprintf(stderr,"\n");
	    }
	  else if ( operfunc == func_roc ) 
	    {
	      if ( DEBUG_ROC ) {
		fprintf(stderr, "#             :     TP     FP     FN     TN         TPR        FPR\n");
		
		for ( binID=0; binID<= nbins; binID++ )  {
		  int p = ctg_tab[binID][TP] + ctg_tab[binID][FN];
		  int n = ctg_tab[binID][FP] + ctg_tab[binID][TN];
		  double tpr = ctg_tab[binID][TP]/(double) p;
		  double fpr = ctg_tab[binID][FP]/(double) n;
		  chksum += ctg_tab[binID][0] + ctg_tab[binID][1] + ctg_tab[binID][2] + ctg_tab[binID][3];
		  
		  roc[binID][TPR] = tpr;
		  roc[binID][FPR] = fpr;
		  
		  fprintf(stderr, "%3i %10.4g: %6i %6i %6i %6i: %10.4g %10.4g\n",
			  binID,binID<nbins?lThresh[binID]:1,
			  ctg_tab[binID][0],ctg_tab[binID][1],ctg_tab[binID][2],ctg_tab[binID][3],
			  tpr,fpr);
		}
		fprintf(stderr,"nbins %10i\n",nbins);
		fprintf(stderr,"#ROC CurveArea: %10.6f\n",
			roc_curve_integrate((const double **)roc,nbins));
	      } // if ( DEBUG_ROC )
	    }   // else if (operfunc == func_roc )
	}       // for ( recID=0; recID<nrecs; recID++ )
      tsID++;
    }           // do [...]
  while ( nrecs0 > 0 );



  if ( operfunc == func_rank )
    {
      double *tmpdoub;
      int osize = gridsize;

      if ( datafunc == TIME ) osize = 1;
      tmpdoub = (double*) Malloc(osize*sizeof(double));

      for ( binID = 0; binID < nfiles; binID++ )
	{
	  for ( i = 0; i < osize; i++ )
	    tmpdoub[i] = (double) array2[binID][i];

	  pstreamDefRecord(streamID2, varID2[varID], binID);
	  pstreamWriteRecord(streamID2, tmpdoub, nmiss);
	}

      Free(tmpdoub);
    }
  else if ( operfunc == func_roc ) {
    fprintf(stdout, "#             :     TP     FP     FN     TN         TPR        FPR\n");
    
    for ( i=0; i<= nbins; i++ )  {
      int p = ctg_tab[i][TP] + ctg_tab[i][FN];
      int n = ctg_tab[i][FP] + ctg_tab[i][TN];
      double tpr = ctg_tab[i][TP]/(double) p;
      double fpr = ctg_tab[i][FP]/(double) n;
      chksum += ctg_tab[i][0] + ctg_tab[i][1] + ctg_tab[i][2] + ctg_tab[i][3];
      
      roc[i][TPR] = tpr;
      roc[i][FPR] = fpr;
      
      int sum = ctg_tab[i][TP] + ctg_tab[i][TN] + ctg_tab[i][FP] + ctg_tab[i][FN];

      fprintf(stdout, "%3i %10.4g: %6i %6i %6i %6i (%6i): %10.4g %10.4g\n",
	      i,i<nbins?lThresh[i]:1,
	      ctg_tab[i][0],ctg_tab[i][1],ctg_tab[i][2],ctg_tab[i][3],sum,
	      tpr,fpr);
    }
 
    fprintf(stdout,"#ROC CurveArea: %10.6f\n",
	    roc_curve_integrate((const double **)roc,nbins));
  }

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = ef[fileID].streamID;
      pstreamClose(streamID);
    }

  if ( operfunc != func_roc ) 
    pstreamClose(streamID2);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) Free(ef[fileID].array);

  if ( ef ) Free(ef);
  if ( array2 ) {
    for (binID=0; binID<nfiles; binID++ )
      Free(array2[binID]);
    Free(array2);
  }

  for ( i = 0; i < ompNumThreads; i++ )
    {
      if ( field[i].ptr    ) Free(field[i].ptr);
      if ( field[i].weight ) Free(field[i].weight);
    }
  
  if ( field ) Free(field);
  
  cdoFinish();
  
  return 0;
}


double roc_curve_integrate(const double **roc, const int n) {
  double y1, y0, x1,x0, dx, dy;
  double step_area;
  double area = 0;

  for ( int i = 1; i <= n; ++i ) {
    x1 = roc[i][FPR]; x0 = roc[i-1][FPR];
    y1 = roc[i][TPR]; y0 = roc[i-1][TPR];
    dx = x1-x0;
    dy = y1-y0;

    step_area = -0.5*dx*dy - dx*y0;
    area += step_area;
  }

  return area-0.5;

}
