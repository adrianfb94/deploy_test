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

*/
#include <time.h> // clock()

#include <cdi.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "grid_search.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "clipping/geometry.h"
#ifdef __cplusplus
}
#endif


void fillmiss(field_type *field1, field_type *field2, int nfill)
{
  int nx, ny, i, j;
  int nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  int globgrid = FALSE,gridtype;
  double s1, s2;
  double xr, xu, xl, xo;
  double **matrix1, **matrix2;

  int gridID = field1->grid;
  int nmiss1 = field1->nmiss;
  double missval  = field1->missval;
  double *array1  = field1->ptr;
  double *array2  = field2->ptr;

  nx       = gridInqXsize(gridID);
  ny       = gridInqYsize(gridID);
  globgrid = gridIsCircular(gridID);

  gridtype = gridInqType(gridID);
  if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) )
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  matrix1 = (double **) Malloc(ny * sizeof(double *));
  matrix2 = (double **) Malloc(ny * sizeof(double *));

  for ( j = 0; j < ny; j++ )
    {
      matrix1[j] = array1 + j*nx;
      matrix2[j] = array2 + j*nx;
    }

  for ( j = 0; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      {
	if ( DBL_IS_EQUAL(matrix1[j][i], missval) )
	  {
	    nmiss2++;

	    kr = ku = kl = ko = 0;
	    xr = xu = xl = xo = 0.;

	    for ( ir = i + 1; ir < nx; ir++ )
	      if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
		{ kr = ir - i; xr = matrix1[j][ir]; break; }

	    if ( globgrid && ir == nx )
	      {
		for ( ir = 0; ir < i; ir++ )
		  if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
		    { kr = nx + ir - i; xr = matrix1[j][ir]; break; }
	      }

	    for ( il = i-1; il >= 0; il-- )
	      if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
		{ kl = i - il; xl = matrix1[j][il]; break; }

	    if ( globgrid && il == -1 )
	      {
		for ( il = nx-1; il > i; il-- )
		  if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
		    { kl = nx + i - il; xl = matrix1[j][il]; break; }
	      }

	    for ( iu = j + 1; iu < ny; iu++ )
	      if ( !DBL_IS_EQUAL(matrix1[iu][i], missval) )
		{ ku = iu - j; xu = matrix1[iu][i]; break; }
	    
	    for ( io = j - 1; io >= 0; io-- )
	      if ( !DBL_IS_EQUAL(matrix1[io][i], missval) )
		{ ko = j - io; xo = matrix1[io][i]; break; }
	    
	    /*  printf("%d %d %d %d %d %d %g %g %g %g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo);*/

	    kh = kl + kr;
	    kv = ko + ku;
	    if      ( kh == 0 ) { s1 = 0.; k1 = 0; }
	    else if ( kl == 0 ) { s1 = xr; k1 = 1; }
	    else if ( kr == 0 ) { s1 = xl; k1 = 1; }
	    else { s1 = xr*kl/kh + xl*kr/kh; k1 = 2; }

	    if      ( kv == 0 ) { s2 = 0.; k2 = 0; }
	    else if ( ku == 0 ) { s2 = xo; k2 = 1; }
	    else if ( ko == 0 ) { s2 = xu; k2 = 1; }
	    else { s2 = xu*ko/kv + xo*ku/kv; k2 = 2; }

	    kk = k1 + k2;
	    if ( kk >= nfill )
	      {
		if      ( kk == 0 ) cdoAbort("no point found!");
		else if ( k1 == 0 ) matrix2[j][i] = s2;
		else if ( k2 == 0 ) matrix2[j][i] = s1;
		else  matrix2[j][i] = s1*k2/kk + s2*k1/kk;
	      }
	    else
	      matrix2[j][i] = matrix1[j][i];

	    /* matrix1[j][i] = matrix2[j][i]; */
	  }
	else
	  {
	    matrix2[j][i] = matrix1[j][i];
	  }
      }

  if ( nmiss1 != nmiss2 ) cdoAbort("found only %d of %d missing values!", nmiss2, nmiss1);

  Free(matrix2);
  Free(matrix1);
}


void fillmiss_one_step(field_type *field1, field_type *field2, int maxfill)
{
  int gridID, nx, ny, i, j;
  int nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  double s1, s2;
  double xr, xu, xl, xo;
  double missval;
  double *array1, *array2;
  double **matrix1, **matrix2;

  gridID  = field1->grid;
  missval = field1->missval;
  array1  = field1->ptr;
  array2  = field2->ptr;

  nx  = gridInqXsize(gridID);
  ny  = gridInqYsize(gridID);

  matrix1 = (double **) Malloc(ny * sizeof(double *));
  matrix2 = (double **) Malloc(ny * sizeof(double *));

  for ( j = 0; j < ny; j++ ) { matrix1[j] = array1 + j*nx; matrix2[j] = array2 + j*nx; }

  for (int fill_iterations=0; fill_iterations < maxfill; fill_iterations++) {
  for ( j = 0; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      {
        if ( DBL_IS_EQUAL(matrix1[j][i], missval) )
          {
            nmiss2++;

            kr = ku = kl = ko = 0;
            xr = xu = xl = xo = 0.;

            for ( ir = i + 1; ir < nx; ir++ )
              if ( !DBL_IS_EQUAL(matrix1[j][ir], missval) )
                { kr = ir - i; xr = matrix1[j][ir]; break; }

            for ( il = i-1; il >= 0; il-- )
              if ( !DBL_IS_EQUAL(matrix1[j][il], missval) )
                { kl = i - il; xl = matrix1[j][il]; break; }


            for ( iu = j + 1; iu < ny; iu++ )
              if ( !DBL_IS_EQUAL(matrix1[iu][i], missval) )
                { ku = iu - j; xu = matrix1[iu][i]; break; }

            for ( io = j - 1; io >= 0; io-- )
              if ( !DBL_IS_EQUAL(matrix1[io][i], missval) )
                { ko = j - io; xo = matrix1[io][i]; break; }


            kh = kl + kr;
            kv = ko + ku;
            if      ( kh == 0 ) { s1 = 0.; k1 = 0; }
            else if ( kl == 0 ) { s1 = xr; k1 = kr; }
            else if ( kr == 0 ) { s1 = xl; k1 = kl; }
            else
              {
                if ( kl < kr )
                {
                  s1 = xl;
                  k1 = kl;
                }
              else
                {
                  s1 = xr;
                  k1 = kr;
                }
              }

            if      ( kv == 0 ) { s2 = 0.; k2 = 0; }
            else if ( ku == 0 ) { s2 = xo; k2 = ko; }
            else if ( ko == 0 ) { s2 = xu; k2 = ku; }
            else
              {
                if ( ku < ko )
                  {
                    s2 = xu;
                    k2 = ku;
                  }
                else
                  {
                    s2 = xo;
                    k2 = ko;
                  }
              }

            kk = k1 + k2;
            if      ( kk == 0 ) matrix2[j][i] = matrix1[j][i];
            else if ( k1 == 0 ) matrix2[j][i] = s2;
            else if ( k2 == 0 ) matrix2[j][i] = s1;
            else
              {
                if ( k1 <= k2 )
                {
                  matrix2[j][i] = s1;
                }
                else
                {
                  matrix2[j][i] = s2;
                }

              }

            //printf("%d %d %2d %2d %2d %2d %2g %2g %2g %2g %2g %2g %2g\n", j,i,kr,kl,ku,ko,xr,xl,xu,xo,s1,s2,matrix2[j][i]);
            /* matrix1[j][i] = matrix2[j][i]; */
          }
        else
          {
            matrix2[j][i] = matrix1[j][i];
          }
      }
  for ( j = 0; j < ny; j++ ) for ( i = 0; i < nx; i++ ) matrix1[j][i] = matrix2[j][i];
  }

  Free(matrix2);
  Free(matrix1);
}


int grid_search_nbr(struct gridsearch *gs, size_t num_neighbors, size_t *restrict nbr_add, double *restrict nbr_dist, double plon, double plat);
double nbr_compute_weights(size_t num_neighbors, const int *restrict src_grid_mask, bool *restrict nbr_mask, const size_t *restrict nbr_add, double *restrict nbr_dist);
size_t nbr_normalize_weights(size_t num_neighbors, double dist_tot, const bool *restrict nbr_mask, size_t *restrict nbr_add, double *restrict nbr_dist);

static
void setmisstodis(field_type *field1, field_type *field2, int num_neighbors)
{
  int gridID = field1->grid;
  int gridID0 = gridID;
  double missval = field1->missval;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;

  unsigned gridsize = gridInqSize(gridID);

  unsigned nmiss = field1->nmiss;
  unsigned nvals = gridsize - nmiss;

  double *xvals = (double*) Malloc(gridsize*sizeof(double));
  double *yvals = (double*) Malloc(gridsize*sizeof(double));

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 0);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_radian(units, gridsize, xvals, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_radian(units, gridsize, yvals, "grid center lat");

  size_t *mindex = nmiss ? (size_t *) Calloc(1, nmiss*sizeof(size_t)) : NULL;
  size_t *vindex = nvals ? (size_t *) Calloc(1, nvals*sizeof(size_t)) : NULL;
  double *lons = nvals ? (double *) Malloc(nvals*sizeof(double)) : NULL;
  double *lats = nvals ? (double *) Malloc(nvals*sizeof(double)) : NULL;
  
  unsigned nv = 0, nm = 0;
  for ( unsigned i = 0; i < gridsize; ++i ) 
    {
      array2[i] = array1[i];
      if ( DBL_IS_EQUAL(array1[i], missval) )
        {
          mindex[nm] = i;
          nm++;
        }
      else
        {
          lons[nv] = xvals[i];
          lats[nv] = yvals[i];
          vindex[nv] = i;
          nv++;
        }
    }

  if ( nv != nvals ) cdoAbort("Internal problem, number of valid values differ!");
  

  NEW_2D(bool, nbr_mask, ompNumThreads, num_neighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, num_neighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, num_neighbors); // angular distance four nearest neighbors

  clock_t start, finish;
  start = clock();

  struct gridsearch *gs = NULL;

  if ( nmiss )
    {
      if ( num_neighbors == 1 )
        gs = gridsearch_create_nn(nvals, lons, lats);
      else
        gs = gridsearch_create(nvals, lons, lats);
    }
  
  finish = clock();

  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);

  progressInit();

  start = clock();

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nbr_mask, nbr_add, nbr_dist)  \
  shared(findex, mindex, vindex, array1, array2, xvals, yvals, gs, nmiss, num_neighbors)
#endif
  for ( unsigned i = 0; i < nmiss; ++i )
    {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/nmiss);

      int ompthID = cdo_omp_get_thread_num();

      grid_search_nbr(gs, num_neighbors, nbr_add[ompthID], nbr_dist[ompthID], xvals[mindex[i]], yvals[mindex[i]]);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */
      double dist_tot = nbr_compute_weights(num_neighbors, NULL, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      /* Normalize weights and store the link */
      size_t nadds = nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);
      if ( nadds )
        {
          double result = 0;
          for ( size_t n = 0; n < nadds; ++n ) result += array1[vindex[nbr_add[ompthID][n]]]*nbr_dist[ompthID][n];
          array2[mindex[i]] = result;
        }
    }

  if ( mindex ) Free(mindex);
  if ( vindex ) Free(vindex);

  finish = clock();

  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

  if ( gridID0 != gridID ) gridDestroy(gridID);

  if ( lons ) Free(lons);
  if ( lats ) Free(lats);
  Free(xvals);
  Free(yvals);
}


void *Fillmiss(void *argument)
{
  int nmiss;
  int nrecs, varID, levelID;
  void (*fill_method) (field_type *fin , field_type *fout , int) = NULL;

  cdoInitialize(argument);

  // clang-format off
  int FILLMISS        = cdoOperatorAdd("fillmiss"   ,   0, 0, "nfill");
  int FILLMISSONESTEP = cdoOperatorAdd("fillmiss2"  ,   0, 0, "nfill");
  int SETMISSTONN     = cdoOperatorAdd("setmisstonn" ,  0, 0, "");
  int SETMISSTODIS    = cdoOperatorAdd("setmisstodis" , 0, 0, "number of neighbors");
  // clang-format on

  int operatorID      = cdoOperatorID();

  int nfill = 1;
  if ( operatorID == FILLMISS )
     {
       fill_method = &fillmiss;
     }
  else if ( operatorID == FILLMISSONESTEP )
     {
       fill_method = &fillmiss_one_step;
     }
  else if ( operatorID == SETMISSTONN )
     {
       fill_method = &setmisstodis;
     }
  else if ( operatorID == SETMISSTODIS )
     {
       nfill = 4;
       fill_method = &setmisstodis;
     }

  /* Argument handling */
  {
    int oargc = operatorArgc();
    char **oargv = operatorArgv();

    if ( oargc == 1 )
      {
        nfill = parameter2int(oargv[0]);
        if ( operatorID == FILLMISS ) 
          {
            if ( nfill < 1 || nfill > 4 ) cdoAbort("nfill out of range!");
          }
      }
    else if ( oargc > 1 )
      cdoAbort("Too many arguments!");
  }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));

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

	  pstreamDefRecord(streamID2, varID, levelID);

          if ( field1.nmiss == 0 )
            {
              pstreamWriteRecord(streamID2, field1.ptr, 0);
            }
          else
            {
              int gridID = vlistInqVarGrid(vlistID1, varID);

              if ( (operatorID == FILLMISS || operatorID == FILLMISSONESTEP) && 
                   (gridInqType(gridID) == GRID_GME || gridInqType(gridID) == GRID_UNSTRUCTURED) )
                cdoAbort("%s data unsupported!", gridNamePtr(gridInqType(gridID)) );
                
              field1.grid    = gridID;
              field1.missval = vlistInqVarMissval(vlistID1, varID);

              field2.grid    = field1.grid;
              field2.nmiss   = 0;
              field2.missval = field1.missval;

              fill_method(&field1, &field2, nfill);

              int gridsize = gridInqSize(field2.grid);
              int nmiss = 0;
              for ( int i = 0; i < gridsize; ++i )
                if ( DBL_IS_EQUAL(field2.ptr[i], field2.missval) ) nmiss++;
              
              pstreamWriteRecord(streamID2, field2.ptr, nmiss);
            }
        }
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( field2.ptr ) Free(field2.ptr);
  if ( field1.ptr ) Free(field1.ptr);

  cdoFinish();

  return 0;
}
