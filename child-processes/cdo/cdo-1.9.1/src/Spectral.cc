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

      Spectral   sp2gp           Spectral to gridpoint
      Spectral   sp2gpl          Spectral to gridpoint linear
      Spectral   gp2sp           Gridpoint to spectral
      Spectral   gp2spl          Gridpoint to spectral linear
      Spectral   sp2sp           Spectral to spectral
      Spectral   spcut           Cut spectral wave number
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"
#include "specspace.h"
#include "listarray.h"


void *Spectral(void *argument)
{
  int nrecs, nvars;
  int varID, levelID;
  int index;
  int gridIDsp = -1, gridIDgp = -1;
  int gridID1 = -1, gridID2 = -1;
  int gridID;
  int nmiss;
  int ncut = 0;
  int *wnums = NULL, *waves = NULL;
  int *vars;
  double *array2 = NULL;
  int nlon, nlat, ntr;
  SPTRANS *sptrans = NULL;
  lista_t *ilista = lista_new(INT_LISTA);

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int GP2SP  = cdoOperatorAdd("gp2sp",  0, 0, NULL);
  int GP2SPL = cdoOperatorAdd("gp2spl", 0, 0, NULL);
  int SP2GP  = cdoOperatorAdd("sp2gp",  0, 0, NULL);
  int SP2GPL = cdoOperatorAdd("sp2gpl", 0, 0, NULL);
  int SP2SP  = cdoOperatorAdd("sp2sp",  0, 0, NULL);
  int SPCUT  = cdoOperatorAdd("spcut",  0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  /* find first spectral grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  gridIDsp = gridID;
	  break;
	}
    }
  /* find first gaussian grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_GAUSSIAN )
	{
	  gridIDgp = gridID;
	  break;
	}
    }

  /* define output grid */
  if ( operatorID == GP2SP || operatorID == GP2SPL )
    {
      if ( gridIDgp == -1 ) cdoWarning("No data on Gaussian grid found!");

      gridID1 = gridIDgp;

      if ( gridID1 != -1 )
	{
	  if ( operatorID == GP2SP )
	    ntr = nlat_to_ntr(gridInqYsize(gridID1));
	  else
	    ntr = nlat_to_ntr_linear(gridInqYsize(gridID1));

	  if ( gridIDsp != -1 )
	    if ( ntr != gridInqTrunc(gridIDsp) ) gridIDsp = -1;

	  if ( gridIDsp == -1 )
	    {
	      gridIDsp = gridCreate(GRID_SPECTRAL, (ntr+1)*(ntr+2));
	      gridDefTrunc(gridIDsp, ntr);
	      gridDefComplexPacking(gridIDsp, 1);
	    }

	  if ( gridIDsp == -1 && gridInqType(vlistGrid(vlistID1, 0)) == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

	  if ( gridIDsp == -1 ) cdoAbort("Computation of spherical harmonics failed!");

	  gridID2 = gridIDsp;

	  nlon = gridInqXsize(gridID1);
	  nlat = gridInqYsize(gridID1);
	  ntr  = gridInqTrunc(gridID2);

	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	}
    }
  else if ( operatorID == SP2GP || operatorID == SP2GPL )
    {   
      if ( gridIDsp == -1 ) cdoWarning("No spectral data found!");

      gridID1 = gridIDsp;

      if ( gridID1 != -1 )
	{
	  if ( gridIDgp != -1 )
	    {
	      if ( operatorID == SP2GP )
		ntr = nlat_to_ntr(gridInqYsize(gridIDgp));
	      else
		ntr = nlat_to_ntr_linear(gridInqYsize(gridIDgp));

	      if ( gridInqTrunc(gridIDsp) != ntr ) gridIDgp = -1;
	    }

	  if ( gridIDgp == -1 )
	    {
	      char gridname[20];
	      if ( operatorID == SP2GP )
		snprintf(gridname, sizeof(gridname), "t%dgrid", gridInqTrunc(gridIDsp));
	      else
		snprintf(gridname, sizeof(gridname), "tl%dgrid", gridInqTrunc(gridIDsp));

	      gridIDgp = grid_from_name(gridname);
	    }

	  gridID2 = gridIDgp;

	  ntr  = gridInqTrunc(gridID1);
	  nlon = gridInqXsize(gridID2);
	  nlat = gridInqYsize(gridID2);
      
	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	}
    }
  else if ( operatorID == SP2SP )
    {
      gridID1 = gridIDsp;

      operatorInputArg("truncation");
      if ( gridID1 != -1 )
	{
	  if ( !isdigit(operatorArgv()[0][0]) ) cdoAbort("parameter truncation must comprise only digits [0-9]!");
	  int ntr = parameter2int(operatorArgv()[0]);
	  int nsp = (ntr+1)*(ntr+2);
	  gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
	  gridDefTrunc(gridIDsp, ntr);
	  gridDefComplexPacking(gridIDsp, 1);
	}
      else
	cdoAbort("No spectral data found!");

      gridID2 = gridIDsp;
    }
  else if ( operatorID == SPCUT )
    {
      long i, j, maxntr;
      gridID1 = gridIDsp;

      operatorInputArg("wave numbers");
      if ( gridID1 != -1 )
	{
	  maxntr = 1+gridInqTrunc(gridID1);
	  ncut = args2int_lista(operatorArgc(), operatorArgv(), ilista);
	  wnums = (int *) lista_dataptr(ilista);
	  waves = (int*) Malloc(maxntr*sizeof(int));
	  for ( i = 0; i < maxntr; i++ ) waves[i] = 1;
	  for ( i = 0; i < ncut; i++ )
	    {
	      j = wnums[i] - 1;
	      if ( j < 0 || j >= maxntr )
		cdoAbort("wave number %d out of range (min=1, max=%d)!", wnums[i], maxntr);
	      waves[j] = 0;
	    }
	}
      else
	cdoAbort("No spectral data found!");

      gridID2 = gridIDsp;
    }

  nvars = vlistNvars(vlistID2);
  vars  = (int*) Malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  if ( gridID1 != -1 ) vlistChangeGrid(vlistID2, gridID1, gridID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  if ( gridID2 != -1 )
    {
      gridsize = gridInqSize(gridID2);
      array2 = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      pstreamReadRecord(streamID1, array1, &nmiss);
	      if ( nmiss ) cdoAbort("Missing values unsupported for spectral data!");

	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      if ( operatorID == GP2SP || operatorID == GP2SPL )
		grid2spec(sptrans, gridID1, array1, gridID2, array2);	      
	      else if ( operatorID == SP2GP || operatorID == SP2GPL )
		spec2grid(sptrans, gridID1, array1, gridID2, array2);
	      else if ( operatorID == SP2SP )
		spec2spec(gridID1, array1, gridID2, array2);
	      else if ( operatorID == SPCUT )
		speccut(gridID1, array1, array2, waves);

	      pstreamDefRecord(streamID2, varID, levelID);
	      pstreamWriteRecord(streamID2, array2, nmiss);  
	    }   
	  else
	    {
	      pstreamDefRecord(streamID2, varID, levelID);
	      if ( lcopy )
		{
		  pstreamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  pstreamReadRecord(streamID1, array1, &nmiss);
		  pstreamWriteRecord(streamID2, array1, nmiss);
		}
	    }    
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);
  if ( vars )   Free(vars);
  if ( waves )  Free(waves);

  lista_destroy(ilista);

  sptrans_delete(sptrans);

  cdoFinish();

  return 0;
}
