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

      Selindex     selindex    Select grid indices
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "listarray.h"
#include "pstream.h"


int gengridcell(int gridID1, int gridsize2, int *cellidx);

static
int genindexgrid(int gridID1, int gridsize2, int *cellidx)
{
  int gridID0 = gridID1;
  int gridtype1 = gridInqType(gridID1);

  if ( gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN || gridtype1 == GRID_PROJECTION )
    {
      gridID1 = gridToCurvilinear(gridID1, 0);
      gridtype1 = GRID_CURVILINEAR;
    }

  int gridID2 = -1;
  if ( gridtype1 == GRID_UNSTRUCTURED || gridtype1 == GRID_CURVILINEAR )
    gridID2 = gengridcell(gridID1, gridsize2, cellidx);

  if ( gridID0 != gridID1 ) gridDestroy(gridID1);

  return gridID2;
}

static
void sel_index(double *array1, double *array2, int nind, int *indarr)
{
  for ( int i = 0; i < nind; ++i )
    {
      array2[i] = array1[indarr[i]];
    }
}


void *Selgridcell(void *argument)
{
  int nrecs;
  int varID;
  int gridID1 = -1, gridID2;
  int index, gridtype = -1;
  typedef struct {
    int gridID1, gridID2;
  } sindex_t;
  lista_t *ilista = lista_new(INT_LISTA);

  cdoInitialize(argument);

                    cdoOperatorAdd("selgridcell", 0, 0, "grid cell indices (1-N)");
  int DELGRIDCELL = cdoOperatorAdd("delgridcell", 0, 0, "grid cell indices (1-N)");

  operatorInputArg(cdoOperatorEnter(0));

  int operatorID = cdoOperatorID();

  int nind;
  int *indarr = NULL;
  if ( operatorArgc() == 1 && fileExists(operatorArgv()[0]) )
    {
      bool *cdo_read_mask(const char *maskfile, int *n);
      int n = 0;
      bool *mask = cdo_read_mask(operatorArgv()[0], &n);
      nind = 0;
      for ( int i = 0; i < n; ++i ) if ( mask[i] ) nind++;
      if ( nind == 0 ) cdoAbort("Mask is empty!");
      else
        {
          indarr = (int*) Malloc(nind*sizeof(double));
          nind = 0;
          for ( int i = 0; i < n; ++i ) if ( mask[i] ) indarr[nind++] = i;
        }
      if ( mask ) Free(mask);
    }
 else
    {
      nind = args2int_lista(operatorArgc(), operatorArgv(), ilista);
      indarr = (int*) lista_dataptr(ilista);

      if ( cdoVerbose )
        for ( int i = 0; i < nind; i++ )
          cdoPrint("int %d = %d", i+1, indarr[i]);

      for ( int i = 0; i < nind; i++ ) indarr[i] -= 1;
    }

  int indmin = indarr[0];
  int indmax = indarr[0];
  for ( int i = 1; i < nind; i++ )
    {
      if ( indmax < indarr[i] ) indmax = indarr[i];
      if ( indmin > indarr[i] ) indmin = indarr[i];
    }

  if ( indmin < 0 ) cdoAbort("Index < 1 not allowed!");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  bool *vars = (bool*) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;

  int ngrids = vlistNgrids(vlistID1);
  sindex_t *sindex = (sindex_t *) Malloc(ngrids*sizeof(sindex_t));

  int ncells = nind;
  int *cellidx = indarr;
  if ( operatorID == DELGRIDCELL )
    {
      int gridsize = vlistGridsizeMax(vlistID1);
      ncells = gridsize - nind;
      cellidx = (int*) Malloc(gridsize*sizeof(int));
      for ( int i = 0; i < gridsize; ++i ) cellidx[i] = 1;
      for ( int i = 0; i < nind; ++i ) cellidx[indarr[i]] = 0;
      int j = 0;
      for ( int i = 0; i < gridsize; ++i )
        if ( cellidx[i] == 1 ) cellidx[j++] = i;
      if ( j != ncells ) cdoAbort("Internal error; number of cells differ");
    }
  if ( ncells == 0 ) cdoAbort("Mask is empty!");

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      int gridsize = gridInqSize(gridID1);
      if ( gridsize == 1 ) continue;
      if ( indmax >= gridsize )
        {
          cdoWarning("Max grid index is greater than grid size, skipped grid %d!", index+1);
          continue;
        }

      gridID2 = genindexgrid(gridID1, ncells, cellidx);

      if ( gridID2 == -1 )
        {
          cdoWarning("Unsupported grid type >%s<, skipped grid %d!", gridNamePtr(gridtype), index+1);
          continue;
        }
	  
      sindex[index].gridID1 = gridID1;
      sindex[index].gridID2 = gridID2;

      vlistChangeGridIndex(vlistID2, index, gridID2);

      for ( varID = 0; varID < nvars; varID++ )
        if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
          vars[varID] = true;
    }

  for ( varID = 0; varID < nvars; varID++ )
    if ( vars[varID] ) break;

  if ( varID >= nvars ) cdoAbort("No variables selected!");

  
  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  int gridsize2 = vlistGridsizeMax(vlistID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double*) Malloc(gridsize2*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
          int nmiss, varID, levelID;
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss);

	  pstreamDefRecord(streamID2, varID, levelID);

	  if ( vars[varID] )
	    {	      
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      for ( index = 0; index < ngrids; index++ )
		if ( gridID1 == sindex[index].gridID1 ) break;

	      if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

	      gridsize2 = gridInqSize(sindex[index].gridID2);

              sel_index(array1, array2, ncells, cellidx);

	      if ( nmiss )
		{
		  nmiss = 0;
                  double missval = vlistInqVarMissval(vlistID2, varID);
		  for ( int i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      pstreamWriteRecord(streamID2, array2, nmiss);
	    }
	  else
	    {
	      pstreamWriteRecord(streamID2, array1, nmiss);
	    }
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  if ( vars   ) Free(vars);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);
  if ( sindex ) Free(sindex);

  lista_destroy(ilista);

  if ( operatorID == DELGRIDCELL ) Free(cellidx);

  cdoFinish();

  return 0;
}
