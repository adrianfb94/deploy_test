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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"

#define  MAX_BLOCKS  65536

static
void genGrids(int gridID1, int *gridIDs, int nxvals, int nyvals, int nxblocks, int nyblocks,
	      int **gridindex, int *ogridsize, int nsplit)
{
  double *xpvals = NULL, *ypvals = NULL;
  int gridtype = gridInqType(gridID1);
  bool lregular = true;
  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC )
    lregular = true;
  else if ( gridtype == GRID_CURVILINEAR )
    lregular = false;
  else
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);

  bool lxcoord = gridInqXvals(gridID1, NULL) > 0;
  bool lycoord = gridInqYvals(gridID1, NULL) > 0;

  double *xvals = NULL, *yvals = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;

  if ( lxcoord )
    {
      if ( lregular )
        {
          xvals = (double*) Malloc(nx*sizeof(double));
        }
      else
        {
          xvals = (double*) Malloc(nx*ny*sizeof(double));
          xvals2 = (double*) Malloc(nxvals*nyvals*sizeof(double));
        }
      
      gridInqXvals(gridID1, xvals);
    }

  if ( lycoord )
    {
      if ( lregular )
        {
          yvals = (double*) Malloc(ny*sizeof(double));
        }
      else
        {
          yvals = (double*) Malloc(nx*ny*sizeof(double));
          yvals2 = (double*) Malloc(nxvals*nyvals*sizeof(double));
        }

      gridInqYvals(gridID1, yvals);
    }
 
  int *xlsize = (int*) Malloc(nxblocks*sizeof(int));
  int *ylsize = (int*) Malloc(nyblocks*sizeof(int));

  for ( int ix = 0; ix < nxblocks; ++ix ) xlsize[ix] = nxvals;
  if ( nx%nxblocks != 0 ) xlsize[nxblocks-1] = nx - (nxblocks-1)*nxvals;
  if ( cdoVerbose ) for ( int ix = 0; ix < nxblocks; ++ix ) cdoPrint("xblock %d: %d", ix, xlsize[ix]);

  for ( int iy = 0; iy < nyblocks; ++iy ) ylsize[iy] = nyvals;
  if ( ny%nyblocks != 0 ) ylsize[nyblocks-1] = ny - (nyblocks-1)*nyvals;
  if ( cdoVerbose ) for ( int iy = 0; iy < nyblocks; ++iy ) cdoPrint("yblock %d: %d", iy, ylsize[iy]);

  int index = 0;
  for ( int iy = 0; iy < nyblocks; ++iy )
    for ( int ix = 0; ix < nxblocks; ++ix )
      {
	int offset = iy*nyvals*nx + ix*nxvals;

	int gridsize2 = xlsize[ix]*ylsize[iy];
	gridindex[index] = (int*) Malloc(gridsize2*sizeof(int));

	gridsize2 = 0;
        // printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
	for ( int j = 0; j < ylsize[iy]; ++j )
          for ( int i = 0; i < xlsize[ix]; ++i )
            {
              // printf(">> %d %d %d\n", j, i, offset + j*nx + i);
              if ( !lregular )
                {
                  if ( lxcoord ) xvals2[gridsize2] = xvals[offset + j*nx + i];
                  if ( lycoord ) yvals2[gridsize2] = yvals[offset + j*nx + i];
                }
              gridindex[index][gridsize2++] = offset + j*nx + i;
            }
	// printf("gridsize2 %d\n", gridsize2);

	int gridID2 = gridCreate(gridtype, gridsize2);
	gridDefXsize(gridID2, xlsize[ix]);
	gridDefYsize(gridID2, ylsize[iy]);

        gridDefNP(gridID2, gridInqNP(gridID1));
        gridDefDatatype(gridID2, gridInqDatatype(gridID1));

        grid_copy_attributes(gridID1, gridID2);

        if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridID1, gridID2);

        if ( lregular )
          {
            if ( lxcoord ) gridDefXvals(gridID2, xvals+ix*nxvals);
            if ( lycoord ) gridDefYvals(gridID2, yvals+iy*nyvals);
          }
        else
          {
            if ( lxcoord ) gridDefXvals(gridID2, xvals2);
            if ( lycoord ) gridDefYvals(gridID2, yvals2);
          }

        int projID1 = gridInqProj(gridID1);
        if ( projID1 != CDI_UNDEFID && gridInqType(projID1) == GRID_PROJECTION )
          {
            int projID2 = gridCreate(GRID_PROJECTION, gridsize2);
            gridDefXsize(projID2, xlsize[ix]);
            gridDefYsize(projID2, ylsize[iy]);

            grid_copy_attributes(projID1, projID2);
            grid_copy_mapping(projID1, projID2);

            bool lxpcoord = gridInqXvals(projID1, NULL) > 0;
            if ( lxpcoord )
              {
                if ( !xpvals )
                  {
                    xpvals = (double*) Malloc(nx*sizeof(double));
                    gridInqXvals(projID1, xpvals);
                  }
                gridDefXvals(projID2, xpvals+ix*nxvals);
              }
            bool lypcoord = gridInqYvals(projID1, NULL) > 0;
            if ( lypcoord )
              {
                if ( !ypvals )
                  {
                    ypvals = (double*) Malloc(ny*sizeof(double));
                    gridInqYvals(projID1, ypvals);
                  }
                gridDefYvals(projID2, ypvals+iy*nyvals);
              }

            gridDefProj(gridID2, projID2);
          }

	gridIDs[index] = gridID2;
	ogridsize[index] = gridsize2;

	index++;
	if ( index > nsplit )
	  cdoAbort("Internal problem, index exceeded bounds!");
      }

  if ( xvals2 ) Free(xvals2);
  if ( yvals2 ) Free(yvals2);
  if ( xvals ) Free(xvals);
  if ( yvals ) Free(yvals);
  if ( xpvals ) Free(xpvals);
  if ( ypvals ) Free(ypvals);
  if ( xlsize ) Free(xlsize);
  if ( ylsize ) Free(ylsize);
}

static
void window_cell(double *array1, double *array2, long gridsize2, int *cellidx)
{
  for ( long i = 0; i < gridsize2; ++i )
    array2[i] = array1[cellidx[i]];
}

typedef struct
{
  int gridID;
  int *gridIDs;
  int *gridsize;
  int **gridindex;
} sgrid_t;


void *Distgrid(void *argument)
{
  int gridID1;
  int varID, levelID;
  int nrecs;
  char filesuffix[32];
  char filename[8192];
  int index;
  int gridtype = -1;
  int nmiss;
  int i;

  cdoInitialize(argument);

  operatorInputArg("nxblocks, [nyblocks]");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
  if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");
  int nxblocks = parameter2int(operatorArgv()[0]);
  int nyblocks = 1;
  if ( operatorArgc() == 2 ) nyblocks = parameter2int(operatorArgv()[1]);

  if ( nxblocks <= 0 ) cdoAbort("nxblocks has to be greater than 0!");
  if ( nyblocks <= 0 ) cdoAbort("nyblocks has to be greater than 0!");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT   || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
	  (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	   break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian, curvilinear or generic grid found (%s data unsupported)!", gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  int gridsize = gridInqSize(gridID1);
  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);
  for ( int i = 1; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      if ( gridsize != gridInqSize(gridID1) )
	cdoAbort("Gridsize must not change!");
    }

  if ( nxblocks > nx )
    {
      cdoPrint("nxblocks (%d) greater than nx (%d), set to %d!", nxblocks, nx, nx);
      nxblocks = nx;
    }
  if ( nyblocks > ny )
    {
      cdoPrint("nyblocks (%d) greater than ny (%d), set to %d!", nyblocks, ny, ny);
      nyblocks = ny;
    }

  int xinc = nx/nxblocks;
  int yinc = ny/nyblocks;

  if ( nx%xinc && nx%(xinc+1) ) xinc++;
  if ( ny%yinc && ny%(yinc+1) ) yinc++;

  int nsplit = nxblocks*nyblocks;
  if ( nsplit > MAX_BLOCKS ) cdoAbort("Too many blocks (max = %d)!", MAX_BLOCKS);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  int *vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
  int *streamIDs = (int*) Malloc(nsplit*sizeof(int));

  sgrid_t *grids = (sgrid_t*) Malloc(ngrids*sizeof(sgrid_t));
  for ( int i = 0; i < ngrids; i++ )
    {  
      grids[i].gridID    = vlistGrid(vlistID1, i);
      grids[i].gridIDs   = (int*) Malloc(nsplit*sizeof(int));
      grids[i].gridsize  = (int*) Malloc(nsplit*sizeof(int));
      grids[i].gridindex = (int**) Malloc(nsplit*sizeof(int*));

      for ( int index = 0; index < nsplit; index++ ) grids[i].gridindex[index] = NULL;
    }

  for ( int index = 0; index < nsplit; index++ )
    vlistIDs[index] = vlistDuplicate(vlistID1);

  if ( cdoVerbose ) cdoPrint("ngrids=%d  nsplit=%d", ngrids, nsplit);

  for ( int i = 0; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, grids[i].gridIDs, xinc, yinc, nxblocks, nyblocks, grids[i].gridindex, grids[i].gridsize, nsplit);
      /*
      if ( cdoVerbose )
	for ( int index = 0; index < nsplit; index++ )
	  cdoPrint("Block %d,  gridID %d,  gridsize %d", index+1, grids[i].gridIDs[index], gridInqSize(grids[i].gridIDs[index]));
      */
      for ( int index = 0; index < nsplit; index++ )
	vlistChangeGridIndex(vlistIDs[index], i, grids[i].gridIDs[index]);
    }

  int gridsize2max = 0;
  for ( int index = 0; index < nsplit; index++ )
    if ( grids[0].gridsize[index] > gridsize2max ) gridsize2max = grids[0].gridsize[index];

  double *array2 = (double*) Malloc(gridsize2max*sizeof(double));

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  const char *refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  for ( int index = 0; index < nsplit; index++ )
    {
      sprintf(filename+nchars, "%05d", index);
      if ( filesuffix[0] )
	sprintf(filename+nchars+5, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      pstreamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  if ( ngrids > 1 ) cdoPrint("Baustelle: number of different grids > 1!");
  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      for ( int index = 0; index < nsplit; index++ )
	pstreamDefTimestep(streamIDs[index], tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss);

	  double missval = vlistInqVarMissval(vlistID1, varID);

	  for ( int index = 0; index < nsplit; index++ )
	    {
	      i = 0;
	      window_cell(array1, array2, grids[i].gridsize[index], grids[i].gridindex[index]);
	      pstreamDefRecord(streamIDs[index], varID, levelID);
	      if ( nmiss > 0 )
		{
		  nmiss = 0;
		  for ( int k = 0; k < grids[i].gridsize[index]; ++k )
		    if ( DBL_IS_EQUAL(array2[k], missval) ) nmiss++;
		}
	      pstreamWriteRecord(streamIDs[index], array2, nmiss);
	    }
	}

      tsID++;
    }

  pstreamClose(streamID1);

  for ( int index = 0; index < nsplit; index++ )
    {
      pstreamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);

  if ( vlistIDs  ) Free(vlistIDs);
  if ( streamIDs ) Free(streamIDs);

  for ( int i = 0; i < ngrids; i++ )
    {
      for ( int index = 0; index < nsplit; index++ )
	gridDestroy(grids[i].gridIDs[index]);
      Free(grids[i].gridIDs);
      Free(grids[i].gridsize);

      for ( int index = 0; index < nsplit; index++ )
        Free(grids[i].gridindex[index]);
      Free(grids[i].gridindex);
    }
  Free(grids);

  cdoFinish();

  return 0;
}
