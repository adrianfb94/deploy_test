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
#include "pstream.h"
#include "grid.h"
#include "util.h"


typedef struct
{
  int streamID;
  int vlistID;
  int gridID;
  int nmiss;
  int gridsize;
  int *gridindex;
  double *array;
} ens_file_t;


typedef struct
{
  double x, y;
  int id;
}
xyinfo_t;

static
int cmpx(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->x < xy2->x ) cmp = -1;
  else if ( xy1->x > xy2->x ) cmp =  1;

  return cmp;
}

static
int cmpxy_lt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return cmp;
}

static
int cmpxy_gt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return cmp;
}

static
int genGrid(int ngrids, int nfiles, ens_file_t *ef, bool ginit, int igrid, int nxblocks)
{
  bool lsouthnorth = true;
  bool lregular = false;
  bool lcurvilinear = false;
  int gridID2 = -1;
  double *xvals2 = NULL, *yvals2 = NULL;

  int nx = -1;
  if ( nxblocks != -1 ) nx = nxblocks;

  int gridID   = vlistGrid(ef[0].vlistID, igrid);
  int gridtype = gridInqType(gridID);
  if ( ngrids > 1 && gridtype == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
    return gridID2;

  int *xsize = (int*) Malloc(nfiles*sizeof(int));
  int *ysize = (int*) Malloc(nfiles*sizeof(int));
  xyinfo_t *xyinfo = (xyinfo_t*) Malloc(nfiles*sizeof(xyinfo_t));
  double **xvals = (double**) Malloc(nfiles*sizeof(double*));
  double **yvals = (double**) Malloc(nfiles*sizeof(double*));

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      gridID   = vlistGrid(ef[fileID].vlistID, igrid);
      int gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
        lregular = true;
      else if ( gridtype == GRID_CURVILINEAR )
        lcurvilinear = true;
      else if ( gridtype == GRID_GENERIC /*&& gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0*/ )
        ;
      else
	cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      xsize[fileID] = gridInqXsize(gridID);
      ysize[fileID] = gridInqYsize(gridID);
      if ( xsize[fileID] == 0 ) xsize[fileID] = 1;
      if ( ysize[fileID] == 0 ) ysize[fileID] = 1;

      if ( lregular )
        {
          xvals[fileID] = (double*) Malloc(xsize[fileID]*sizeof(double));
          yvals[fileID] = (double*) Malloc(ysize[fileID]*sizeof(double));
        }
      else if ( lcurvilinear )
        {
          xvals[fileID] = (double*) Malloc(xsize[fileID]*ysize[fileID]*sizeof(double));
          yvals[fileID] = (double*) Malloc(xsize[fileID]*ysize[fileID]*sizeof(double));
        }
      else
        {
          xvals[fileID] = NULL;
          yvals[fileID] = NULL;
        }
        
      if ( lregular || lcurvilinear )
        {
          gridInqXvals(gridID, xvals[fileID]);
          gridInqYvals(gridID, yvals[fileID]);
        }
      // printf("fileID %d, gridID %d\n", fileID, gridID);

      if ( lregular )
        {
          xyinfo[fileID].x  = xvals[fileID][0];
          xyinfo[fileID].y  = yvals[fileID][0];
          xyinfo[fileID].id = fileID;

          if ( ysize[fileID] > 1 )
            {
              if ( yvals[fileID][0] > yvals[fileID][ysize[fileID]-1] ) lsouthnorth = false;
            }
        }
      else
        {
          xyinfo[fileID].x  = 0;
          xyinfo[fileID].y  = 0;
          xyinfo[fileID].id = fileID;
        }
    }

  if ( cdoVerbose && lregular )
    for ( int fileID = 0; fileID < nfiles; fileID++ )
      printf("1 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

  if ( lregular )
    {
      qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpx);  	      

      if ( cdoVerbose )
        for ( int fileID = 0; fileID < nfiles; fileID++ )
          printf("2 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      if ( lsouthnorth )
        qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_lt);  
      else
        qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_gt);  	      

      if ( cdoVerbose )
        for ( int fileID = 0; fileID < nfiles; fileID++ )
          printf("3 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      if ( nx <= 0 )
        {
          nx = 1;
          for ( int fileID = 1; fileID < nfiles; fileID++ )
            {
              if ( DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y) ) nx++;
              else break;
            }
        }
    }
  else
    {
      if ( nx <= 0 ) nx = nfiles;
    }

  int ny = nfiles/nx;
  if ( nx*ny != nfiles ) cdoAbort("Number of input files (%d) and number of blocks (%dx%d) differ!", nfiles, nx, ny);
 
  int xsize2 = 0;
  for ( int i = 0; i < nx; ++i ) xsize2 += xsize[xyinfo[i].id];
  int ysize2 = 0;
  for ( int j = 0; j < ny; ++j ) ysize2 += ysize[xyinfo[j*nx].id];
  if ( cdoVerbose ) cdoPrint("xsize2 %d  ysize2 %d", xsize2, ysize2);

  if ( lregular )
    {
      xvals2 = (double*) Malloc(xsize2*sizeof(double));
      yvals2 = (double*) Malloc(ysize2*sizeof(double));
    }
  else if ( lcurvilinear )
    {
      xvals2 = (double*) Malloc(xsize2*ysize2*sizeof(double));
      yvals2 = (double*) Malloc(xsize2*ysize2*sizeof(double));
    }
    
  int *xoff = (int*) Malloc((nx+1)*sizeof(int));
  int *yoff = (int*) Malloc((ny+1)*sizeof(int));

  xoff[0] = 0;
  for ( int i = 0; i < nx; ++i )
    {
      int idx = xyinfo[i].id;
      if ( lregular ) memcpy(xvals2+xoff[i], xvals[idx], xsize[idx]*sizeof(double));
      xoff[i+1] = xoff[i] + xsize[idx];
    }

  yoff[0] = 0;
  for ( int j = 0; j < ny; ++j )
    {
      int idx = xyinfo[j*nx].id;
      if ( lregular ) memcpy(yvals2+yoff[j], yvals[idx], ysize[idx]*sizeof(double));
      yoff[j+1] = yoff[j] + ysize[idx];
    }

  if ( ginit == false )
    {
      for ( int fileID = 0; fileID < nfiles; fileID++ )
	{
	  int idx = xyinfo[fileID].id;
	  int iy = fileID/nx;
	  int ix = fileID - iy*nx;

          int offset = yoff[iy]*xsize2 + xoff[ix];
	  /*
	  printf("fileID %d %d, iy %d, ix %d, offset %d\n",
		 fileID, xyinfo[fileID].id, iy, ix, offset);
	  */
	  int ij = 0;
	  for ( int j = 0; j < ysize[idx]; ++j )
	    for ( int i = 0; i < xsize[idx]; ++i )
	      {
                if ( lcurvilinear )
                  {
                    xvals2[offset+j*xsize2+i] = xvals[idx][ij];
                    yvals2[offset+j*xsize2+i] = yvals[idx][ij];
                  }
                ef[idx].gridindex[ij++] = offset+j*xsize2+i;
	      }
	}
    }

  gridID2 = gridCreate(gridtype, xsize2*ysize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);
  if ( lregular || lcurvilinear )
    {
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);
    }

  Free(xoff);
  Free(yoff);
  Free(xsize);
  Free(ysize);
  if ( xvals2 ) Free(xvals2);
  if ( yvals2 ) Free(yvals2);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( xvals[fileID] ) Free(xvals[fileID]);
      if ( yvals[fileID] ) Free(yvals[fileID]);
    }
  Free(xvals);
  Free(yvals);
  Free(xyinfo);

  gridID = vlistGrid(ef[0].vlistID, igrid);

  grid_copy_attributes(gridID, gridID2);

  if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridID, gridID2);

  return gridID2;
}


void *Collgrid(void *argument)
{
  int nxblocks = -1;
  int varID, levelID;
  int nrecs0;

  cdoInitialize(argument);
    
  int nfiles = cdoStreamCnt() - 1;
  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename) )
    cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t*) Malloc(nfiles*sizeof(ens_file_t));

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].streamID = pstreamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID  = pstreamInqVlist(ef[fileID].streamID);
    }

  int vlistID1 = ef[0].vlistID;
  vlistClearFlag(vlistID1);

  /* check that the contents is always the same */
  for ( int fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(vlistID1, ef[fileID].vlistID, CMP_NAME | CMP_NLEVEL);

  int nvars = vlistNvars(vlistID1);
  bool *vars = (bool*) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;
  bool *vars1 = (bool*) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars1[varID] = false;

  int nsel = operatorArgc();
  int noff = 0;

  if ( nsel > 0 )
    {
      int len = (int) strlen(operatorArgv()[0]);
      while ( --len >= 0 && isdigit(operatorArgv()[0][len]) ) ;

      if ( len == -1 )
        {
          nsel--;
          noff++;
          nxblocks = parameter2int(operatorArgv()[0]);
        }
    }

  if ( nsel == 0 )
    {
      for ( varID = 0; varID < nvars; varID++ ) vars1[varID] = true;
    }
  else
    {
      char **argnames = operatorArgv() + noff;

      if ( cdoVerbose )
	for ( int i = 0; i < nsel; i++ )
	  fprintf(stderr, "name %d = %s\n", i+1, argnames[i]);

      bool *selfound = (bool*) Malloc(nsel*sizeof(bool));
      for ( int i = 0; i < nsel; i++ ) selfound[i] = false;

      char varname[CDI_MAX_NAME];
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  for ( int isel = 0; isel < nsel; isel++ )
	    {
	      if ( strcmp(argnames[isel], varname) == 0 )
		{
		  selfound[isel] = true;
		  vars1[varID] = true;
		}
	    }
	}

      for ( int isel = 0; isel < nsel; isel++ )
	if ( selfound[isel] == false )
	  cdoAbort("Variable name %s not found!", argnames[isel]);

      Free(selfound);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vars1[varID] )
	{
	  int zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  int nlevs    = zaxisInqSize(zaxisID);
	  for ( int levID = 0; levID < nlevs; levID++ )
	    vlistDefFlag(vlistID1, varID, levID, TRUE);
	}
    }

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      int gridsize = vlistGridsizeMax(ef[fileID].vlistID);
      ef[fileID].gridsize = gridsize;
      ef[fileID].gridindex = (int*) Malloc(gridsize*sizeof(int));
      ef[fileID].array = (double*) Malloc(gridsize*sizeof(double));
    }

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);
  /*
  if ( cdoVerbose )
    {
      vlistPrint(vlistID1);
      vlistPrint(vlistID2);
    }
  */
  //int vlistID2 = vlistDuplicate(vlistID1);
  int nvars2 = vlistNvars(vlistID2);
  // int *vars  = (int*) Malloc(nvars*sizeof(int));
  //for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;

  int ngrids1 = vlistNgrids(vlistID1);
  int ngrids2 = vlistNgrids(vlistID2);

  int *gridIDs = (int*) Malloc(ngrids2*sizeof(int));

  bool ginit = false;
  for ( int i2 = 0; i2 < ngrids2; ++i2 )
    {
      int i1;
      for ( i1 = 0; i1 < ngrids1; ++i1 )
	if ( vlistGrid(vlistID1, i1) == vlistGrid(vlistID2, i2) ) break;

      //   printf("i1 %d i2 %d\n", i1, i2);

      if ( !ginit )
	{
	  gridIDs[i2] = genGrid(ngrids2, nfiles, ef, ginit, i1, nxblocks);
	  if ( gridIDs[i2] != -1 ) ginit = true;
	}
      else
	gridIDs[i2] = genGrid(ngrids2, nfiles, ef, ginit, i1, nxblocks);
    }


  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsize2 = 0;
  for ( int i = 0; i < ngrids2; ++i )
    {
      if ( gridIDs[i] != -1 ) 
	{
	  if ( gridsize2 == 0 ) gridsize2 = gridInqSize(gridIDs[i]);
	  if ( gridsize2 != gridInqSize(gridIDs[i]) ) cdoAbort("gridsize differ!");
	  vlistChangeGridIndex(vlistID2, i, gridIDs[i]);
	}
    }

  for ( varID = 0; varID < nvars2; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID2, varID);

      for ( int i = 0; i < ngrids2; ++i )
	{
	  if ( gridIDs[i] != -1 ) 
	    {
	      if ( gridID == vlistGrid(vlistID2, i) ) vars[varID] = true;
	      break;
	    }
	}
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);
	  
  double *array2 = (gridsize2 > 0) ? (double*) Malloc(gridsize2*sizeof(double)) : NULL;

  int tsID = 0;
  do
    {
      nrecs0 = pstreamInqTimestep(ef[0].streamID, tsID);
      for ( int fileID = 1; fileID < nfiles; fileID++ )
	{
	  int nrecs = pstreamInqTimestep(ef[fileID].streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records at time step %d of %s and %s differ!", tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      if ( nrecs0 > 0 ) pstreamDefTimestep(streamID2, tsID);
      
      for ( int recID = 0; recID < nrecs0; recID++ )
	{
	  pstreamInqRecord(ef[0].streamID, &varID, &levelID);
	  if ( cdoVerbose && tsID == 0 ) printf(" tsID, recID, varID, levelID %d %d %d %d\n", tsID, recID, varID, levelID);

	  for ( int fileID = 1; fileID < nfiles; fileID++ )
	    {
	      int varIDx, levelIDx;
	      pstreamInqRecord(ef[fileID].streamID, &varIDx, &levelIDx);
	    }

	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      int varID2   = vlistFindVar(vlistID2, varID);
	      int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
	      if ( cdoVerbose && tsID == 0 ) printf("varID %d %d levelID %d %d\n", varID, varID2, levelID, levelID2);

	      double missval = vlistInqVarMissval(vlistID2, varID2);
	      for ( int i = 0; i < gridsize2; i++ ) array2[i] = missval;

#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( int fileID = 0; fileID < nfiles; fileID++ )
		{
		  pstreamReadRecord(ef[fileID].streamID, ef[fileID].array, &ef[fileID].nmiss);

		  if ( vars[varID2] )
		    {
                      int gridsize = ef[fileID].gridsize;
		      for ( int i = 0; i < gridsize; ++i )
			array2[ef[fileID].gridindex[i]] = ef[fileID].array[i];
		    }
		}

	      pstreamDefRecord(streamID2, varID2, levelID2);

	      if ( vars[varID2] )
		{
		  int nmiss = 0;
		  for ( int i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

		  pstreamWriteRecord(streamID2, array2, nmiss);
		}
	      else
		pstreamWriteRecord(streamID2, ef[0].array, ef[0].nmiss);
	    }
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    pstreamClose(ef[fileID].streamID);

  pstreamClose(streamID2);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( ef[fileID].gridindex ) Free(ef[fileID].gridindex);
      if ( ef[fileID].array ) Free(ef[fileID].array);
    }

  if ( ef ) Free(ef);
  if ( array2 ) Free(array2);

  Free(gridIDs);
  if ( vars   ) Free(vars);
  if ( vars1  ) Free(vars1);

  cdoFinish();

  return 0;
}
