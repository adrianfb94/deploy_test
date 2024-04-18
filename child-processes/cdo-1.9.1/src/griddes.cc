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

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>

#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"
#include "griddes.h"
#include "error.h"


int grid_read(FILE *gfp, const char *dname);
int grid_read_pingo(FILE *gfp, const char *dname);


void gridInit(griddes_t *grid)
{
  grid->mask          = NULL;
  grid->xvals         = NULL;
  grid->yvals         = NULL;
  grid->xbounds       = NULL;
  grid->ybounds       = NULL;
  grid->area          = NULL;
  grid->type          = CDI_UNDEFID;
  grid->datatype      = CDI_UNDEFID;
  grid->size          = 0;
  grid->xsize         = 0;
  grid->ysize         = 0;
  grid->np            = 0;
  grid->lcomplex      = 1;
  grid->ntr           = 0;
  grid->nvertex       = 0;
  grid->genBounds     = false;
  grid->def_xfirst    = false;
  grid->def_yfirst    = false;
  grid->def_xlast     = false;
  grid->def_ylast     = false;
  grid->def_xinc      = false;
  grid->def_yinc      = false;
  grid->xfirst        = 0;
  grid->yfirst        = 0;
  grid->xlast         = 0;
  grid->ylast         = 0;
  grid->xinc          = 0;
  grid->yinc          = 0;
  grid->nd            = 0;
  grid->ni            = 0;
  grid->ni2           = 0;
  grid->ni3           = 0;
  grid->number        = 0;
  grid->position      = 0;
  grid->uuid[0]       = 0;
  grid->path[0]       = 0;
  grid->xname[0]      = 0;
  grid->xlongname[0]  = 0;
  grid->xunits[0]     = 0;
  grid->yname[0]      = 0;
  grid->ylongname[0]  = 0;
  grid->yunits[0]     = 0;
  grid->xdimname[0]   = 0;
  grid->ydimname[0]   = 0;
  grid->vdimname[0]   = 0;
  grid->uvRelativeToGrid = false;
  grid->scanningMode  = 64;
  /* scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
               64  = 128 * 0                + 64 *        1         + 32 * 0
               00  = 128 * 0                + 64 *        0         + 32 * 0
               96  = 128 * 0                + 64 *        1         + 32 * 1
     Default / implicit scanning mode is 64:
                        i and j scan positively, i points are consecutive (row-major)        */
}


int getoptname(char *optname, const char *optstring, int nopt)
{
  int nerr = 0;
  const char *pname = optstring;
  const char *pend  = optstring;

  for ( int i = 0; i < nopt; i++ )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	break;
      else
	pname = pend + 1;
    }

  if ( pend )
    {
      pend = strchr(pname, ',');
      size_t namelen;
      if ( pend == NULL )
	namelen = strlen(pname);
      else
	namelen = pend - pname;

      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return nerr;
}


int gridDefine(griddes_t grid)
{
  int gridID = CDI_UNDEFID;

  switch ( grid.type )
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_PROJECTION:
      {
	if ( grid.size != 1 )
	  {
	    if ( grid.xsize == 0 ) Error("xsize undefined!");
	    if ( grid.ysize == 0 ) Error("ysize undefined!");
	  }

	if ( grid.size == 0 ) grid.size = (long)grid.xsize*grid.ysize;

	if ( grid.size != (long)grid.xsize*grid.ysize )
	  Error("Inconsistent grid declaration: xsize*ysize!=gridsize (xsize=%d ysize=%d gridsize=%d)",
		grid.xsize, grid.ysize, grid.size);

	if ( grid.size < 0 || grid.size > INT_MAX ) Error("grid size (%ld) out of bounds (0 - %d)!", grid.size, INT_MAX);

	gridID = gridCreate(grid.type, grid.size);

	if ( grid.xsize > 0 ) gridDefXsize(gridID, grid.xsize);
	if ( grid.ysize > 0 ) gridDefYsize(gridID, grid.ysize);
	if ( grid.np    > 0 ) gridDefNP(gridID, grid.np);

        if ( grid.uvRelativeToGrid ) gridDefUvRelativeToGrid(gridID, 1);
	if ( grid.nvertex ) gridDefNvertex(gridID, grid.nvertex);

	if ( (grid.def_xfirst || grid.def_xlast || grid.def_xinc) && grid.xvals == NULL )
	  {
	    grid.xvals = (double*) Malloc(grid.xsize*sizeof(double));
	    gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, grid.xvals);

	    if ( grid.genBounds && grid.xbounds == NULL && grid.xsize > 1 )
	      {
		grid.nvertex = 2;
		grid.xbounds = (double*) Malloc(grid.xsize*grid.nvertex*sizeof(double));
		for ( int i = 0; i < (int) grid.xsize-1; i++ )
		  {
		    grid.xbounds[2*i+1]   = 0.5*(grid.xvals[i] + grid.xvals[i+1]);
		    grid.xbounds[2*(i+1)] = 0.5*(grid.xvals[i] + grid.xvals[i+1]);
		  }
		grid.xbounds[0] = 2*grid.xvals[0] - grid.xbounds[1];
		grid.xbounds[2*grid.xsize-1] = 2*grid.xvals[grid.xsize-1] - grid.xbounds[2*(grid.xsize-1)];
	      }
	  }

	if ( (grid.def_yfirst || grid.def_ylast || grid.def_yinc) && grid.yvals == NULL )
	  {
	    if ( ! grid.def_ylast ) grid.ylast = grid.yfirst;
	    grid.yvals = (double*) Malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, grid.yvals);

	    if ( grid.genBounds && grid.ybounds == NULL && grid.ysize > 1 )
	      {
		grid.nvertex = 2;
		grid.ybounds = (double*) Malloc(grid.ysize*grid.nvertex*sizeof(double));
		for ( int i = 0; i < (int) grid.ysize-1; i++ )
		  {
		    grid.ybounds[2*i+1]   = 0.5*(grid.yvals[i] + grid.yvals[i+1]);
		    grid.ybounds[2*(i+1)] = 0.5*(grid.yvals[i] + grid.yvals[i+1]);
		  }

		if ( grid.yvals[0] > grid.yvals[grid.ysize-1] )
		  {
		    grid.ybounds[0] = 90;
		    grid.ybounds[grid.ysize*grid.nvertex-1] = -90;
		  }
		else
		  {
		    grid.ybounds[0] = -90;
		    grid.ybounds[grid.ysize*grid.nvertex-1] = 90;
		  }
	      }
	  }

	if ( grid.xvals )   { gridDefXvals(gridID, grid.xvals); Free(grid.xvals); }
	if ( grid.yvals )   { gridDefYvals(gridID, grid.yvals); Free(grid.yvals); }
	if ( grid.xbounds ) { gridDefXbounds(gridID, grid.xbounds); Free(grid.xbounds); }
	if ( grid.ybounds ) { gridDefYbounds(gridID, grid.ybounds); Free(grid.ybounds); }
	if ( grid.mask )    { gridDefMask(gridID, grid.mask); Free(grid.mask); }

	break;
      }
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
	if ( grid.size == 0 )
          grid.size = (grid.type == GRID_CURVILINEAR) ? grid.xsize*grid.ysize : grid.xsize;

	gridID = gridCreate(grid.type, grid.size);

	if ( grid.type == GRID_CURVILINEAR )
	  {
	    if ( grid.xsize == 0 ) Error("xsize undefined!");
	    if ( grid.ysize == 0 ) Error("ysize undefined!");
	    gridDefXsize(gridID, grid.xsize);
	    gridDefYsize(gridID, grid.ysize);
	  }
	else
	  {
	    if ( grid.nvertex > 0 ) gridDefNvertex(gridID, grid.nvertex);
	    if ( grid.number > 0 )
	      {
		gridDefNumber(gridID, grid.number);
		if ( grid.position >= 0 ) gridDefPosition(gridID, grid.position);
	      }
	    if ( *grid.path ) gridDefReference(gridID, grid.path);
	  }

	if ( grid.xvals )   { gridDefXvals(gridID, grid.xvals); Free(grid.xvals); }
	if ( grid.yvals )   { gridDefYvals(gridID, grid.yvals); Free(grid.yvals); }
	if ( grid.area )    { gridDefArea(gridID, grid.area); Free(grid.area); }
	if ( grid.xbounds ) { gridDefXbounds(gridID, grid.xbounds); Free(grid.xbounds); }
	if ( grid.ybounds ) { gridDefYbounds(gridID, grid.ybounds); Free(grid.ybounds); }
	if ( grid.mask )    { gridDefMask(gridID, grid.mask); Free(grid.mask); }

	break;
      }
    case GRID_SPECTRAL:
      {
	if ( grid.ntr == 0 )
	  Error("truncation undefined!");
	if ( grid.size == 0 )
	  grid.size = (grid.ntr+1) * (grid.ntr+2);

	gridID = gridCreate(grid.type, grid.size);

	gridDefTrunc(gridID, grid.ntr);
	gridDefComplexPacking(gridID, grid.lcomplex);

	break;
      }
    case GRID_GME:
      {
	if ( grid.nd   == 0 ) Error("nd undefined!");
	if ( grid.ni   == 0 ) Error("ni undefined!");
	if ( grid.size == 0 ) Error("size undefined!");

	gridID = gridCreate(grid.type, grid.size);

	gridDefParamGME(gridID, grid.nd, grid.ni, grid.ni2, grid.ni3);
	
	if ( grid.mask ) { gridDefMask(gridID, grid.mask); Free(grid.mask); }

	break;
      }
    default:
      {
	if ( grid.type == -1 )
	  Error("Undefined grid type!");
	else
	  Error("Unsupported grid type: %s", gridNamePtr(grid.type));

	break;
      }
    }

  if ( grid.datatype != CDI_UNDEFID ) gridDefDatatype(gridID, grid.datatype);

  if ( grid.uuid[0] )      gridDefUUID(gridID, grid.uuid);

  if ( grid.xname[0]     ) cdiGridDefKeyStr(gridID, CDI_KEY_XNAME,     strlen(grid.xname)+1, grid.xname);
  if ( grid.xlongname[0] ) cdiGridDefKeyStr(gridID, CDI_KEY_XLONGNAME, strlen(grid.xlongname)+1, grid.xlongname);
  if ( grid.xunits[0]    ) cdiGridDefKeyStr(gridID, CDI_KEY_XUNITS,    strlen(grid.xunits)+1, grid.xunits);
  if ( grid.yname[0]     ) cdiGridDefKeyStr(gridID, CDI_KEY_YNAME,     strlen(grid.yname)+1, grid.yname);
  if ( grid.ylongname[0] ) cdiGridDefKeyStr(gridID, CDI_KEY_YLONGNAME, strlen(grid.ylongname)+1, grid.ylongname);
  if ( grid.yunits[0]    ) cdiGridDefKeyStr(gridID, CDI_KEY_YUNITS,    strlen(grid.yunits)+1, grid.yunits);
  if ( grid.xdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_XDIMNAME,  strlen(grid.xdimname)+1, grid.xdimname);
  if ( grid.ydimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_YDIMNAME,  strlen(grid.ydimname)+1, grid.ydimname);
  if ( grid.vdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_VDIMNAME,  strlen(grid.vdimname)+1, grid.vdimname);

  return gridID;
}


int cdoDefineGrid(const char *gridfile)
{
  int gridID = -1;
  size_t len;
  bool isreg = false;
  bool lalloc = false;

  char *filename = expand_filename(gridfile);
  if ( filename )
    lalloc = true;
  else
    filename = (char *) gridfile;

  int fileno = open(filename, O_RDONLY);
  if ( fileno >= 0 )
    {
      struct stat filestat;
      if ( fstat(fileno, &filestat) == 0 )
	isreg = S_ISREG(filestat.st_mode);
    }

  if ( fileno == -1 || !isreg )
    {
      if ( isreg ) close(fileno);

      gridID = grid_from_name(gridfile);

      if ( gridID == -1 ) cdoAbort("Open failed on %s!", gridfile);
    }
  else
    {
      char buffer[4];
      if ( read(fileno, buffer, 4) != 4 )
	SysError("Read grid from %s failed!", filename);

      close(fileno);

      if ( cmpstrlen(buffer, "CDF", len) == 0 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from NetCDF file");
	  gridID = gridFromNCfile(filename);
	}

      if ( gridID == -1 )
	{
	  if ( cmpstrlen(buffer+1, "HDF", len) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from HDF5 file");
	      gridID = gridFromH5file(filename);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cmpstrlen(buffer+1, "HDF", len) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from NetCDF4 file");
	      gridID = gridFromNCfile(filename);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from CDI file");
	  openLock();
	  int streamID = streamOpenRead(filename);
	  openUnlock();
	  if ( streamID >= 0 )
	    {
	      int vlistID = streamInqVlist(streamID);
	      gridID  = vlistGrid(vlistID, 0);
	      streamClose(streamID);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from ASCII file");
	  FILE *gfp = fopen(filename, "r");
	  //size_t buffersize = 20*1024*1024;
	  //char *buffer = (char*) Malloc(buffersize);
	  //setvbuf(gfp, buffer, _IOFBF, buffersize);
	  gridID = grid_read(gfp, filename);
	  fclose(gfp);
	  //free(buffer);
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from PINGO file");
	  FILE *gfp = fopen(filename, "r");
	  gridID = grid_read_pingo(gfp, filename);
	  fclose(gfp);
	}

      if ( gridID == -1 ) cdoAbort("Invalid grid description file %s!", filename);
    }

  if ( lalloc ) Free(filename);

  return gridID;
}


void cdo_set_grids(const char *gridarg)
{
  char gridfile[4096];
  int nfile = 0;

  while ( getoptname(gridfile, gridarg, nfile++) == 0 )
    {      
      (void) cdoDefineGrid(gridfile);
    }
}
