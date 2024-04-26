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
#include "cdo_int.h"
#include "grid.h"
#include "griddes.h"


static
void skip_nondigit_lines(FILE *gfp)
{
  int c;

  if ( feof(gfp) ) return;

  while (1)
    {
      do
	c = fgetc(gfp);
      while ( (isspace(c) || c == ',') && c != EOF );

      if ( c == EOF || isdigit (c) || c == '.' || c == '+' || c == '-' ) break;
      else
	while ( c != '\n' && c != EOF )
	  c = fgetc(gfp);
    }

  ungetc(c, gfp);
}

static
int input_ival(FILE *gfp, int *ival)
{
  skip_nondigit_lines(gfp);

  if ( feof(gfp) ) return 0;

  *ival = 0;
  int read_items = fscanf(gfp, "%d", ival);

  return read_items;
}


int input_darray(FILE *gfp, int n_values, double *array)
{
  if ( n_values <= 0 ) return 0;

  int read_items = 0;
  for ( int i = 0; i < n_values; i++ )
    {
      skip_nondigit_lines(gfp);

      if ( feof(gfp) ) break;

      read_items += fscanf(gfp, "%lg", &array[i]);

      if ( feof(gfp) ) break;
    }

  return read_items;
}


int grid_read_pingo(FILE *gfp, const char *dname)
{
  UNUSED(dname);
  int gridID = -1;
  int i;

  griddes_t grid;
  gridInit(&grid);

  int nlon, nlat;
  if ( ! input_ival(gfp, &nlon) ) return gridID;
  if ( ! input_ival(gfp, &nlat) ) return gridID;

  if ( nlon > 0 && nlon < 9999 && nlat > 0 && nlat < 9999 )
    {
      grid.xsize = nlon;
      grid.ysize = nlat;

      grid.xvals = (double*) Malloc(grid.xsize*sizeof(double));
      grid.yvals = (double*) Malloc(grid.ysize*sizeof(double));

      if ( ! input_ival(gfp, &nlon) ) return gridID;
      if ( nlon == 2 )
	{
	  if ( input_darray(gfp, 2, grid.xvals) != 2 ) return gridID;
	  grid.xvals[1] -= 360 * floor((grid.xvals[1] - grid.xvals[0]) / 360);

	  if ( grid.xsize > 1 )
	    if ( IS_EQUAL(grid.xvals[0], grid.xvals[1]) )
	      grid.xvals[1] += 360;

	  for ( i = 0; i < (int)grid.xsize; i++ )
	    grid.xvals[i] = grid.xvals[0] + i*(grid.xvals[1] - grid.xvals[0]);
	}
      else if ( nlon == (int)grid.xsize )
	{
	  if ( input_darray(gfp, nlon, grid.xvals) != nlon ) return gridID;
	  for ( i = 0; i < nlon - 1; i++ )
	    if ( grid.xvals[i+1] <= grid.xvals[i] ) break;

	  for ( i++; i < nlon; i++ )
	    {
	      grid.xvals[i] += 360;
	      if ( i < nlon - 1 && grid.xvals[i+1] + 360 <= grid.xvals[i] )
		{
		  cdoPrint("Longitudes are not in ascending order!");
		  return gridID;
		}
	    }
	}
      else
	return gridID;

      if ( ! input_ival(gfp, &nlat) ) return gridID;
      if ( nlat == 2 )
	{
	  if ( input_darray(gfp, 2, grid.yvals) != 2 ) return gridID;
	  for ( i = 0; i < (int)grid.ysize; i++ )
	    grid.yvals[i] = grid.yvals[0] + i*(grid.yvals[1] - grid.yvals[0]);
	}
      else if ( nlat == (int)grid.ysize )
	{
	  if ( input_darray(gfp, nlat, grid.yvals) != nlat ) return gridID;
	}
      else
	return gridID;

      if ( grid.yvals[0]      >  90.001  || 
	   grid.yvals[nlat-1] >  90.001  || 
	   grid.yvals[0]      < -90.001  || 
	   grid.yvals[nlat-1] < -90.001 )
	{
	  cdoPrint("Latitudes must be between 90 and -90!");
	  return gridID;
	}

      for ( i = 0; i < nlat - 1; i++ )
	if ( IS_EQUAL(grid.yvals[i+1], grid.yvals[i]) || (i < nlat - 2 &&
	    ((grid.yvals[i+1] > grid.yvals[i]) != (grid.yvals[i+2] > grid.yvals[i+1]))) )
	  {
	    cdoPrint("Latitudes must be in descending or ascending order!");
	    return gridID;
	  }
		    
      bool lgauss = false;
      if ( nlat > 2 ) /* check if gaussian */
	{
	  double *yvals, *yw;
	  yvals = (double*) Malloc(grid.ysize*sizeof(double));
	  yw    = (double*) Malloc(grid.ysize*sizeof(double));
	  gaussaw(yvals, yw, grid.ysize);
	  Free(yw);
	  for ( i = 0; i < (int) grid.ysize; i++ )
	    yvals[i] = asin(yvals[i])*RAD2DEG;

	  for ( i = 0; i < (int) grid.ysize; i++ )
	    if ( fabs(yvals[i] - grid.yvals[i]) > ((yvals[0] - yvals[1])/500) ) break;
		      
	  if ( i == (int) grid.ysize ) lgauss = true;

	  Free(yvals);
	}

      if ( lgauss )
	grid.type = GRID_GAUSSIAN;
      else
	grid.type = GRID_LONLAT;
    }
  
  if ( grid.type != CDI_UNDEFID ) gridID = gridDefine(grid);

  return gridID;
}
