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

      Mrotuvb     mrotuvb          Backward rotation for MPIOM data
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"

/*
!----------------------------------------------------------------------
!
!     rotation of vectors: in ocean models with rotated grids velocity
!     vectors are given in the direction of grid lines and rows. they
!     have to be rotated in latitudinal and longitudinal direction.
!
!     note: this routine assumes positive meridional flow for a flow
!           from grid point(i,j) to grid point(i,j+1) and positive
!           zonal flow for a flow from grid point(i,j) to point(i+1,j).
!           this is not the case for mpi-om!
!
!           if this routine is used to rotate data of mpi-om, the
!           logical change_sign_v needs to be true.
!j. jungclaus: 22.01.04:
!note here for the coupling fields u-i,v_j are on the non-verlapping
! (ie-2=ix) grid, furthermore, the velocity fields were previously
! interpolated onto the scalar points !
!
!h.haak: 07.10.2005 vectorisation and omp directives
!malte: use outside mpiom 02.06.2006     
!----------------------------------------------------------------------
*/

void rotate_uv2(double *u_i, double *v_j, int ix, int iy,
		double *lon, double *lat, double *u_lon, double *v_lat)
{
  /*
 real,intent(in)      :: u_i(ix,iy,iz),v_j(ix,iy,iz)                  ! vector components in i-j-direction
 real,intent(out)     :: u_lon(ix,iy,iz),v_lat(ix,iy,iz)              ! vector components in lon-lat direction
 real,intent(in)      :: lat(ix,iy),lon(ix,iy)                        ! latitudes and longitudes
  */
  double dlat_i, dlat_j,dlon_i,dlon_j,dist_i,dist_j;
  double lat_factor;
  double absold, absnew;  /* velocity vector lengths */
  int  i, j, ip1, im1, jp1, jm1;
  double pi = 3.14159265359;


  /* specification whether change in sign is needed for the input arrays */
  bool change_sign_u = false;
  bool change_sign_v = true;

  /* initialization */
  for ( i = 0; i < ix*iy; i++ )
    {
      v_lat[i] = 0;
      u_lon[i] = 0;
    }

  /* change sign */
  if ( change_sign_u )
    for ( i = 0; i < ix*iy; i++ ) u_i[i] *= -1;

  if ( change_sign_v )
    for ( i = 0; i < ix*iy; i++ ) v_j[i] *= -1;

  /* rotation */
  for ( j = 0; j < iy; j++ )
    for ( i = 0; i < ix; i++ )
      {
	ip1 = i + 1;
	im1 = i - 1;
	jp1 = j + 1;
	jm1 = j - 1;
	if ( ip1 >= ix ) ip1 = 0;   /* the 0-meridian */
	if ( im1 <   0 ) im1 = ix-1;
	if ( jp1 >= iy ) jp1 = j;   /* treatment of the last..  */
	if ( jm1 <   0 ) jm1 = j;   /* .. and the fist grid-row */

	/* difference in latitudes */
	dlat_i = lat[IX2D(j,ip1,ix)] - lat[IX2D(j,im1,ix)];
	dlat_j = lat[IX2D(jp1,i,ix)] - lat[IX2D(jm1,i,ix)];

	/* difference in longitudes */
	dlon_i = lon[IX2D(j,ip1,ix)] - lon[IX2D(j,im1,ix)];
	if ( dlon_i >   pi  ) dlon_i -= 2*pi;
	if ( dlon_i < (-pi) ) dlon_i += 2*pi;
	dlon_j = lon[IX2D(jp1,i,ix)] - lon[IX2D(jm1,i,ix)];
	if ( dlon_j >   pi  ) dlon_j -= 2*pi;
	if ( dlon_j < (-pi) ) dlon_j += 2*pi;

	lat_factor = cos(lat[IX2D(j,i,ix)]);
	dlon_i = dlon_i * lat_factor;
	dlon_j = dlon_j * lat_factor;

	/* projection by scalar product */
	u_lon[IX2D(j,i,ix)] = u_i[IX2D(j,i,ix)]*dlon_i + v_j[IX2D(j,i,ix)]*dlat_i;
	v_lat[IX2D(j,i,ix)] = u_i[IX2D(j,i,ix)]*dlon_j + v_j[IX2D(j,i,ix)]*dlat_j;

	dist_i = sqrt(dlon_i*dlon_i + dlat_i*dlat_i);
	dist_j = sqrt(dlon_j*dlon_j + dlat_j*dlat_j);

	if ( fabs(dist_i) > 0 && fabs(dist_j) > 0 )
	  {
	    u_lon[IX2D(j,i,ix)] /= dist_i;
	    v_lat[IX2D(j,i,ix)] /= dist_j;
	  }
	else
	  {
	    u_lon[IX2D(j,i,ix)] = 0.0;
	    v_lat[IX2D(j,i,ix)] = 0.0;
	  }

	if ( cdoVerbose )
	  {
	    absold = sqrt(u_i[IX2D(j,i,ix)]*u_i[IX2D(j,i,ix)] + v_j[IX2D(j,i,ix)]*v_j[IX2D(j,i,ix)]);
	    absnew = sqrt(u_lon[IX2D(j,i,ix)]*u_lon[IX2D(j,i,ix)] + v_lat[IX2D(j,i,ix)]*v_lat[IX2D(j,i,ix)]);

	    if ( i%20 == 0 && j%20 == 0 && absold > 0 )
	      {
		printf("(absold,absnew) %d %d %g %g %g %g %g %g\n",
		       j+1, i+1, absold, absnew, u_i[IX2D(j,i,ix)], v_j[IX2D(j,i,ix)], u_lon[IX2D(j,i,ix)], v_lat[IX2D(j,i,ix)]);
		
		/* test orthogonality */
		if ( (dlon_i*dlon_j + dlat_j*dlat_i) > 0.1 )            
		  fprintf(stderr, "orthogonal? %d %d %g\n", j+1, i+1, (dlon_i*dlon_j + dlat_j*dlat_i));
	      }
	  }
      }
}


static
void uv_to_p_grid(int nlon, int nlat, double *grid1x, double *grid1y, 
		  double *grid2x, double *grid2y, double *grid3x, double *grid3y)
{
  double gx, gy;
  double gx2, gy2;

  int gridsizex = (nlon+2)*nlat;
  double *gxhelp = (double*) Malloc(gridsizex*sizeof(double));
  double *gyhelp = (double*) Malloc(gridsizex*sizeof(double));

  /* load to a help field */
  for ( int j = 0; j < nlat; j++ )
    for ( int i = 0; i < nlon; i++ )
      {
	gxhelp[IX2D(j,i+1,nlon+2)] = grid1x[IX2D(j,i,nlon)];
	gyhelp[IX2D(j,i+1,nlon+2)] = grid1y[IX2D(j,i,nlon)];
      }

  /* make help field cyclic */
  for ( int j = 0; j < nlat; j++ )
    {
      gxhelp[IX2D(j,0,nlon+2)]      = gxhelp[IX2D(j,nlon,nlon+2)];
      gxhelp[IX2D(j,nlon+1,nlon+2)] = gxhelp[IX2D(j,1,nlon+2)];
      gyhelp[IX2D(j,0,nlon+2)]      = gyhelp[IX2D(j,nlon,nlon+2)];
      gyhelp[IX2D(j,nlon+1,nlon+2)] = gyhelp[IX2D(j,1,nlon+2)];
    }

  /* interpolate u to scalar points */
  for ( int j = 0; j < nlat; j++ )
    for ( int i = 0; i < nlon; i++ )
      {
	grid3x[IX2D(j,i,nlon)] = (gxhelp[IX2D(j,i,nlon+2)]+gxhelp[IX2D(j,i+1,nlon+2)])*0.5;
	if ( (gxhelp[IX2D(j,i,nlon+2)] > 340 && gxhelp[IX2D(j,i+1,nlon+2)] <  20) ||
             (gxhelp[IX2D(j,i,nlon+2)] < 20  && gxhelp[IX2D(j,i+1,nlon+2)] > 340) )
	  {
	    if ( grid3x[IX2D(j,i,nlon)] < 180 )
	      grid3x[IX2D(j,i,nlon)] += 180;
	    else
	      grid3x[IX2D(j,i,nlon)] -= 180;
	  }
	    
	grid3y[IX2D(j,i,nlon)] = (gyhelp[IX2D(j,i,nlon+2)]+gyhelp[IX2D(j,i+1,nlon+2)])*0.5;
      }

  /* load to a help field */
  for ( int j = 0; j < nlat; j++ )
    for ( int i = 0; i < nlon; i++ )
      {
	gxhelp[IX2D(j,i+1,nlon+2)] = grid2x[IX2D(j,i,nlon)];
	gyhelp[IX2D(j,i+1,nlon+2)] = grid2y[IX2D(j,i,nlon)];
      }

  /* make help field cyclic */
  for ( int j = 0; j < nlat; j++ )
    {
      gxhelp[IX2D(j,0,nlon+2)]      = gxhelp[IX2D(j,nlon,nlon+2)];
      gxhelp[IX2D(j,nlon+1,nlon+2)] = gxhelp[IX2D(j,1,nlon+2)];
      gyhelp[IX2D(j,0,nlon+2)]      = gyhelp[IX2D(j,nlon,nlon+2)];
      gyhelp[IX2D(j,nlon+1,nlon+2)] = gyhelp[IX2D(j,1,nlon+2)];
    }

  /* interpolate v to scalar points */
  for ( int j = 1; j < nlat-1; j++ )
    for ( int i = 0; i < nlon; i++ )
      {
	gx = (gxhelp[IX2D(j,i+1,nlon+2)]+gxhelp[IX2D(j-1,i+1,nlon+2)])*0.5;
	if ( (gxhelp[IX2D(j,i+1,nlon+2)] > 340 && gxhelp[IX2D(j-1,i+1,nlon+2)] <  20) ||
             (gxhelp[IX2D(j,i+1,nlon+2)] < 20  && gxhelp[IX2D(j-1,i+1,nlon+2)] > 340) )
	  {
	    if ( gx < 180 )
	      gx += 180;
	    else
	      gx -= 180;
	  }
	    
	gy = (gyhelp[IX2D(j,i+1,nlon+2)]+gyhelp[IX2D(j-1,i+1,nlon+2)])*0.5;
       
	/* printf("%d %d %g %g %g %g \n", j, i, gx, gy, grid3x[IX2D(j,i,nlon)], grid3y[IX2D(j,i,nlon)]); */

	gx2 = (gx+grid3x[IX2D(j,i,nlon)])*0.5;
	if ( (gx > 340 && grid3x[IX2D(j,i,nlon)] <  20) ||
             (gx < 20  && grid3x[IX2D(j,i,nlon)] > 340) )
	  {
	    if ( gx2 < 180 )
	      gx2 += 180;
	    else
	      gx2 -= 180;
	  }
	    
	gy2 = (gy+grid3y[IX2D(j,i,nlon)])*0.5;

	grid3x[IX2D(j,i,nlon)] = gx2;
	grid3y[IX2D(j,i,nlon)] = gy2;

	/* printf("%d %d %g %g %g %g \n", j, i, gx2, gy2, grid3x[IX2D(j,i,nlon)], grid3y[IX2D(j,i,nlon)]); */
      }

  Free(gxhelp);
  Free(gyhelp);
}


void *Mrotuvb(void *argument)
{
  int nrecs, nrecs2;
  int levelID;
  int varID1, varID2;
  int nmiss1, nmiss2;
  bool gpint = true;

  cdoInitialize(argument);

  if ( operatorArgc() == 1 )
    if ( strcmp(operatorArgv()[0], "noint") == 0 ) gpint = false;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);

  int nvars = vlistNvars(vlistID1);
  if ( nvars > 1 ) cdoAbort("More than one variable found in %s",  cdoStreamName(0)->args);
  nvars = vlistNvars(vlistID2);
  if ( nvars > 1 ) cdoAbort("More than one variable found in %s",  cdoStreamName(1)->args);

  int gridID1 = vlistGrid(vlistID1, 0);
  int gridID2 = vlistGrid(vlistID2, 0);
  int gridsize = gridInqSize(gridID1);
  if ( gpint == TRUE  && gridID1 == gridID2 ) cdoAbort("Input grids are the same, use parameter >noint< to disable interpolation!");
  if ( gpint == FALSE && gridID1 != gridID2 ) cdoAbort("Input grids are not the same!");
  if ( gridsize != gridInqSize(gridID2) ) cdoAbort("Grids have different size!");

  if ( gridInqType(gridID1) != GRID_LONLAT      &&
       gridInqType(gridID1) != GRID_GAUSSIAN    &&
       gridInqType(gridID1) != GRID_CURVILINEAR )
    cdoAbort("Grid %s unsupported!", gridNamePtr(gridInqType(gridID1)));
  
  if ( gridInqType(gridID1) != GRID_CURVILINEAR )
    gridID1 = gridToCurvilinear(gridID1, 1);

  if ( gridsize != gridInqSize(gridID1) ) cdoAbort("Internal problem: gridsize changed!");

  if ( gridInqType(gridID2) != GRID_CURVILINEAR )
    gridID2 = gridToCurvilinear(gridID2, 1);

  if ( gridsize != gridInqSize(gridID2) ) cdoAbort("Internal problem: gridsize changed!");

  int nlon    = gridInqXsize(gridID1);
  int nlat    = gridInqYsize(gridID1);

  double *grid1x = (double*) Malloc(gridsize*sizeof(double));
  double *grid1y = (double*) Malloc(gridsize*sizeof(double));
  double *grid2x = (double*) Malloc(gridsize*sizeof(double));
  double *grid2y = (double*) Malloc(gridsize*sizeof(double));
  double *grid3x = (double*) Malloc(gridsize*sizeof(double));
  double *grid3y = (double*) Malloc(gridsize*sizeof(double));

  gridInqXvals(gridID1, grid1x);
  gridInqYvals(gridID1, grid1y);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID1, units);
    grid_to_degree(units, gridsize, grid1x, "grid1 center lon");
    gridInqYunits(gridID1, units);
    grid_to_degree(units, gridsize, grid1y, "grid1 center lat");
  }

  gridInqXvals(gridID2, grid2x);
  gridInqYvals(gridID2, grid2y);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID2, units);
    grid_to_degree(units, gridsize, grid2x, "grid2 center lon");
    gridInqYunits(gridID2, units);
    grid_to_degree(units, gridsize, grid2y, "grid2 center lat");
  }

  if ( gpint )
    {
      uv_to_p_grid(nlon, nlat, grid1x, grid1y, grid2x, grid2y, grid3x, grid3y);
    }
  else
    {
      memcpy(grid3x, grid1x, gridsize*sizeof(double));
      memcpy(grid3y, grid1y, gridsize*sizeof(double));
    }

  if ( grid1x ) Free(grid1x);
  if ( grid1y ) Free(grid1y);
  if ( grid2x ) Free(grid2x);
  if ( grid2y ) Free(grid2y);

  int gridID3 = gridCreate(GRID_CURVILINEAR, gridsize);
  gridDefDatatype(gridID3, gridInqDatatype(gridID1));
  gridDefXsize(gridID3, nlon);
  gridDefYsize(gridID3, nlat);
  gridDefXvals(gridID3, grid3x);
  gridDefYvals(gridID3, grid3y);

  for ( int i = 0; i < gridsize; i++ )
    {
      grid3x[i] *= DEG2RAD;
      grid3y[i] *= DEG2RAD;
    }

  int vlistID3 = vlistCreate();
  vlistCopy(vlistID3, vlistID1);
  vlistCat(vlistID3, vlistID2);

  int code1 = vlistInqVarCode(vlistID1, 0);
  int code2 = vlistInqVarCode(vlistID2, 0);

  if ( code1 == code2 ) vlistDefVarCode(vlistID3, 1, code1+1);
  
  vlistChangeGrid(vlistID3, gridID1, gridID3);
  vlistChangeGrid(vlistID3, gridID2, gridID3);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  if ( cdoVerbose ) vlistPrint(vlistID3);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  double missval1 = vlistInqVarMissval(vlistID1, 0);
  double missval2 = vlistInqVarMissval(vlistID2, 0);

  double *ufield  = (double*) Malloc(gridsize*sizeof(double));
  double *vfield  = (double*) Malloc(gridsize*sizeof(double));
  double *urfield = (double*) Malloc(gridsize*sizeof(double));
  double *vrfield = (double*) Malloc(gridsize*sizeof(double));

  double *uhelp = NULL, *vhelp = NULL;
  if ( gpint )
    {
      int gridsizex = (nlon+2)*nlat;
      uhelp   = (double*) Malloc(gridsizex*sizeof(double));
      vhelp   = (double*) Malloc(gridsizex*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID3, tsID);

      nrecs2 = pstreamInqTimestep(streamID2, tsID);

      if ( nrecs != nrecs2 ) cdoAbort("Input streams have different number of levels!");
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID1, &levelID);
	  pstreamInqRecord(streamID2, &varID2, &levelID);

	  pstreamReadRecord(streamID1, ufield, &nmiss1);
	  pstreamReadRecord(streamID2, vfield, &nmiss2);

	  /* remove missing values */
	  if ( nmiss1 || nmiss2 )
	    {
	      for ( int i = 0; i < gridsize; i++ )
		{
		  if ( DBL_IS_EQUAL(ufield[i], missval1) ) ufield[i] = 0;
		  if ( DBL_IS_EQUAL(vfield[i], missval2) ) vfield[i] = 0;
		}
	    }

	  if ( gpint )
	    {
	      /* load to a help field */
	      for ( int j = 0; j < nlat; j++ )
		for ( int i = 0; i < nlon; i++ )
		  {
		    uhelp[IX2D(j,i+1,nlon+2)] = ufield[IX2D(j,i,nlon)];
		    vhelp[IX2D(j,i+1,nlon+2)] = vfield[IX2D(j,i,nlon)];
		  }

	      /* make help field cyclic */
	      for ( int j = 0; j < nlat; j++ )
		{
		  uhelp[IX2D(j,0,nlon+2)]      = uhelp[IX2D(j,nlon,nlon+2)];
		  uhelp[IX2D(j,nlon+1,nlon+2)] = uhelp[IX2D(j,1,nlon+2)];
		  vhelp[IX2D(j,0,nlon+2)]      = vhelp[IX2D(j,nlon,nlon+2)];
		  vhelp[IX2D(j,nlon+1,nlon+2)] = vhelp[IX2D(j,1,nlon+2)];
		}

	      /* interpolate on pressure points */
	      for ( int j = 1; j < nlat; j++ )
		for ( int i = 0; i < nlon; i++ )
		  {
		    ufield[IX2D(j,i,nlon)] = (uhelp[IX2D(j,i,nlon+2)]+uhelp[IX2D(j,i+1,nlon+2)])*0.5;
		    vfield[IX2D(j,i,nlon)] = (vhelp[IX2D(j-1,i+1,nlon+2)]+vhelp[IX2D(j,i+1,nlon+2)])*0.5;
		  }
	    }

	  for ( int i = 0; i < nlon; i++ )
	    {
	      ufield[IX2D(0,i,nlon)] = 0;
	      vfield[IX2D(0,i,nlon)] = 0;
	    }

	  /* rotate*/
	  
	  rotate_uv2(ufield, vfield, nlon, nlat, grid3x, grid3y, urfield, vrfield);

	  /* calc lat, lon, Auv and alpha */
	  /*
          {
	  double lat, lon, auv, alpha;
	  for ( int j = 1; j < nlat-1; j += 3 )
	    for ( int i = 0; i < nlon; i += 3 )
	      {
		lat = grid3y[IX2D(j,i,nlon)]*RAD2DEG;
		lon = grid3x[IX2D(j,i,nlon)]*RAD2DEG; 
		auv = sqrt(urfield[IX2D(j,i,nlon)]*urfield[IX2D(j,i,nlon)] +
			   vrfield[IX2D(j,i,nlon)]*vrfield[IX2D(j,i,nlon)]);
		alpha = atan2(vrfield[IX2D(j,i,nlon)], urfield[IX2D(j,i,nlon)]);
		alpha = 90. - alpha*RAD2DEG;

		if ( alpha <   0 ) alpha += 360.;
		if ( alpha > 360 ) alpha -= 360.;

		printf("%g %g %g %g\n", lon, lat, alpha, auv);
	      }
          }
	  */
	  nmiss1 = 0;
	  nmiss2 = 0;
	  pstreamDefRecord(streamID3, 0, levelID);
	  pstreamWriteRecord(streamID3, urfield, nmiss1);     
	  pstreamDefRecord(streamID3, 1, levelID);
	  pstreamWriteRecord(streamID3, vrfield, nmiss2);     
	}

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( ufield  ) Free(ufield);
  if ( vfield  ) Free(vfield);
  if ( urfield ) Free(urfield);
  if ( vrfield ) Free(vrfield);
  if ( gpint )
    {
      if ( uhelp   ) Free(uhelp);
      if ( vhelp   ) Free(vhelp);
    }
  if ( grid3x  ) Free(grid3x);
  if ( grid3y  ) Free(grid3y);

  cdoFinish();

  return 0;
}
