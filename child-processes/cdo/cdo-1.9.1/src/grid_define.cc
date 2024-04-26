#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"

/*
@Function  cdo_define_destagered_grid
@Title     Define a de-staggered grid for U and V

@Prototype int cdo_define_destagered_grid(int gridID_u_stag, int gridID_v_stag, double *destagGridOffsets)
@Parameter
    @Item  grid_u_stag       Staggered grid of u-wind component
    @Item  grid_v_stag       Staggered grid of v-wind component
    @Item  grid_uv_destag    Destaggered grid of uv-wind

@Description
The function @func{cdo_define_destagered_grid} defines a de-staggered grid for U and V

@EndFunction
*/
int cdo_define_destagered_grid(int gridID_u_stag, int gridID_v_stag, double *destagGridOffsets)
{
/* Example of horizontal grids (Hirlam LAMH_D11):
     U : lonlat       > size      : dim = 399300  nlon = 726  nlat = 550
                        rlon      : first = -30.15  last = 42.35  inc = 0.1  degrees
                        rlat      : first = -30.8  last = 24.1  inc = 0.1  degrees
                        northpole : lon = -195  lat = 30
     V : lonlat       > size      : dim = 399300  nlon = 726  nlat = 550
                        rlon      : first = -30.2  last = 42.3  inc = 0.1  degrees
                        rlat      : first = -30.75  last = 24.15  inc = 0.1  degrees
                        northpole : lon = -195  lat = 30
=>   RESULT:
     R : lonlat       > size      : dim = 399300  nlon = 726  nlat = 550
                        rlon      : first = -30.2  last = 42.3  inc = 0.1  degrees
                        rlat      : first = -30.8  last = 24.1  inc = 0.1  degrees
                        northpole : lon = -195  lat = 30
*/
  if ( cdoDebugExt )
    cdoPrint("%s(gridID_u=%d,gridID_v=%d,destagGridOffsets(%02.1f,%02.1f)) ...\n",
             __func__, gridID_u_stag, gridID_v_stag, destagGridOffsets[0],destagGridOffsets[1]);

  if ( cdoDebugExt > 1 )
    {
      cdo_print_grid(gridID_u_stag, 1);
      cdo_print_grid(gridID_v_stag, 1);
    }

  int gridtype = gridInqType(gridID_u_stag);
  int xsize = gridInqXsize(gridID_u_stag);
  int ysize = gridInqYsize(gridID_u_stag);

  double xfirst_U = gridInqXval(gridID_u_stag, 0);   // staggered grid of u-wind
  double yfirst_U = gridInqYval(gridID_u_stag, 0);
  double xlast_U  = gridInqXval(gridID_u_stag, xsize-1);
  double ylast_U  = gridInqYval(gridID_u_stag, ysize-1);
  double xfirst_V = gridInqXval(gridID_v_stag, 0);   // staggered grid of v-wind
  double yfirst_V = gridInqYval(gridID_v_stag, 0);
  double xlast_V  = gridInqXval(gridID_v_stag, xsize-1);
  double ylast_V  = gridInqYval(gridID_v_stag, ysize-1);
  double xinc     = gridInqXinc(gridID_u_stag);
  double yinc     = gridInqYinc(gridID_u_stag);

  int gridID_uv_destag = gridDuplicate(gridID_u_stag);

  if ( cdoDebugExt )
    {
      cdo_print_grid(gridID_uv_destag, 1);

      cdoPrint("%s(): (gridXsize=%d, gridYsize=%d)", __func__, xsize, ysize);
      cdoPrint("%s(): (xfirst_U = %3.2f; yfirst_U = %3.2f); (xfirst_V = %3.2f; yfirst_V = %3.2f)", __func__, xfirst_U, yfirst_U, xfirst_V, yfirst_V);
      cdoPrint("%s(): (xlast_U  = %3.2f; ylast_U  = %3.2f); (xlast_V  = %3.2f; ylast_V  = %3.2f)", __func__, xlast_U, ylast_U, xlast_V, ylast_V);
    }

  double xfirst = 0, xlast = 0, yfirst = 0, ylast = 0;
  if ( IS_EQUAL(destagGridOffsets[0], -0.5) && IS_EQUAL(destagGridOffsets[1], -0.5) )
    {
      xfirst =  xfirst_V;
      xlast  =  xlast_V;
      
      yfirst =  yfirst_U;
      ylast  =  ylast_U;
    }
  else if ( IS_EQUAL(destagGridOffsets[0], 0.5) && IS_EQUAL(destagGridOffsets[1], 0.5) )
    {
      xfirst =  xfirst_V + xinc*destagGridOffsets[0];
      xlast  =  xlast_V  + xinc*destagGridOffsets[0];

      yfirst =  yfirst_U + yinc*destagGridOffsets[1];
      ylast  =  ylast_U  + yinc*destagGridOffsets[1];
    }
  else
    cdoAbort("%s() Unsupported destaggered grid offsets! We support only: (-0.5,-0.5) or (0.5,0.5)", __func__);

  double *xvals = (double *) Malloc(xsize*sizeof(double));
  gridGenXvals(xsize, xfirst, xlast, xinc, xvals);
  gridDefXvals(gridID_uv_destag, xvals);
  Free(xvals);
  
  double *yvals = (double *) Malloc(ysize*sizeof(double));
  gridGenYvals(gridtype, ysize, yfirst, ylast, yinc, yvals);
  gridDefYvals(gridID_uv_destag, yvals);
  Free(yvals);
  
  if ( cdoDebugExt )
    {
      cdoPrint("%s():", __func__);
      cdo_print_grid(gridID_uv_destag, 1);
    }

  return gridID_uv_destag;
}

/*
@Function  cdo_define_sample_grid
@Title     Define a sampled grid of another grid

@Prototype int cdo_define_sample_grid(int gridSrcID, int sampleFactor)
@Parameter
    @Item  gridSrcID       Source grid
    @Item  sampleFactor    sampleFactor; typically 2,3,4 ...

@Description
The function @func{cdo_define_sample_grid} defines a sampled grid of another grid

@EndFunction
*/
int cdo_define_sample_grid(int gridSrcID, int sampleFactor)
{
/* Example of horizontal grids (Harmonie HARM36_L25):
            #
            # gridID 2
            #
            gridtype  = projection
            gridsize  = 622521
            xsize     = 789
            ysize     = 789
            xunits    = "m"
            yunits    = "m"
            xfirst    = 0
            xinc      = 2500
            yfirst    = 0
            yinc      = 2500
            grid_mapping = Lambert_Conformal
            grid_mapping_name = lambert_conformal_conic
            standard_parallel = 52.5
            longitude_of_central_meridian = 0.
            latitude_of_projection_origin = 52.5
            longitudeOfFirstGridPointInDegrees = -7.89
            latitudeOfFirstGridPointInDegrees = 42.935
=>   RESULT:
            #
            # gridID 2
            #
            gridtype  = projection
            gridsize  = 156025
            xsize     = 395
            ysize     = 395
            xunits    = "m"
            yunits    = "m"
            xfirst    = 0
            xinc      = 5000
            yfirst    = 0
            yinc      = 5000
            grid_mapping = Lambert_Conformal
            grid_mapping_name = lambert_conformal_conic
            standard_parallel = 52.5
            longitude_of_central_meridian = 0.
            latitude_of_projection_origin = 52.5
            longitudeOfFirstGridPointInDegrees = -7.89
            latitudeOfFirstGridPointInDegrees = 42.935
*/
  if ( cdoDebugExt )
    cdoPrint("%s(gridSrcID=%d, sampleFactor=%d) ...", __func__, gridSrcID, sampleFactor);

  int gridtype = gridInqType(gridSrcID);
  if ( ! (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION ||
          gridtype == GRID_CURVILINEAR || gridtype == GRID_GENERIC) )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridtype));
    
  int gridXsize = gridInqXsize(gridSrcID);
  int gridYsize = gridInqYsize(gridSrcID);

  if ( (sampleFactor<1) || (gridXsize<1) || (gridYsize<1) || (sampleFactor > (gridXsize/4) ) || (sampleFactor > (gridYsize/4)) )
    cdoAbort("%s(): Unsupported sampleFactor (%d)! Note that: gridXsize = %d, gridYsize = %d",
             __func__, sampleFactor, gridXsize, gridYsize);

  if ( cdoDebugExt>20 ) cdo_print_grid(gridSrcID, 1);

  int xsize = (gridXsize + (sampleFactor-1)) / sampleFactor; // HARM36_L25: (789 + 2-1) / 2 = 395
  int ysize = (gridYsize + (sampleFactor-1)) / sampleFactor;

  int gridID_sampled = gridCreate(gridtype, xsize*ysize);

  gridDefXsize(gridID_sampled, xsize);
  gridDefYsize(gridID_sampled, ysize);

  gridDefNP(gridID_sampled, gridInqNP(gridSrcID));
  gridDefDatatype(gridID_sampled, gridInqDatatype(gridSrcID));
  if ( gridInqUvRelativeToGrid(gridSrcID) ) gridDefUvRelativeToGrid(gridID_sampled, 1);

  grid_copy_attributes(gridSrcID, gridID_sampled);
  
  if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridSrcID, gridID_sampled);

  if ( gridInqXvals(gridSrcID, NULL) && gridInqYvals(gridSrcID, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
        {
          double *vals = (double *) Malloc(gridXsize*gridYsize*sizeof(double));
          gridInqXvals(gridSrcID, vals);
          double *pvals = vals;
          for ( int j = 0; j < gridYsize; j += sampleFactor )
            for ( int i = 0; i < gridXsize; i += sampleFactor )
              *pvals++ = vals[j*gridXsize+i];
          gridDefXvals(gridID_sampled, vals);
          
          gridInqYvals(gridSrcID, vals);
          pvals = vals;
          for ( int j = 0; j < gridYsize; j += sampleFactor )
            for ( int i = 0; i < gridXsize; i += sampleFactor )
              *pvals++ = vals[j*gridXsize+i];
          gridDefYvals(gridID_sampled, vals);
          Free(vals);
        }
      else
        {
          double *xvals = (double *) Malloc(gridXsize*sizeof(double));
          gridInqXvals(gridSrcID, xvals);
          for ( int i = 0, j = 0; i < gridXsize; i += sampleFactor ) xvals[j++] = xvals[i];
          gridDefXvals(gridID_sampled, xvals);
          Free(xvals);

          double *yvals = (double *) Malloc(gridYsize*sizeof(double));
          gridInqYvals(gridSrcID, yvals);
          for ( int i = 0, j = 0; i < gridYsize; i += sampleFactor ) yvals[j++] = yvals[i];
          gridDefYvals(gridID_sampled, yvals);
          Free(yvals);
        }
    }

  if ( cdoDebugExt>20 )
    {
      cdoPrint("cdo SampleGrid: define_sample_grid(): ");
      cdo_print_grid(gridID_sampled, 1);
    }

  return gridID_sampled;
}


/*
@Function  cdo_define_subgrid_grid
@Title     Define a sub-grid of another grid (LCC)

@Prototype int cdo_define_subgrid_grid(int gridIDsrc, int subI0, int subI1, int subJ0, int subJ1)
@Parameter
    @Item  gridSrcID                    Source grid
    @Item  subI0,subI1, subJ0, subJ1    Sub-grid indices

@Description
The function @func{cdo_define_subgrid_grid} defines a sub-grid of another grid (LCC)

@EndFunction
*/
int cdo_define_subgrid_grid(int gridSrcID, int subI0, int subI1, int subJ0, int subJ1)
{
/* Example of horizontal grids (Harmonie HARM36_L25):
            #
            # gridID 2
            #
            gridtype  = projection
            gridsize  = 622521
            xsize     = 789
            ysize     = 789
            xunits    = "m"
            yunits    = "m"
            xfirst    = 0
            xinc      = 2500
            yfirst    = 0
            yinc      = 2500
            grid_mapping = Lambert_Conformal
            grid_mapping_name = lambert_conformal_conic
            standard_parallel = 52.5
            longitude_of_central_meridian = 0.
            latitude_of_projection_origin = 52.5
            longitudeOfFirstGridPointInDegrees = -7.89
            latitudeOfFirstGridPointInDegrees = 42.935
=>   RESULT:
            #
            # gridID 2
            #
            gridtype  = projection
            gridsize  = 156025
            xsize     = 350
            ysize     = 350
            xunits    = "m"
            yunits    = "m"
            xfirst    = 0
            xinc      = 2500
            yfirst    = 0
            yinc      = 2500
            grid_mapping = Lambert_Conformal
            grid_mapping_name = lambert_conformal_conic
            standard_parallel = 52.5
            longitude_of_central_meridian = 0.
            latitude_of_projection_origin = 52.5
            longitudeOfFirstGridPointInDegrees = ...
            latitudeOfFirstGridPointInDegrees = ...
*/
  if ( cdoDebugExt )
    cdoPrint("%s(gridSrcID=%d, (subI0,subI1,subJ0,subJ1) = (%d,%d,%d,%d) ...",
             __func__, gridSrcID, subI0,subI1, subJ0, subJ1 );

  int gridXsize = gridInqXsize(gridSrcID);
  int gridYsize = gridInqYsize(gridSrcID);
  int maxIndexI = gridXsize-1;
  int maxIndexJ = gridYsize-1;

  if ( (subI0<0)      || (subI0>maxIndexI) ||
       (subI1<=subI0) || (subI1>maxIndexI) ||
       (subJ0<0)      || (subJ0>maxIndexJ) ||
       (subJ1<=subJ0) || (subJ1>maxIndexJ) )
    cdoAbort("%s() Incorrect subgrid specified!  (subI0,subI1,subJ0,subJ1) =(%d,%d,%d,%d) Note that: gridXsize = %d, gridYsize = %d", __func__, subI0,subI1, subJ0, subJ1, gridXsize, gridYsize);

  int gridtype = gridInqType(gridSrcID);
  if ( !(gridtype == GRID_PROJECTION && gridInqProjType(gridSrcID) == CDI_PROJ_LCC) )
    cdoAbort("%s() Error; Only LCC grid is supported; use selindexbox!", __func__);

  double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;

  gridInqParamLCC(gridSrcID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);
  x_0 = grid_missval;
  y_0 = grid_missval;

  if ( cdoDebugExt>20 ) cdo_print_grid(gridSrcID, 1);

  if ( cdoDebugExt )
    {
      cdoPrint("%s() Original LCC grid:", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridXsize, gridYsize);
      cdoPrint("xval_0 %4.3f, yval_0 %4.3f", xval_0, yval_0);
    }
  
  int gridIDcurvl = gridToCurvilinear(gridSrcID, 1);

  xval_0 = gridInqXval(gridIDcurvl, 0);
  yval_0 = gridInqYval(gridIDcurvl, 0);

  if ( cdoDebugExt )
    {
      cdoPrint("%s() Original LCC grid as curvilinear (with lats-lons computed):", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridIDcurvl), gridInqYsize(gridIDcurvl));
      cdoPrint("grid Xfirst  %4.3f, grid Yfirst  %4.3f", gridInqXval(gridIDcurvl, 0), gridInqYval(gridIDcurvl, 0));
      cdoPrint("grid Xlast   %4.3f, grid Ylast   %4.3f", gridInqXval(gridIDcurvl, gridInqSize(gridIDcurvl) -1), gridInqYval(gridIDcurvl, gridInqSize(gridIDcurvl) -1));
      cdoPrint("xval_0 %4.3f, yval_0 %4.3f", xval_0, yval_0);
    }

  int xsize = subI1 - subI0 + 1;
  int ysize = subJ1 - subJ0 + 1;

  int gridID_sampled = gridCreate(gridtype, xsize*ysize);
  
  gridDefXsize(gridID_sampled, xsize);
  gridDefYsize(gridID_sampled, ysize);
  
  if ( gridInqXvals(gridSrcID, NULL) && gridInqYvals(gridSrcID, NULL) )
    {
      double *xvals = (double*) Malloc(gridXsize*sizeof(double));
      double *yvals = (double*) Malloc(gridYsize*sizeof(double));
      gridInqXvals(gridSrcID, xvals);
      gridInqYvals(gridSrcID, yvals);
      gridDefXvals(gridID_sampled, xvals);
      gridDefYvals(gridID_sampled, yvals);
      Free(xvals);
      Free(yvals);
    }

  gridDefNP(gridID_sampled, gridInqNP(gridSrcID));
  gridDefDatatype(gridID_sampled, gridInqDatatype(gridSrcID));
  if ( gridInqUvRelativeToGrid(gridSrcID) ) gridDefUvRelativeToGrid(gridID_sampled, 1);

  grid_copy_attributes(gridSrcID, gridID_sampled);

  xval_0 = gridInqXval(gridIDcurvl, subJ0*gridXsize + subI0);
  yval_0 = gridInqYval(gridIDcurvl, subJ0*gridXsize + subI0);

  if ( cdoDebugExt )
    {
      cdoPrint("%s()  Sub-grid:", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridID_sampled), gridInqYsize(gridID_sampled));
      cdoPrint("xval_0 %4.3f, yval_0 %4.3f", xval_0, yval_0);
    }

  gridDefParamLCC(gridID_sampled, grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0);
    
  gridDestroy(gridIDcurvl);

  if ( cdoDebugExt>20 )
    {
      cdoPrint("%s(): ", __func__);
      cdo_print_grid(gridID_sampled, 1);
    }
    
  return gridID_sampled;
}
