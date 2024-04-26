#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

#include <cdi.h>
#include "cdo_int.h"
#include "griddes.h"


#if defined(HAVE_LIBNETCDF)
static void nce(int istat)
{
  // This routine provides a simple interface to NetCDF error message routine.
  if ( istat != NC_NOERR ) cdoAbort(nc_strerror(istat));
}
#endif


int cdf_openread(const char *filename)
{
  int fileID = -1;
#if defined(HAVE_LIBNETCDF)
  int nc_file_id;      /* NetCDF grid file id           */

  openLock();
  int istat = nc_open(filename, NC_NOWRITE, &nc_file_id);
  openUnlock();
  if ( istat != NC_NOERR )      
    cdoAbort("Open failed on %s! %s", filename, nc_strerror(istat));
  fileID = nc_file_id;
#else
  cdoWarning("NetCDF support not compiled in!");
#endif

  return fileID;
}


int gridFromNCfile(const char *gridfile)
{
  int gridID = -1;
#if defined(HAVE_LIBNETCDF)
  int nc_file_id;      /* NetCDF grid file id           */
  int nc_gridsize_id;  /* NetCDF grid size dim id       */
  int nc_gridcorn_id;  /* NetCDF grid corner dim id     */
  int nc_gridrank_id;  /* NetCDF grid rank dim id       */
  int nc_griddims_id;  /* NetCDF grid dimension size id */
  int nc_gridclat_id;  /* NetCDF grid corner lat var id */
  int nc_gridclon_id;  /* NetCDF grid corner lon var id */
  int nc_gridlat_id;   /* NetCDF grid center lat var id */
  int nc_gridlon_id;   /* NetCDF grid center lon var id */
  int nc_gridmask_id;  /* NetCDF grid mask id           */

  nc_type xtype;
  size_t attlen;
  size_t grid_rank, grid_size, grid_nvertex;
  int grid_dims[2];
  griddes_t grid;


  gridInit(&grid);

  /* open grid file and read grid size/name data */
  
  nc_file_id = cdf_openread(gridfile);

  if ( nc_inq_dimid(nc_file_id, "grid_size", &nc_gridsize_id)    == NC_NOERR &&
       nc_inq_dimid(nc_file_id, "grid_rank", &nc_gridrank_id)    == NC_NOERR &&
       nc_inq_dimid(nc_file_id, "grid_corners", &nc_gridcorn_id) == NC_NOERR )
    {
      nce(nc_inq_dimlen(nc_file_id, nc_gridsize_id, &grid_size)); grid.size = (int) grid_size;
      nce(nc_inq_dimlen(nc_file_id, nc_gridrank_id, &grid_rank));
      nce(nc_inq_dimlen(nc_file_id, nc_gridcorn_id, &grid_nvertex)); grid.nvertex = (int) grid_nvertex;
  
      /* check variables */
      if ( nc_inq_varid(nc_file_id, "grid_dims", &nc_griddims_id)       != NC_NOERR || 
	   nc_inq_varid(nc_file_id, "grid_center_lat", &nc_gridlat_id)  != NC_NOERR || 
	   nc_inq_varid(nc_file_id, "grid_center_lon", &nc_gridlon_id)  != NC_NOERR || 
	   nc_inq_varid(nc_file_id, "grid_corner_lat", &nc_gridclat_id) != NC_NOERR || 
	   nc_inq_varid(nc_file_id, "grid_corner_lon", &nc_gridclon_id) != NC_NOERR ) return gridID;

      nce(nc_get_var_int(nc_file_id, nc_griddims_id, grid_dims));

      if ( grid_rank == 1 )
	{
	  grid.type = GRID_UNSTRUCTURED;
	  if ( (size_t)grid_dims[0] != grid_size ) return gridID;
	}
      else
	{
	  grid.type = GRID_CURVILINEAR;
	  if ( grid.nvertex != 4 )
	    cdoAbort("curvilinear grid with %d corners unsupported", grid.nvertex);

	  grid.xsize = grid_dims[0];
	  grid.ysize = grid_dims[1];
	  if ( (size_t)grid_dims[0]*grid_dims[1] != grid_size ) return gridID;
	}

      /* allocate grid coordinates and read data */

      grid.xvals   = (double*) Malloc(grid.size*sizeof(double));
      grid.yvals   = (double*) Malloc(grid.size*sizeof(double));
      grid.xbounds = (double*) Malloc(grid.nvertex*grid.size*sizeof(double));
      grid.ybounds = (double*) Malloc(grid.nvertex*grid.size*sizeof(double));

      nce(nc_inq_vartype(nc_file_id, nc_gridlat_id, &xtype));
      grid.datatype = (xtype == NC_FLOAT) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;

      nce(nc_get_var_double(nc_file_id, nc_gridlon_id, grid.xvals));
      nce(nc_get_var_double(nc_file_id, nc_gridlat_id, grid.yvals));
      nce(nc_get_var_double(nc_file_id, nc_gridclon_id, grid.xbounds));
      nce(nc_get_var_double(nc_file_id, nc_gridclat_id, grid.ybounds));

      nce(nc_inq_attlen(nc_file_id, nc_gridlon_id, "units", &attlen));
      nce(nc_get_att_text(nc_file_id, nc_gridlon_id, "units", grid.xunits));
      grid.xunits[attlen] = 0;
      nce(nc_inq_attlen(nc_file_id, nc_gridlat_id, "units", &attlen));
      nce(nc_get_att_text(nc_file_id, nc_gridlat_id, "units", grid.yunits));
      grid.yunits[attlen] = 0;

      if ( nc_inq_varid(nc_file_id, "grid_imask", &nc_gridmask_id) == NC_NOERR )
	{
	  int i;
	  grid.mask = (int*) Malloc(grid.size*sizeof(int));
	  nce(nc_get_var_int(nc_file_id, nc_gridmask_id, grid.mask));
	  for ( i = 0; i < grid.size; ++i )
	    if ( grid.mask[i] != 1 ) break;

	  if ( i == grid.size )
	    {
	      Free(grid.mask);
	      grid.mask = NULL;
	    }
	}

      gridID = gridDefine(grid);
    }

  nce(nc_close(nc_file_id));

#else
  cdoWarning("NetCDF support not compiled in!");
#endif

  return gridID;
}


void writeNCgrid(const char *gridfile, int gridID, int *grid_imask)
{
#if defined(HAVE_LIBNETCDF)
  int nc_file_id;      /* NetCDF grid file id           */
  int nc_gridsize_id;  /* NetCDF grid size dim id       */
  int nc_gridcorn_id;  /* NetCDF grid corner dim id     */
  int nc_gridrank_id;  /* NetCDF grid rank dim id       */
  int nc_griddims_id;  /* NetCDF grid dimension size id */
  int nc_gridclat_id;  /* NetCDF grid corner lat var id */
  int nc_gridclon_id;  /* NetCDF grid corner lon var id */
  int nc_gridlat_id;   /* NetCDF grid center lat var id */
  int nc_gridlon_id;   /* NetCDF grid center lon var id */
  int nc_gridxsize_id, nc_gridysize_id, nc_grdimask_id;

  size_t grid_rank = 0, len;
  int grid_dims[2];
  int nc_dims_id[3], ndims;
  double *vals;
  char units[CDI_MAX_NAME];


  int gridtype = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);

  nc_type xtype = (gridInqDatatype(gridID) == CDI_DATATYPE_FLT64) ? NC_DOUBLE : NC_FLOAT;

  if ( gridtype == GRID_CURVILINEAR )
    {
      grid_rank = 2;
      grid_dims[0] = gridInqXsize(gridID);
      grid_dims[1] = gridInqYsize(gridID);
    }
  else if ( gridtype == GRID_UNSTRUCTURED )
    {
      grid_rank = 1;
      grid_dims[0] = gridInqSize(gridID);
    }
  else
    {
    }

  gridInqYunits(gridID, units);
  if      ( memcmp(units, "degrees", 7) == 0 )
    units[7] = 0;
  else if ( memcmp(units, "radian", 6) == 0 )
    units[6] = 0;
  else
    cdoWarning("Unknown units supplied for grid!");

  /* create NetCDF dataset for this grid */
  
  nce(nc_create(gridfile, NC_CLOBBER, &nc_file_id));

  len = strlen(gridfile);
  if ( gridfile[len-2] == 'n' && gridfile[len-1] == 'c' ) len -= 3;
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "title", len, gridfile));

  if ( CDO_Version_Info )
    nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "CDO", (int)strlen(cdoComment())+1, cdoComment()));

  /* define grid size dimension */

  nce(nc_def_dim (nc_file_id, "grid_size", gridsize, &nc_gridsize_id));
  if ( gridtype == GRID_CURVILINEAR )
    {
      nce(nc_def_dim(nc_file_id, "grid_xsize", gridInqXsize(gridID), &nc_gridxsize_id));
      nce(nc_def_dim(nc_file_id, "grid_ysize", gridInqYsize(gridID), &nc_gridysize_id));
    }

  /* define grid corner dimension */

  nce(nc_def_dim(nc_file_id, "grid_corners", gridInqNvertex(gridID), &nc_gridcorn_id));

  /* define grid rank dimension */

  nce(nc_def_dim(nc_file_id, "grid_rank", grid_rank, &nc_gridrank_id));

  /* define grid dimension size array */

  nce(nc_def_var(nc_file_id, "grid_dims", NC_INT, 1, &nc_gridrank_id, &nc_griddims_id));

  /* define grid center latitude array */

  if ( gridtype == GRID_CURVILINEAR )
    {
      ndims = 2;
      nc_dims_id[0] = nc_gridysize_id;
      nc_dims_id[1] = nc_gridxsize_id;
    }
  else
    {
      ndims = 1;
      nc_dims_id[0] = nc_gridsize_id;
    }

  nce(nc_def_var(nc_file_id, "grid_center_lat", xtype, ndims, nc_dims_id, &nc_gridlat_id));

  nce(nc_put_att_text(nc_file_id, nc_gridlat_id, "units", strlen(units), units));
  nce(nc_put_att_text(nc_file_id, nc_gridlat_id, "bounds", 15, "grid_corner_lat"));

  /* define grid center longitude array */

  nce(nc_def_var(nc_file_id, "grid_center_lon", xtype, ndims, nc_dims_id, &nc_gridlon_id));
 
  nce(nc_put_att_text(nc_file_id, nc_gridlon_id, "units", strlen(units), units));
  nce(nc_put_att_text(nc_file_id, nc_gridlon_id, "bounds", 15, "grid_corner_lon"));

  /* define grid mask */

  nce(nc_def_var(nc_file_id, "grid_imask", NC_INT, ndims, nc_dims_id, &nc_grdimask_id));

  nce(nc_put_att_text(nc_file_id, nc_grdimask_id, "units", 8, "unitless"));
  nce(nc_put_att_text(nc_file_id, nc_grdimask_id, "coordinates", 31, "grid_center_lon grid_center_lat"));

  /* define grid corner latitude array */

  if ( gridtype == GRID_CURVILINEAR )
    {
      ndims = 3;
      nc_dims_id[0] = nc_gridysize_id;
      nc_dims_id[1] = nc_gridxsize_id;
      nc_dims_id[2] = nc_gridcorn_id;
    }
  else
    {
      ndims = 2;
      nc_dims_id[0] = nc_gridsize_id;
      nc_dims_id[1] = nc_gridcorn_id;
    }

  nce(nc_def_var(nc_file_id, "grid_corner_lat", xtype, ndims, nc_dims_id, &nc_gridclat_id));

  nce(nc_put_att_text(nc_file_id, nc_gridclat_id, "units", strlen(units), units));

  /* define grid corner longitude array */

  nce(nc_def_var(nc_file_id, "grid_corner_lon", xtype, ndims, nc_dims_id, &nc_gridclon_id));

  nce(nc_put_att_text(nc_file_id, nc_gridclon_id, "units", strlen(units), units));

  /* end definition stage */

  nce(nc_enddef(nc_file_id));

  /*  write grid data */

  nce(nc_put_var_int(nc_file_id, nc_griddims_id, grid_dims));

  nce(nc_put_var_int(nc_file_id, nc_grdimask_id, grid_imask));

  vals = (double*) Malloc(gridInqNvertex(gridID)*gridsize*sizeof(double));

  gridInqYvals(gridID, vals);
  nce(nc_put_var_double(nc_file_id, nc_gridlat_id, vals));

  gridInqXvals(gridID, vals);
  nce(nc_put_var_double(nc_file_id, nc_gridlon_id, vals));

  gridInqYbounds(gridID, vals);
  nce(nc_put_var_double(nc_file_id, nc_gridclat_id, vals));

  gridInqXbounds(gridID, vals);
  nce(nc_put_var_double(nc_file_id, nc_gridclon_id, vals));

  Free(vals);

  nce(nc_close(nc_file_id));

#else
  cdoAbort("NetCDF support not compiled in!");
#endif
}
