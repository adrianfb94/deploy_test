#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"

void check_lon_range(long nlons, double *lons);
void check_lat_range(long nlats, double *lats);

typedef struct
{
  int gridID;
  long size;
  long num_cell_corners;
  double *cell_corner_lon;
  double *cell_corner_lat;
} grid_type;

typedef struct
{
  grid_type *src_grid;
  grid_type *tgt_grid;
  float *src_cell_bound_box;
} cellsearch_type;


static
grid_type *grid_new(int gridID, const char *txt)
{
  bool lgrid_destroy = false;

  if ( gridInqType(gridID) == GRID_GME )
    {
      lgrid_destroy = true;
      int gridID_gme = gridToUnstructured(gridID, 1);
      gridCompress(gridID_gme);
      gridID = gridID_gme;
    }

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    {
      lgrid_destroy = true;
      gridID = gridToCurvilinear(gridID, 1);
    }

  if ( gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0 )
    cdoAbort("%s grid corner missing!", txt);

  grid_type *grid = (grid_type*) Malloc(sizeof(grid_type));

  grid->gridID = gridID;
  grid->size = gridInqSize(grid->gridID);
  grid->num_cell_corners = (gridInqType(grid->gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(grid->gridID) : 4;

  printf("%s grid size %ld nv %ld\n", txt, grid->size, grid->num_cell_corners);
  grid->cell_corner_lon = (double*) Malloc(grid->num_cell_corners*grid->size*sizeof(double));
  grid->cell_corner_lat = (double*) Malloc(grid->num_cell_corners*grid->size*sizeof(double));
  gridInqXbounds(grid->gridID, grid->cell_corner_lon);
  gridInqYbounds(grid->gridID, grid->cell_corner_lat);

  char xunits[CDI_MAX_NAME]; xunits[0] = 0;
  char yunits[CDI_MAX_NAME]; yunits[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  grid_to_radian(xunits, grid->num_cell_corners*grid->size, grid->cell_corner_lon, "grid corner lon"); 
  grid_to_radian(yunits, grid->num_cell_corners*grid->size, grid->cell_corner_lat, "grid corner lat"); 

  // check_lon_range(grid->num_cell_corners*grid->size, grid->cell_corner_lon);
  // check_lat_range(grid->num_cell_corners*grid->size, grid->cell_corner_lat);

  if ( lgrid_destroy ) gridDestroy(gridID);

  return grid;
}

static
void grid_delete(grid_type *grid)
{
  if ( grid->cell_corner_lon ) Free(grid->cell_corner_lon);
  if ( grid->cell_corner_lat ) Free(grid->cell_corner_lat);
  if ( grid ) Free(grid);
}

void boundbox_from_corners1r(long ic, long nc, const double *restrict corner_lon,
			     const double *restrict corner_lat, float *restrict bound_box)
{
  long inc = ic*nc;

  float clat = corner_lat[inc];
  float clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for ( long j = 1; j < nc; ++j )
    {
      clat = corner_lat[inc+j];
      clon = corner_lon[inc+j];

      if ( clat < bound_box[0] ) bound_box[0] = clat;
      if ( clat > bound_box[1] ) bound_box[1] = clat;
      if ( clon < bound_box[2] ) bound_box[2] = clon;
      if ( clon > bound_box[3] ) bound_box[3] = clon;
    }

  /*
  if ( fabs(bound_box[3] - bound_box[2]) > PI )
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }
  */
}

void boundbox_from_corners(long size, long nc, const double *restrict corner_lon,
			   const double *restrict corner_lat, float *restrict bound_box)
{
  long i4, inc, j;
  double clon, clat;

  for ( long i = 0; i < size; ++i )
    {
      i4 = i<<2; // *4
      inc = i*nc;
      clat = corner_lat[inc];
      clon = corner_lon[inc];
      bound_box[i4  ] = clat;
      bound_box[i4+1] = clat;
      bound_box[i4+2] = clon;
      bound_box[i4+3] = clon;
      for ( j = 1; j < nc; ++j )
	{
	  clat = corner_lat[inc+j];
	  clon = corner_lon[inc+j];
	  if ( clat < bound_box[i4  ] ) bound_box[i4  ] = clat;
	  if ( clat > bound_box[i4+1] ) bound_box[i4+1] = clat;
	  if ( clon < bound_box[i4+2] ) bound_box[i4+2] = clon;
	  if ( clon > bound_box[i4+3] ) bound_box[i4+3] = clon;
	}
    }
}

static
cellsearch_type *cellsearch_new(grid_type *src_grid, grid_type *tgt_grid)
{
  cellsearch_type *cellsearch = (cellsearch_type*) Malloc(sizeof(cellsearch_type));

  cellsearch->src_grid = src_grid;
  cellsearch->tgt_grid = tgt_grid;

  float *src_cell_bound_box = (float*) Malloc(4*src_grid->size*sizeof(double));

  boundbox_from_corners(src_grid->size, src_grid->num_cell_corners, 
                        src_grid->cell_corner_lon, src_grid->cell_corner_lat, src_cell_bound_box);

  cellsearch->src_cell_bound_box = src_cell_bound_box;

  return cellsearch;
}

static
void cellsearch_delete(cellsearch_type *cellsearch)
{
  if ( cellsearch ) Free(cellsearch);
}

static
long search_cells(cellsearch_type *cellsearch, long tgt_cell_add, long *srch_add)
{
  grid_type *src_grid = cellsearch->src_grid;
  grid_type *tgt_grid = cellsearch->tgt_grid;
  float *src_cell_bound_box = cellsearch->src_cell_bound_box;
  float tgt_cell_bound_box[4];

  boundbox_from_corners1r(tgt_cell_add, tgt_grid->num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat, tgt_cell_bound_box);

  long src_cell_addm4;

  float bound_box_lat1 = tgt_cell_bound_box[0];
  float bound_box_lat2 = tgt_cell_bound_box[1];
  float bound_box_lon1 = tgt_cell_bound_box[2];
  float bound_box_lon2 = tgt_cell_bound_box[3];

  long num_srch_cells = 0;
  for ( long src_cell_add = 0; src_cell_add < src_grid->size; ++src_cell_add )
    {
      src_cell_addm4 = src_cell_add<<2;
      if ( (src_cell_bound_box[src_cell_addm4+2] <= bound_box_lon2)  &&
	   (src_cell_bound_box[src_cell_addm4+3] >= bound_box_lon1) )
	{
	  if ( (src_cell_bound_box[src_cell_addm4  ] <= bound_box_lat2)  &&
	       (src_cell_bound_box[src_cell_addm4+1] >= bound_box_lat1) )
	    {
	      srch_add[num_srch_cells] = src_cell_add;
	      num_srch_cells++;
	    }
	}
    }

  return num_srch_cells;
}      

static
void cell_search(int gridIDsrc, int gridIDtgt)
{
  grid_type *src_grid = grid_new(gridIDsrc, "source");
  grid_type *tgt_grid = grid_new(gridIDtgt, "target");

  long *srch_add = (long*) Malloc(src_grid->size*sizeof(long));

  cellsearch_type *cellsearch = cellsearch_new(src_grid, tgt_grid);
  
  for ( long tgt_cell_add = 0; tgt_cell_add < tgt_grid->size; ++tgt_cell_add )
    {
      long num_srch_cells = search_cells(cellsearch, tgt_cell_add, srch_add);

      if ( cdoVerbose && num_srch_cells > 0 )
        {
          printf("tgt cell %ld: found %ld src cells\n", tgt_cell_add, num_srch_cells);
          for ( long n = 0; n < num_srch_cells; ++n ) printf("   %ld: %ld\n", n+1, srch_add[n]);
        }
    }

  cellsearch_delete(cellsearch);
  grid_delete(src_grid);
  grid_delete(tgt_grid);
  Free(srch_add);
}

static
void point_search(int gridIDsrc, int gridIDtgt)
{
}


void *Gridsearch(void *argument)
{
  cdoInitialize(argument);

  // clang-format off
  int PSEARCH = cdoOperatorAdd("testpointsearch",  0,   0, NULL);
  int CSEARCH = cdoOperatorAdd("testcellsearch",   0,   0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg("source and target grid description file or name");
  operatorCheckArgc(2);
  
  int gridID1 = cdoDefineGrid(operatorArgv()[0]);
  int gridID2 = cdoDefineGrid(operatorArgv()[1]);

  if      ( operatorID == CSEARCH ) cell_search(gridID1, gridID2);
  else if ( operatorID == CSEARCH ) point_search(gridID1, gridID2);

  cdoFinish();

  return 0;
}
