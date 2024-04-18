#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <string.h>
#include <float.h>  /* FLT_EPSILON */
#include <limits.h> /* INT_MAX     */

#include "dmemory.h"
#include "cdi.h"
#include "cdi_cksum.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "grid.h"
#include "gaussgrid.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"
#include "vlist.h"

double grid_missval = -9999.;
int (*proj_lonlat_to_lcc_func)() = NULL;
int (*proj_lcc_to_lonlat_func)() = NULL;

/* the value in the second pair of brackets must match the length of
 * the longest string (including terminating NUL) */
static const char Grids[][17] = {
  /*  0 */  "undefined",
  /*  1 */  "generic",
  /*  2 */  "gaussian",
  /*  3 */  "gaussian reduced",
  /*  4 */  "lonlat",
  /*  5 */  "spectral",
  /*  6 */  "fourier",
  /*  7 */  "gme",
  /*  8 */  "trajectory",
  /*  9 */  "unstructured",
  /* 10 */  "curvilinear",
  /* 11 */  "lcc",
  /* 12 */  "projection",
  /* 13 */  "characterXY",
};

/* must match table below */
enum xystdname_idx {
  grid_xystdname_grid_latlon,
  grid_xystdname_latlon,
  grid_xystdname_projection,
  grid_xystdname_char,
};
static const char xystdname_tab[][2][24] = {
  [grid_xystdname_grid_latlon] = { "grid_longitude",
                                   "grid_latitude" },
  [grid_xystdname_latlon] = { "longitude",
                              "latitude" },
  [grid_xystdname_projection] = { "projection_x_coordinate",
                                  "projection_y_coordinate" },
  [grid_xystdname_char] = { "region",
                            "region" },
};



static int    gridCompareP    ( void * gridptr1, void * gridptr2 );
static void   gridDestroyP    ( void * gridptr );
static void   gridPrintP      ( void * gridptr, FILE * fp );
static int    gridGetPackSize ( void * gridptr, void *context);
static void   gridPack        ( void * gridptr, void * buff, int size, int *position, void *context);
static int    gridTxCode      ( void );

static const resOps gridOps = {
  gridCompareP,
  gridDestroyP,
  gridPrintP,
  gridGetPackSize,
  gridPack,
  gridTxCode
};

static int  GRID_Debug = 0;   /* If set to 1, debugging */


grid_t *grid_to_pointer(int gridID)
{
  return (grid_t *)reshGetVal(gridID, &gridOps);
}

#define gridMark4Update(gridID) reshSetStatus(gridID, &gridOps, RESH_DESYNC_IN_USE)

static
bool cdiInqAttConvertedToFloat(int gridID, int atttype, const char *attname, int attlen, double *attflt)
{
  bool status = true;

  if ( atttype == CDI_DATATYPE_INT32 )
    {
      int attint[attlen];
      cdiInqAttInt(gridID, CDI_GLOBAL, attname, attlen, attint);
      for ( int i = 0; i < attlen; ++i ) attflt[i] = (double)attint[i];
    }
  else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
    {
      cdiInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, attflt);
    }
  else
    {
      status = false;
    }

  return status;
}


void grid_init(grid_t *gridptr)
{
  gridptr->self          = CDI_UNDEFID;
  gridptr->type          = CDI_UNDEFID;
  gridptr->proj          = CDI_UNDEFID;
  gridptr->projtype      = CDI_UNDEFID;
  gridptr->mask          = NULL;
  gridptr->mask_gme      = NULL;
  gridptr->x.vals        = NULL;
  gridptr->x.cvals       = NULL;
  gridptr->x.clength     = 0;
  gridptr->y.vals        = NULL;
  gridptr->y.cvals       = NULL;
  gridptr->y.clength     = 0;
  gridptr->x.bounds      = NULL;
  gridptr->y.bounds      = NULL;
  gridptr->area          = NULL;
  gridptr->rowlon        = NULL;
  gridptr->nrowlon       = 0;

  gridptr->x.first       = 0.0;
  gridptr->x.last        = 0.0;
  gridptr->x.inc         = 0.0;
  gridptr->y.first       = 0.0;
  gridptr->y.last        = 0.0;
  gridptr->y.inc         = 0.0;

  gridptr->gme.nd        = 0;
  gridptr->gme.ni        = 0;
  gridptr->gme.ni2       = 0;
  gridptr->gme.ni3       = 0;

  gridptr->trunc         = 0;
  gridptr->nvertex       = 0;
  gridptr->number        = 0;
  gridptr->position      = 0;
  gridptr->reference     = NULL;
  gridptr->datatype      = 0;
  gridptr->size          = 0;
  gridptr->x.size        = 0;
  gridptr->y.size        = 0;
  gridptr->np            = 0;
  gridptr->x.flag        = 0;
  gridptr->y.flag        = 0;
  gridptr->isCyclic      = CDI_UNDEFID;

  gridptr->lcomplex      = false;
  gridptr->hasdims       = true;
  gridptr->x.dimname[0]  = 0;
  gridptr->y.dimname[0]  = 0;
  gridptr->x.name[0]     = 0;
  gridptr->y.name[0]     = 0;
  gridptr->x.longname[0] = 0;
  gridptr->y.longname[0] = 0;
  gridptr->x.units[0]    = 0;
  gridptr->y.units[0]    = 0;
  gridptr->x.stdname     = NULL;
  gridptr->y.stdname     = NULL;
  gridptr->vdimname[0]   = 0;
  gridptr->mapname[0]    = 0;
  gridptr->mapping[0]    = 0;
  memset(gridptr->uuid, 0, CDI_UUID_SIZE);
  gridptr->name          = NULL;
  gridptr->vtable        = &cdiGridVtable;
  gridptr->atts.nalloc   = MAX_ATTRIBUTES;
  gridptr->atts.nelems   = 0;
  gridptr->uvRelativeToGrid      = 0;   // Some models deliver wind U,V relative to the grid-cell
  gridptr->iScansNegatively      = 0;
  gridptr->jScansPositively      = 1;
  gridptr->jPointsAreConsecutive = 0;
  gridptr->scanningMode          = 128*gridptr->iScansNegatively + 64*gridptr->jScansPositively + 32*gridptr->jPointsAreConsecutive;
  /* scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
               64  = 128 * 0                + 64 *        1         + 32 * 0
               00  = 128 * 0                + 64 *        0         + 32 * 0
               96  = 128 * 0                + 64 *        1         + 32 * 1
     Default / implicit scanning mode is 64:
                        i and j scan positively, i points are consecutive (row-major)        */
}


static
void grid_free_components(grid_t *gridptr)
{
  void *p2free[] = { gridptr->mask, gridptr->mask_gme,
                     gridptr->x.vals, gridptr->y.vals,
                     gridptr->x.cvals, gridptr->y.cvals,
                     gridptr->x.bounds, gridptr->y.bounds,
                     gridptr->rowlon, gridptr->area,
                     gridptr->reference, gridptr->name};

  for ( size_t i = 0; i < sizeof(p2free)/sizeof(p2free[0]); ++i )
    if ( p2free[i] ) Free(p2free[i]);
}

void grid_free(grid_t *gridptr)
{
  grid_free_components(gridptr);
  grid_init(gridptr);
}

static
grid_t *gridNewEntry(cdiResH resH)
{
  grid_t *gridptr = (grid_t*) Malloc(sizeof(grid_t));
  grid_init(gridptr);

  if ( resH == CDI_UNDEFID )
    gridptr->self = reshPut(gridptr, &gridOps);
  else
    {
      gridptr->self = resH;
      reshReplace(resH, gridptr, &gridOps);
    }

  return gridptr;
}

static
void gridInit(void)
{
  static bool gridInitialized = false;
  if ( gridInitialized ) return;
  gridInitialized = true;

  const char *env = getenv("GRID_DEBUG");
  if ( env ) GRID_Debug = atoi(env);
}

static
void grid_copy_base_scalar_fields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  memcpy(gridptrDup, gridptrOrig, sizeof(grid_t));
  gridptrDup->self = CDI_UNDEFID;
  if ( gridptrOrig->reference )
    gridptrDup->reference = strdupx(gridptrOrig->reference);
}


static
grid_t *grid_copy_base(grid_t *gridptrOrig)
{
  grid_t *gridptrDup = (grid_t *)Malloc(sizeof (*gridptrDup));
  gridptrOrig->vtable->copyScalarFields(gridptrOrig, gridptrDup);
  gridptrOrig->vtable->copyArrayFields(gridptrOrig, gridptrDup);
  return gridptrDup;
}


unsigned cdiGridCount(void)
{
  return reshCountType(&gridOps);
}

static inline
void gridSetString(char *gridstrname, const char *name, size_t len)
{
  if ( len > CDI_MAX_NAME ) len = CDI_MAX_NAME;
  strncpy(gridstrname, name, len);
  gridstrname[len - 1] = 0;
}

static inline
void gridGetString(char *name, const char *gridstrname, size_t len)
{
  if ( len > CDI_MAX_NAME ) len = CDI_MAX_NAME;
  strncpy(name, gridstrname, len);
  name[len - 1] = 0;
}

static inline
void gridSetName(char *gridstrname, const char *name)
{
  strncpy(gridstrname, name, CDI_MAX_NAME);
  gridstrname[CDI_MAX_NAME - 1] = 0;
}


void cdiGridTypeInit(grid_t *gridptr, int gridtype, int size)
{
  gridptr->type = gridtype;
  gridptr->size = size;

  if      ( gridtype == GRID_CURVILINEAR  ) gridptr->nvertex = 4;
  else if ( gridtype == GRID_UNSTRUCTURED ) gridptr->x.size = size;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_GME:
      {
        if ( gridtype == GRID_TRAJECTORY )
          {
            if ( !gridptr->x.name[0] ) gridSetName(gridptr->x.name, "tlon");
            if ( !gridptr->y.name[0] ) gridSetName(gridptr->y.name, "tlat");
          }
        else
          {
            if ( !gridptr->x.name[0] ) gridSetName(gridptr->x.name, "lon");
            if ( !gridptr->y.name[0] ) gridSetName(gridptr->y.name, "lat");
          }

        if ( !gridptr->x.longname[0] ) gridSetName(gridptr->x.longname, "longitude");
        if ( !gridptr->y.longname[0] ) gridSetName(gridptr->y.longname, "latitude");

        if ( !gridptr->x.units[0] ) gridSetName(gridptr->x.units, "degrees_east");
        if ( !gridptr->y.units[0] ) gridSetName(gridptr->y.units, "degrees_north");

        gridptr->x.stdname = xystdname_tab[grid_xystdname_latlon][0];
        gridptr->y.stdname = xystdname_tab[grid_xystdname_latlon][1];

        break;
      }
    case GRID_CHARXY:
      {
        if ( gridptr->x.cvals ) gridptr->x.stdname = xystdname_tab[grid_xystdname_char][0];
        if ( gridptr->y.cvals ) gridptr->y.stdname = xystdname_tab[grid_xystdname_char][0];

        break;
     }
    case GRID_GENERIC:
    case GRID_PROJECTION:
      {
        if ( gridptr->x.name[0] == 0 ) gridSetName(gridptr->x.name, "x");
        if ( gridptr->y.name[0] == 0 ) gridSetName(gridptr->y.name, "y");
        if ( gridtype == GRID_PROJECTION )
          {
            gridSetName(gridptr->mapname, "Projection");

            gridptr->x.stdname = xystdname_tab[grid_xystdname_projection][0];
            gridptr->y.stdname = xystdname_tab[grid_xystdname_projection][1];

            if ( !gridptr->x.units[0] ) gridSetName(gridptr->x.units, "m");
            if ( !gridptr->y.units[0] ) gridSetName(gridptr->y.units, "m");
          }
        break;
      }
    }
}


// used also in CDO
void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *restrict xvals)
{
  if ( (! (fabs(xinc) > 0)) && xsize > 1 )
    {
      if ( xfirst >= xlast )
        {
          while ( xfirst >= xlast ) xlast += 360;
          xinc = (xlast-xfirst)/(xsize);
        }
      else
        {
          xinc = (xlast-xfirst)/(xsize-1);
        }
    }

  for ( int i = 0; i < xsize; ++i )
    xvals[i] = xfirst + i*xinc;
}

static
void calc_gaussgrid(double *restrict yvals, int ysize, double yfirst, double ylast)
{
  double *restrict yw = (double *) Malloc((size_t)ysize * sizeof(double));
  gaussaw(yvals, yw, (size_t)ysize);
  Free(yw);
  for (int i = 0; i < ysize; i++ )
    yvals[i] = asin(yvals[i])/M_PI*180.0;

  if ( yfirst < ylast && yfirst > -90.0 && ylast < 90.0 )
    {
      int yhsize = ysize/2;
      for (int i = 0; i < yhsize; i++ )
        {
          double ytmp = yvals[i];
          yvals[i] = yvals[ysize-i-1];
          yvals[ysize-i-1] = ytmp;
        }
    }
}

// used also in CDO
void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *restrict yvals)
{
  const double deleps = 0.002;

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
    {
      if ( ysize > 2 )
	{
	  calc_gaussgrid(yvals, ysize, yfirst, ylast);

	  if ( ! (IS_EQUAL(yfirst, 0) && IS_EQUAL(ylast, 0)) )
	    if ( fabs(yvals[0] - yfirst) > deleps || fabs(yvals[ysize-1] - ylast) > deleps )
	      {
		double *restrict ytmp = NULL;
		int nstart, lfound = 0;
		int ny = (int) (180./(fabs(ylast-yfirst)/(ysize-1)) + 0.5);
		ny -= ny%2;
		if ( ny > ysize && ny < 4096 )
		  {
		    ytmp = (double *) Malloc((size_t)ny * sizeof (double));
		    calc_gaussgrid(ytmp, ny, yfirst, ylast);
                    {
                      int i;
                      for ( i = 0; i < (ny-ysize); i++ )
                        if ( fabs(ytmp[i] - yfirst) < deleps ) break;
                      nstart = i;
                    }

		    lfound = (nstart+ysize-1) < ny
                      && fabs(ytmp[nstart+ysize-1] - ylast) < deleps;
                    if ( lfound )
                      {
                        for (int i = 0; i < ysize; i++) yvals[i] = ytmp[i+nstart];
                      }
		  }

		if ( !lfound )
		  {
		    Warning("Cannot calculate gaussian latitudes for lat1 = %g latn = %g!", yfirst, ylast);
		    for (int i = 0; i < ysize; i++ ) yvals[i] = 0;
		    yvals[0] = yfirst;
		    yvals[ysize-1] = ylast;
		  }

		if ( ytmp ) Free(ytmp);
	      }
	}
      else
        {
          yvals[0] = yfirst;
          yvals[ysize-1] = ylast;
        }
    }
  /*     else if ( gridtype == GRID_LONLAT || gridtype == GRID_GENERIC ) */
  else
    {
      if ( (! (fabs(yinc) > 0)) && ysize > 1 )
        {
          if ( IS_EQUAL(yfirst, ylast) && IS_NOT_EQUAL(yfirst, 0) ) ylast *= -1;

          if ( yfirst > ylast )
            yinc = (yfirst-ylast)/(ysize-1);
          else if ( yfirst < ylast )
            yinc = (ylast-yfirst)/(ysize-1);
          else
            {
              if ( ysize%2 != 0 )
                {
                  yinc = 180.0/(ysize-1);
                  yfirst = -90;
                }
              else
                {
                  yinc = 180.0/ysize;
                  yfirst = -90 + yinc/2;
                }
            }
        }

      if ( yfirst > ylast && yinc > 0 ) yinc = -yinc;

      for (int i = 0; i < ysize; i++ )
        yvals[i] = yfirst + i*yinc;
    }
  /*
    else
    Error("unable to calculate values for %s grid!", gridNamePtr(gridtype));
  */
}

/*
@Function  gridCreate
@Title     Create a horizontal Grid

@Prototype int gridCreate(int gridtype, int size)
@Parameter
    @Item  gridtype  The type of the grid, one of the set of predefined CDI grid types.
                     The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_GAUSSIAN},
                     @func{GRID_LONLAT}, @func{GRID_PROJECTION}, @func{GRID_SPECTRAL},
                     @func{GRID_GME}, @func{GRID_CURVILINEAR} and @func{GRID_UNSTRUCTURED}.
    @Item  size      Number of gridpoints.

@Description
The function @func{gridCreate} creates a horizontal Grid.

@Result
@func{gridCreate} returns an identifier to the Grid.

@Example
Here is an example using @func{gridCreate} to create a regular lon/lat Grid:

@Source
#include "cdi.h"
   ...
#define  nlon  12
#define  nlat   6
   ...
double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
double lats[nlat] = {-75, -45, -15, 15, 45, 75};
int gridID;
   ...
gridID = gridCreate(GRID_LONLAT, nlon*nlat);
gridDefXsize(gridID, nlon);
gridDefYsize(gridID, nlat);
gridDefXvals(gridID, lons);
gridDefYvals(gridID, lats);
   ...
@EndSource
@EndFunction
*/
int gridCreate(int gridtype, int size)
{
  if ( CDI_Debug ) Message("gridtype=%s  size=%d", gridNamePtr(gridtype), size);

  if ( size < 0 || size > INT_MAX ) Error("Grid size (%d) out of bounds (0 - %d)!", size, INT_MAX);

  gridInit();

  grid_t *gridptr = gridNewEntry(CDI_UNDEFID);
  if ( ! gridptr ) Error("No memory");

  int gridID = gridptr->self;

  if ( CDI_Debug ) Message("gridID: %d", gridID);

  cdiGridTypeInit(gridptr, gridtype, size);

  return gridID;
}

static
void gridDestroyKernel( grid_t * gridptr )
{
  xassert ( gridptr );

  int id = gridptr->self;

  grid_free_components(gridptr);
  Free( gridptr );

  reshRemove ( id, &gridOps );
}

/*
@Function  gridDestroy
@Title     Destroy a horizontal Grid

@Prototype void gridDestroy(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@EndFunction
*/
void gridDestroy(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->destroy(gridptr);
}

static
void gridDestroyP(void * gridptr)
{
  ((grid_t *)gridptr)->vtable->destroy((grid_t *)gridptr);
}


const char *gridNamePtr(int gridtype)
{
  int size = (int) (sizeof(Grids)/sizeof(Grids[0]));

  const char *name = (gridtype >= 0 && gridtype < size) ? Grids[gridtype] : Grids[GRID_GENERIC];

  return name;
}


void gridName(int gridtype, char *gridname)
{
  strcpy(gridname, gridNamePtr(gridtype));
}

static
void *grid_key_to_ptr(grid_t *gridptr, int key)
{
  void *keyptr = NULL;

  switch (key)
    {
    case CDI_KEY_XNAME:      keyptr = (void*)gridptr->x.name; break;
    case CDI_KEY_XLONGNAME:  keyptr = (void*)gridptr->x.longname; break;
    case CDI_KEY_XUNITS:     keyptr = (void*)gridptr->x.units; break;
    case CDI_KEY_YNAME:      keyptr = (void*)gridptr->y.name; break;
    case CDI_KEY_YLONGNAME:  keyptr = (void*)gridptr->y.longname; break;
    case CDI_KEY_YUNITS:     keyptr = (void*)gridptr->y.units; break;
    case CDI_KEY_XDIMNAME:   keyptr = (void*)gridptr->x.dimname; break;
    case CDI_KEY_YDIMNAME:   keyptr = (void*)gridptr->y.dimname; break;
    case CDI_KEY_VDIMNAME:   keyptr = (void*)gridptr->vdimname; break;
    case CDI_KEY_MAPPING:    keyptr = (void*)gridptr->mapname; break;
    case CDI_KEY_MAPNAME:    keyptr = (void*)gridptr->mapping; break;
    }

  return keyptr;
}

/*
@Function  cdiGridDefKeyStr
@Title     Define a CDI grid string value from a key

@Prototype int cdiGridDefKeyStr(int gridID, int key, int size, const char *mesg)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  key      The key to be searched
    @Item  size     The allocated length of the string on input
    @Item  mesg     The address of a string where the data will be read

@Description
The function @func{cdiGridDefKeyStr} defines a CDI grid string value from a key.

@Result
@func{cdiGridDefKeyStr} returns 0 if OK and integer value on error.

@EndFunction
*/
int cdiGridDefKeyStr(int gridID, int key, int size, const char *mesg)
{
  if ( size < 1 || mesg == NULL || *mesg == 0 ) return -1;

  grid_t *gridptr = grid_to_pointer(gridID);

  char *keyptr = (char*)grid_key_to_ptr(gridptr, key);
  if ( keyptr == NULL )
    {
      Warning("CDI grid string key %d not supported!", key);
      return -1;
    }

  gridSetString(keyptr, mesg, (size_t)size);
  gridMark4Update(gridID);

  return 0;
}

/*
@Function  cdiGridInqKeyStr
@Title     Get a CDI grid string value from a key

@Prototype int cdiGridInqKeyStr(int gridID, int key, int size, char *mesg)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  key      The key to be searched.
    @Item  size     The allocated length of the string on input.
    @Item  mesg     The address of a string where the data will be retrieved.
                    The caller must allocate space for the returned string.
                    The maximum possible length, in characters, of the string
                    is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{cdiGridInqKeyStr} return a CDI grid string value from a key.

@Result
@func{cdiGridInqKeyStr} returns 0 if OK and integer value on error.

@EndFunction
*/
int cdiGridInqKeyStr(int gridID, int key, int size, char *mesg)
{
  if ( size < 1 || mesg == NULL ) return -1;

  grid_t *gridptr = grid_to_pointer(gridID);
  const char *keyptr = (const char*)grid_key_to_ptr(gridptr, key);
  if ( keyptr == NULL)
    {
      Warning("CDI grid string key %d not supported!", key);
      return -1;
    }

  gridGetString(mesg, keyptr, (size_t)size);

  return 0;
}

/*
@Function  gridDefXname
@Title     Define the name of a X-axis

@Prototype void gridDefXname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the X-axis.

@Description
The function @func{gridDefXname} defines the name of a X-axis.

@EndFunction
*/
void gridDefXname(int gridID, const char *xname)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, xname);
}

/*
@Function  gridDefXlongname
@Title     Define the longname of a X-axis

@Prototype void gridDefXlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the X-axis.

@Description
The function @func{gridDefXlongname} defines the longname of a X-axis.

@EndFunction
*/
void gridDefXlongname(int gridID, const char *xlongname)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, xlongname);
}

/*
@Function  gridDefXunits
@Title     Define the units of a X-axis

@Prototype void gridDefXunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the X-axis.

@Description
The function @func{gridDefXunits} defines the units of a X-axis.

@EndFunction
*/
void gridDefXunits(int gridID, const char *xunits)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
}

/*
@Function  gridDefYname
@Title     Define the name of a Y-axis

@Prototype void gridDefYname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the Y-axis.

@Description
The function @func{gridDefYname} defines the name of a Y-axis.

@EndFunction
*/
void gridDefYname(int gridID, const char *yname)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, yname);
}

/*
@Function  gridDefYlongname
@Title     Define the longname of a Y-axis

@Prototype void gridDefYlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the Y-axis.

@Description
The function @func{gridDefYlongname} defines the longname of a Y-axis.

@EndFunction
*/
void gridDefYlongname(int gridID, const char *ylongname)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, ylongname);
}

/*
@Function  gridDefYunits
@Title     Define the units of a Y-axis

@Prototype void gridDefYunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the Y-axis.

@Description
The function @func{gridDefYunits} defines the units of a Y-axis.

@EndFunction
*/
void gridDefYunits(int gridID, const char *yunits)
{
  (void)cdiGridDefKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
}

/*
@Function  gridInqXname
@Title     Get the name of a X-axis

@Prototype void gridInqXname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  name     Name of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXname} returns the name of a X-axis.

@Result
@func{gridInqXname} returns the name of the X-axis to the parameter name.

@EndFunction
*/
void gridInqXname(int gridID, char *xname)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, xname);
}

/*
@Function  gridInqXlongname
@Title     Get the longname of a X-axis

@Prototype void gridInqXlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  longname Longname of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXlongname} returns the longname of a X-axis.

@Result
@func{gridInqXlongname} returns the longname of the X-axis to the parameter longname.

@EndFunction
*/
void gridInqXlongname(int gridID, char *xlongname)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, xlongname);
}

/*
@Function  gridInqXunits
@Title     Get the units of a X-axis

@Prototype void gridInqXunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  units    Units of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXunits} returns the units of a X-axis.

@Result
@func{gridInqXunits} returns the units of the X-axis to the parameter units.

@EndFunction
*/
void gridInqXunits(int gridID, char *xunits)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
}


void gridInqXstdname(int gridID, char *xstdname)
{
  if ( xstdname )
    {
      xstdname[0] = 0;
      grid_t *gridptr = grid_to_pointer(gridID);
      if ( gridptr->x.stdname ) strcpy(xstdname, gridptr->x.stdname);
    }
}

/*
@Function  gridInqYname
@Title     Get the name of a Y-axis

@Prototype void gridInqYname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  name     Name of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYname} returns the name of a Y-axis.

@Result
@func{gridInqYname} returns the name of the Y-axis to the parameter name.

@EndFunction
*/
void gridInqYname(int gridID, char *yname)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, yname);
}

/*
@Function  gridInqYlongname
@Title     Get the longname of a Y-axis

@Prototype void gridInqYlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  longname Longname of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYlongname} returns the longname of a Y-axis.

@Result
@func{gridInqYlongname} returns the longname of the Y-axis to the parameter longname.

@EndFunction
*/
void gridInqYlongname(int gridID, char *ylongname)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, ylongname);
}

/*
@Function  gridInqYunits
@Title     Get the units of a Y-axis

@Prototype void gridInqYunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  units    Units of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYunits} returns the units of a Y-axis.

@Result
@func{gridInqYunits} returns the units of the Y-axis to the parameter units.

@EndFunction
*/
void gridInqYunits(int gridID, char *yunits)
{
  (void)cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
}


void gridInqYstdname(int gridID, char *ystdname)
{
  if ( ystdname )
    {
      ystdname[0] = 0;
      grid_t *gridptr = grid_to_pointer(gridID);
      if ( gridptr->y.stdname ) strcpy(ystdname, gridptr->y.stdname);
    }
}


void gridDefProj(int gridID, int projID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->proj = projID;

  if ( gridptr->type == GRID_CURVILINEAR )
    {
      grid_t *projptr = grid_to_pointer(projID);
      if ( projptr->x.name[0] ) strcpy(gridptr->x.dimname, projptr->x.name);
      if ( projptr->y.name[0] ) strcpy(gridptr->y.dimname, projptr->y.name);
    }
}


int gridInqProj(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->proj;
}


int gridInqProjType(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int projtype = gridptr->projtype;

  if ( projtype == -1 )
    {
      char mapping[CDI_MAX_NAME]; mapping[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
      if ( mapping[0] )
        {
          if      ( strcmp(mapping, "rotated_latitude_longitude") == 0 )   projtype = CDI_PROJ_RLL;
          else if ( strcmp(mapping, "lambert_azimuthal_equal_area") == 0 ) projtype = CDI_PROJ_LAEA;
          else if ( strcmp(mapping, "lambert_conformal_conic") == 0 )      projtype = CDI_PROJ_LCC;
          else if ( strcmp(mapping, "sinusoidal") == 0 )                   projtype = CDI_PROJ_SINU;

          gridptr->projtype = projtype;
        }
    }

  return projtype;
}


void gridVerifyProj(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int projtype = gridInqProjType(gridID);

  if ( projtype == CDI_PROJ_RLL )
    {
      gridptr->x.stdname = xystdname_tab[grid_xystdname_grid_latlon][0];
      gridptr->y.stdname = xystdname_tab[grid_xystdname_grid_latlon][1];
      gridSetName(gridptr->x.units, "degrees");
      gridSetName(gridptr->y.units, "degrees");
    }
  else if ( projtype == CDI_PROJ_LCC )
    {
      gridptr->x.stdname = xystdname_tab[grid_xystdname_projection][0];
      gridptr->y.stdname = xystdname_tab[grid_xystdname_projection][1];
      if ( !gridptr->x.units[0] ) gridSetName(gridptr->x.units, "m");
      if ( !gridptr->y.units[0] ) gridSetName(gridptr->y.units, "m");
    }
}

/*
@Function  gridInqType
@Title     Get the type of a Grid

@Prototype int gridInqType(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqType} returns the type of a Grid.

@Result
@func{gridInqType} returns the type of the grid,
one of the set of predefined CDI grid types.
The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_GAUSSIAN},
@func{GRID_LONLAT}, @func{GRID_PROJECTION}, @func{GRID_SPECTRAL}, @func{GRID_GME},
@func{GRID_CURVILINEAR} and @func{GRID_UNSTRUCTURED}.

@EndFunction
*/
int gridInqType(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  return gridptr->type;
}


/*
@Function  gridInqSize
@Title     Get the size of a Grid

@Prototype int gridInqSize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqSize} returns the size of a Grid.

@Result
@func{gridInqSize} returns the number of grid points of a Grid.

@EndFunction
*/
int gridInqSize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int size = gridptr->size;

  if ( size == 0 )
    {
      int xsize = gridptr->x.size;
      int ysize = gridptr->y.size;

      if ( ysize )
        size = xsize * ysize;
      else
        size = xsize;

      gridptr->size = size;
    }

  return size;
}

static
int nsp2trunc(int nsp)
{
  /*  nsp = (trunc+1)*(trunc+1)              */
  /*      => trunc^2 + 3*trunc - (x-2) = 0   */
  /*                                         */
  /*  with:  y^2 + p*y + q = 0               */
  /*         y = -p/2 +- sqrt((p/2)^2 - q)   */
  /*         p = 3 and q = - (x-2)           */
  int trunc = (int) (sqrt(nsp*4 + 1.) - 3) / 2;
  return trunc;
}


int gridInqTrunc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->trunc == 0 )
    {
      if ( gridptr->type == GRID_SPECTRAL )
        gridptr->trunc = nsp2trunc(gridptr->size);
      /*
      else if      ( gridptr->type == GRID_GAUSSIAN )
        gridptr->trunc = nlat2trunc(gridptr->y.size);
      */
    }

  return gridptr->trunc;
}


void gridDefTrunc(int gridID, int trunc)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->trunc != trunc )
    {
      gridMark4Update(gridID);
      gridptr->trunc = trunc;
    }
}

/*
@Function  gridDefXsize
@Title     Define the number of values of a X-axis

@Prototype void gridDefXsize(int gridID, int xsize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xsize    Number of values of a X-axis.

@Description
The function @func{gridDefXsize} defines the number of values of a X-axis.

@EndFunction
*/
void gridDefXsize(int gridID, int xsize)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int gridSize = gridInqSize(gridID);
  if ( xsize > gridSize )
    Error("xsize %d is greater then gridsize %d", xsize, gridSize);

  int gridType = gridInqType(gridID);
  if ( gridType == GRID_UNSTRUCTURED && xsize != gridSize )
    Error("xsize %d must be equal to gridsize %d for gridtype: UNSTRUCTURED", xsize, gridSize);

  if ( gridptr->x.size != xsize )
    {
      gridMark4Update(gridID);
      gridptr->x.size = xsize;
    }

  if ( gridType != GRID_UNSTRUCTURED && gridType != GRID_PROJECTION )
    {
      long axisproduct = gridptr->x.size*gridptr->y.size;
      if ( axisproduct > 0 && axisproduct != gridSize )
        Error("Inconsistent grid declaration! (xsize=%d ysize=%d gridsize=%d)",
              gridptr->x.size, gridptr->y.size, gridSize);
    }
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefDatatype(int gridID, int datatype)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->datatype != datatype )
    {
      gridMark4Update(gridID);
      gridptr->datatype = datatype;
    }
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
int gridInqDatatype(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->datatype;
}

/*
@Function  gridInqXsize
@Title     Get the number of values of a X-axis

@Prototype int gridInqXsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqXsize} returns the number of values of a X-axis.

@Result
@func{gridInqXsize} returns the number of values of a X-axis.

@EndFunction
*/
int gridInqXsize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->x.size;
}

/*
@Function  gridDefYsize
@Title     Define the number of values of a Y-axis

@Prototype void gridDefYsize(int gridID, int ysize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ysize    Number of values of a Y-axis.

@Description
The function @func{gridDefYsize} defines the number of values of a Y-axis.

@EndFunction
*/
void gridDefYsize(int gridID, int ysize)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int gridSize = gridInqSize(gridID);

  if ( ysize > gridSize )
    Error("ysize %d is greater then gridsize %d", ysize, gridSize);

  int gridType = gridInqType(gridID);
  if ( gridType == GRID_UNSTRUCTURED && ysize != gridSize )
    Error("ysize %d must be equal gridsize %d for gridtype: UNSTRUCTURED", ysize, gridSize);

  if ( gridptr->y.size != ysize )
    {
      gridMark4Update(gridID);
      gridptr->y.size = ysize;
    }

  if ( gridType != GRID_UNSTRUCTURED && gridType != GRID_PROJECTION )
    {
      long axisproduct = gridptr->x.size*gridptr->y.size;
      if ( axisproduct > 0 && axisproduct != gridSize )
        Error("Inconsistent grid declaration! (xsize=%d ysize=%d gridsize=%d)",
              gridptr->x.size, gridptr->y.size, gridSize);
    }
}

/*
@Function  gridInqYsize
@Title     Get the number of values of a Y-axis

@Prototype int gridInqYsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqYsize} returns the number of values of a Y-axis.

@Result
@func{gridInqYsize} returns the number of values of a Y-axis.

@EndFunction
*/
int gridInqYsize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->y.size;
}

/*
@Function  gridDefNP
@Title     Define the number of parallels between a pole and the equator

@Prototype void gridDefNP(int gridID, int np)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  np       Number of parallels between a pole and the equator.

@Description
The function @func{gridDefNP} defines the number of parallels between a pole and the equator
of a Gaussian grid.

@EndFunction
*/
void gridDefNP(int gridID, int np)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->np != np )
    {
      gridMark4Update(gridID);
      gridptr->np = np;
    }
}

/*
@Function  gridInqNP
@Title     Get the number of parallels between a pole and the equator

@Prototype int gridInqNP(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqNP} returns the number of parallels between a pole and the equator
of a Gaussian grid.

@Result
@func{gridInqNP} returns the number of parallels between a pole and the equator.

@EndFunction
*/
int gridInqNP(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->np;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefRowlon(int gridID, int nrowlon, const int rowlon[])
{
  grid_t *gridptr = grid_to_pointer(gridID);

  gridptr->rowlon = (int *) Malloc((size_t)nrowlon * sizeof(int));
  gridptr->nrowlon = nrowlon;
  memcpy(gridptr->rowlon, rowlon, (size_t)nrowlon * sizeof(int));
  gridMark4Update(gridID);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridInqRowlon(int gridID, int *rowlon)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->rowlon == 0 )  Error("undefined pointer!");

  memcpy(rowlon, gridptr->rowlon, (size_t)gridptr->nrowlon * sizeof(int));
}

static int
gridInqMaskSerialGeneric(grid_t *gridptr, mask_t **internalMask,
                        int *restrict mask)
{
  long size = gridptr->size;

  if ( CDI_Debug && size == 0 )
    Warning("Size undefined for gridID = %d", gridptr->self);

  const mask_t *restrict mask_src = *internalMask;
  if ( mask_src )
    {
      if (mask && size > 0)
        for (size_t i = 0; i < (size_t)size; ++i)
          mask[i] = (int)mask_src[i];
    }
  else
    size = 0;

  return (int)size;
}

static int
gridInqMaskSerial(grid_t *gridptr, int *mask)
{
  return gridInqMaskSerialGeneric(gridptr, &gridptr->mask, mask);
}


int gridInqMask(int gridID, int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqMask(gridptr, mask);
}

static void
gridDefMaskSerial(grid_t *gridptr, const int *mask)
{
  long size = gridptr->size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridptr->self);

  if ( mask == NULL )
    {
      if ( gridptr->mask )
	{
	  Free(gridptr->mask);
	  gridptr->mask = NULL;
	}
    }
  else
    {
      if ( gridptr->mask == NULL )
	gridptr->mask = (mask_t *) Malloc((size_t)size*sizeof(mask_t));
      else if ( CDI_Debug )
	Warning("grid mask already defined!");

      for (long i = 0; i < size; ++i )
	gridptr->mask[i] = (mask_t)(mask[i] != 0);
    }
}

void gridDefMask(int gridID, const int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defMask(gridptr, mask);
  gridMark4Update(gridID);
}

static int
gridInqMaskGMESerial(grid_t *gridptr, int *mask_gme)
{
  return gridInqMaskSerialGeneric(gridptr, &gridptr->mask_gme, mask_gme);
}

int gridInqMaskGME(int gridID, int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqMaskGME(gridptr, mask);
}

static void
gridDefMaskGMESerial(grid_t *gridptr, const int *mask)
{
  long size = gridptr->size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridptr->self);

  if ( gridptr->mask_gme == NULL )
    gridptr->mask_gme = (mask_t *) Malloc((size_t)size * sizeof (mask_t));
  else if ( CDI_Debug )
    Warning("mask already defined!");

  for (long i = 0; i < size; ++i)
    gridptr->mask_gme[i] = (mask_t)(mask[i] != 0);
}

void gridDefMaskGME(int gridID, const int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defMaskGME(gridptr, mask);
  gridMark4Update(gridID);
}

static
int gridInqXValsSerial(grid_t *gridptr, double *xvals)
{
  long size;
  if ( gridptr->type == GRID_CURVILINEAR || gridptr->type == GRID_UNSTRUCTURED )
    size = gridptr->size;
  else if ( gridptr->type == GRID_GAUSSIAN_REDUCED )
    size = 2;
  else
    size = gridptr->x.size;

  if ( CDI_Debug && size == 0 )
    Warning("size undefined for gridID = %d", gridptr->self);

  if ( gridptr->x.vals )
    {
      if ( size && xvals )
        {
          const double *gridptr_xvals = gridptr->vtable->inqXValsPtr(gridptr);
          memcpy(xvals, gridptr_xvals, (size_t)size * sizeof (double));
        }
    }
  else
    size = 0;

  return (int)size;
}

static
int gridInqXCvalsSerial(grid_t *gridptr, char **xcvals)
{
  if ( gridptr->type != GRID_CHARXY )
    Error("Function only valid for grid type 'GRID_CHARXY'.");

  int size = gridptr->x.size;
  int maxclength = 0;

  const char **gridptr_xcvals = gridptr->vtable->inqXCvalsPtr(gridptr);
  if ( gridptr_xcvals && size && xcvals )
    {
      maxclength = gridptr->x.clength;
      for ( int i = 0; i < size; i++ )
        memcpy(xcvals[i], gridptr_xcvals[i], (size_t)maxclength*sizeof(char));
    }

  return maxclength;
}

static
int gridInqXIscSerial(grid_t *gridptr)
{
  int clen = gridptr->x.clength;
  /*
  if ( gridptr->type != GRID_CHARXY )
    Error("Axis type is 'char' but grid is not type 'GRID_CHARXY'.");
  */
  return clen;
}

/*
@Function  gridInqXvals
@Title     Get all values of a X-axis

@Prototype int gridInqXvals(int gridID, double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  xvals    Pointer to the location into which the X-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXvals} returns all values of the X-axis.

@Result
Upon successful completion @func{gridInqXvals} returns the number of values and
the values are stored in @func{xvals}.
Otherwise, 0 is returned and @func{xvals} is empty.

@EndFunction
*/
int gridInqXvals(int gridID, double *xvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXVals(gridptr, xvals);
}


int gridInqXCvals(int gridID, char **xcvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXCvals(gridptr, xcvals);
}


int gridInqXIsc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXIsc(gridptr);
}

static
void gridDefXValsSerial(grid_t *gridptr, const double *xvals)
{
  int gridtype = gridptr->type;

  long size;
  if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    size = gridptr->size;
  else if ( gridtype == GRID_GAUSSIAN_REDUCED )
    size = 2;
  else
    size = gridptr->x.size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d", gridptr->self);

  if (gridptr->x.vals && CDI_Debug)
    Warning("values already defined!");
  gridptr->x.vals = (double *)Realloc(gridptr->x.vals,
                                      (size_t)size * sizeof(double));
  memcpy(gridptr->x.vals, xvals, (size_t)size * sizeof (double));
}

static
int gridInqYCvalsSerial(grid_t *gridptr, char **ycvals)
{
  if ( gridptr->type != GRID_CHARXY )
    Error("Function only valid for grid type 'GRID_CHARXY'.");

  int size = gridptr->y.size;
  int maxclength = 0;

  const char **gridptr_ycvals = gridptr->vtable->inqYCvalsPtr(gridptr);
  if ( gridptr_ycvals && size && ycvals )
    {
      maxclength = gridptr->y.clength;
      for ( int i = 0; i < size; i++ )
        memcpy(ycvals[i], gridptr_ycvals[i], (size_t)maxclength*sizeof(char));
    }

  return maxclength;
}

static
int gridInqYIscSerial(grid_t *gridptr)
{
  int clen = gridptr->y.clength;
  /*
  if ( gridptr->type != GRID_CHARXY )
    Error("Axis type is 'char' but grid is not type 'GRID_CHARXY'.");
  */
  return clen;
}

/*
@Function  gridDefXvals
@Title     Define the values of a X-axis

@Prototype void gridDefXvals(int gridID, const double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xvals    X-values of the grid.

@Description
The function @func{gridDefXvals} defines all values of the X-axis.

@EndFunction
*/
void gridDefXvals(int gridID, const double *xvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defXVals(gridptr, xvals);
  gridMark4Update(gridID);
}

static
int gridInqYValsSerial(grid_t *gridptr, double *yvals)
{
  int gridtype = gridptr->type;
  long size = (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
            ? gridptr->size : gridptr->y.size;

  if ( CDI_Debug && size == 0 )
    Warning("size undefined for gridID = %d!", gridptr->self);

  if ( gridptr->y.vals )
    {
      if ( size && yvals )
        {
          const double *gridptr_yvals = gridptr->vtable->inqYValsPtr(gridptr);
          memcpy(yvals, gridptr_yvals, (size_t)size * sizeof (double));
        }
    }
  else
    size = 0;

  return (int)size;
}

/*
@Function  gridInqYvals
@Title     Get all values of a Y-axis

@Prototype int gridInqYvals(int gridID, double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  yvals    Pointer to the location into which the Y-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYvals} returns all values of the Y-axis.

@Result
Upon successful completion @func{gridInqYvals} returns the number of values and
the values are stored in @func{yvals}.
Otherwise, 0 is returned and @func{yvals} is empty.

@EndFunction
*/
int gridInqYvals(int gridID, double *yvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYVals(gridptr, yvals);
}


int gridInqYCvals(int gridID, char **ycvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYCvals(gridptr, ycvals);
}


int gridInqYIsc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYIsc(gridptr);
}

static
void gridDefYValsSerial(grid_t *gridptr, const double *yvals)
{
  int gridtype = gridptr->type;
  long size = (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
            ? gridptr->size : gridptr->y.size;

  if ( size == 0 )
    Error("Size undefined for gridID = %d!", gridptr->self);

  if ( gridptr->y.vals && CDI_Debug )
    Warning("Values already defined!");

  gridptr->y.vals = (double *)Realloc(gridptr->y.vals, (size_t)size * sizeof (double));
  memcpy(gridptr->y.vals, yvals, (size_t)size * sizeof (double));
}


/*
@Function  gridDefYvals
@Title     Define the values of a Y-axis

@Prototype void gridDefYvals(int gridID, const double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  yvals    Y-values of the grid.

@Description
The function @func{gridDefYvals} defines all values of the Y-axis.

@EndFunction
*/
void gridDefYvals(int gridID, const double *yvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defYVals(gridptr, yvals);
  gridMark4Update(gridID);
}

static double
gridInqXValSerial(grid_t *gridptr, int index)
{
  double xval = gridptr->x.vals ? gridptr->x.vals[index] : 0;
  return xval;
}


double gridInqXval(int gridID, int index)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXVal(gridptr, index);
}

static double
gridInqYValSerial(grid_t *gridptr, int index)
{
  double yval = gridptr->y.vals ? gridptr->y.vals[index] : 0;
  return yval;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqYval(int gridID, int index)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYVal(gridptr, index);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqXinc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  double xinc = gridptr->x.inc;
  const double *restrict xvals = gridptr->vtable->inqXValsPtr(gridptr);

  if ( (! (fabs(xinc) > 0)) && xvals )
    {
      int xsize = gridptr->x.size;
      if ( xsize > 1 )
        {
          xinc = fabs(xvals[xsize-1] - xvals[0])/(xsize-1);
          for ( int i = 2; i < xsize; i++ )
            if ( fabs(fabs(xvals[i-1] - xvals[i]) - xinc) > 0.01*xinc )
              {
                xinc = 0;
                break;
              }

          gridptr->x.inc = xinc;
        }
    }

  return xinc;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double gridInqYinc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  double yinc = gridptr->y.inc;
  const double *yvals = gridptr->vtable->inqYValsPtr(gridptr);

  if ( (! (fabs(yinc) > 0)) && yvals )
    {
      int ysize = gridptr->y.size;
      if ( ysize > 1 )
        {
          yinc = yvals[1] - yvals[0];
          double abs_yinc = fabs(yinc);
          for ( size_t i = 2; i < (size_t)ysize; i++ )
            if ( fabs(fabs(yvals[i] - yvals[i-1]) - abs_yinc) > (0.01*abs_yinc) )
              {
                yinc = 0;
                break;
              }

          gridptr->y.inc = yinc;
        }
    }

  return yinc;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridInqParamRLL(int gridID, double *xpole, double *ypole, double *angle)
{
  *xpole = 0; *ypole = 0; *angle = 0;

  const char *projection = "rotated_latitude_longitude";
  char mapping[CDI_MAX_NAME]; mapping[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
  if ( mapping[0] && strcmp(mapping, projection) == 0 )
    {
      int atttype, attlen;
      char attname[CDI_MAX_NAME+1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      for ( int iatt = 0; iatt < natts; ++iatt )
        {
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if ( attlen != 1 ) continue;

          double attflt;
          if ( cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, &attflt) )
            {
              if      ( strcmp(attname, "grid_north_pole_longitude") == 0 ) *xpole = attflt;
              else if ( strcmp(attname, "grid_north_pole_latitude")  == 0 ) *ypole = attflt;
              else if ( strcmp(attname, "north_pole_grid_longitude") == 0 ) *angle = attflt;
            }
        }
    }
  else
    Warning("%s mapping parameter missing!", projection);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefParamRLL(int gridID, double xpole, double ypole, double angle)
{
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, CDI_MAX_NAME, "rotated_pole");

  const char *mapping = "rotated_latitude_longitude";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)(strlen(mapping)), mapping);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "grid_north_pole_longitude", CDI_DATATYPE_FLT64, 1, &xpole);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "grid_north_pole_latitude", CDI_DATATYPE_FLT64, 1, &ypole);
  if ( IS_NOT_EQUAL(angle, 0) ) cdiDefAttFlt(gridID, CDI_GLOBAL, "north_pole_grid_longitude", CDI_DATATYPE_FLT64, 1, &angle);

  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->projtype = CDI_PROJ_RLL;

  gridVerifyProj(gridID);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridInqParamGME(int gridID, int *nd, int *ni, int *ni2, int *ni3)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  *nd  = gridptr->gme.nd;
  *ni  = gridptr->gme.ni;
  *ni2 = gridptr->gme.ni2;
  *ni3 = gridptr->gme.ni3;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridDefParamGME(int gridID, int nd, int ni, int ni2, int ni3)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->gme.nd != nd )
    {
      gridptr->gme.nd  = nd;
      gridptr->gme.ni  = ni;
      gridptr->gme.ni2 = ni2;
      gridptr->gme.ni3 = ni3;
      gridMark4Update(gridID);
    }
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void gridChangeType(int gridID, int gridtype)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( CDI_Debug )
    Message("Changed grid type from %s to %s", gridNamePtr(gridptr->type), gridNamePtr(gridtype));

  if (gridptr->type != gridtype)
    {
      gridptr->type = gridtype;
      gridMark4Update(gridID);
    }
}

static
void grid_check_cyclic(grid_t *gridptr)
{
  gridptr->isCyclic = 0;
  enum { numVertices = 4 };
  size_t xsize = gridptr->x.size >= 0 ? (size_t)gridptr->x.size : 0,
         ysize = gridptr->y.size >= 0 ? (size_t)gridptr->y.size : 0;
  const double *xvals = gridptr->vtable->inqXValsPtr(gridptr),
               *yvals = gridptr->vtable->inqYValsPtr(gridptr),
    (*xbounds)[numVertices]
    = (const double (*)[numVertices])gridptr->vtable->inqXBoundsPtr(gridptr);

  if ( gridptr->type == GRID_GAUSSIAN || gridptr->type == GRID_LONLAT )
    {
      if ( xvals && xsize > 1 )
        {
          double xval1 = xvals[0];
          double xval2 = xvals[1];
          double xvaln = xvals[xsize-1];
          if ( xval2 < xval1 ) xval2 += 360;
          if ( xvaln < xval1 ) xvaln += 360;

          if ( IS_NOT_EQUAL(xval1, xvaln) )
            {
              double xinc = xval2 - xval1;
              if ( IS_EQUAL(xinc, 0) ) xinc = (xvaln - xval1)/(xsize-1);

              double x0 = xvaln + xinc - 360;

              if ( fabs(x0 - xval1) < 0.01*xinc ) gridptr->isCyclic = 1;
            }
        }
    }
  else if ( gridptr->type == GRID_CURVILINEAR )
    {
      bool lcheck = true;
      if ( yvals && xvals )
        {
          if ( (fabs(yvals[0] - yvals[xsize-1]) > fabs(yvals[0] - yvals[xsize*ysize-xsize])) &&
               (fabs(yvals[xsize*ysize-xsize] - yvals[xsize*ysize-1]) > fabs(yvals[xsize-1] - yvals[xsize*ysize-1])) )
            lcheck = false;
        }
      else lcheck = false;

      if ( lcheck && xvals && xsize > 1 )
        {
          size_t nc = 0;
          for ( size_t j = 0; j < ysize; ++j )
            {
              size_t i1 = j*xsize,
                     i2 = j*xsize+1,
                     in = j*xsize+(xsize-1);
              double val1 = xvals[i1],
                     val2 = xvals[i2],
                     valn = xvals[in];
              double xinc = fabs(val2-val1);

	      if ( val1 <    1 && valn > 300 ) val1 += 360;
	      if ( valn <    1 && val1 > 300 ) valn += 360;
	      if ( val1 < -179 && valn > 120 ) val1 += 360;
	      if ( valn < -179 && val1 > 120 ) valn += 360;
              if ( fabs(valn-val1) > 180 ) val1 += 360;

              double x0 = valn + copysign(xinc, val1 - valn);

              nc += fabs(x0-val1) < 0.5*xinc;
            }
          gridptr->isCyclic = nc > ysize/2;
        }

      if ( lcheck && xbounds && xsize > 1 )
	{
          bool isCyclic = true;
	  for ( size_t j = 0; j < ysize; ++j )
	    {
	      size_t i1 = j*xsize,
                     i2 = j*xsize+(xsize-1);
	      for (size_t k1 = 0; k1 < numVertices; ++k1 )
		{
		  double val1 = xbounds[i1][k1];
		  for (size_t k2 = 0; k2 < numVertices; ++k2 )
		    {
		      double val2 = xbounds[i2][k2];

		      if ( val1 <    1 && val2 > 300 ) val1 += 360;
		      if ( val2 <    1 && val1 > 300 ) val2 += 360;
		      if ( val1 < -179 && val2 > 120 ) val1 += 360;
		      if ( val2 < -179 && val1 > 120 ) val2 += 360;
                      if ( fabs(val2-val1) > 180 ) val1 += 360;

		      if ( fabs(val1-val2) < 0.001 )
                        goto foundCloseVertices;
		    }
		}
              /* all vertices more than 0.001 degrees apart */
              isCyclic = false;
              break;
              foundCloseVertices:
              ;
	    }
          gridptr->isCyclic = isCyclic;
	}
    }
}


int gridIsCircular(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->isCyclic == CDI_UNDEFID ) grid_check_cyclic(gridptr);

  return gridptr->isCyclic;
}

static
bool compareXYvals(grid_t *gridRef, grid_t *gridTest)
{
  bool differ = false;

  int xsizeTest = gridTest->x.size, ysizeTest = gridTest->y.size;
  if ( !differ && xsizeTest > 0 && xsizeTest == gridRef->vtable->inqXVals(gridRef, NULL) )
    {
      const double *restrict xvalsRef = gridRef->vtable->inqXValsPtr(gridRef),
        *restrict xvalsTest = gridTest->vtable->inqXValsPtr(gridTest);

      for ( size_t i = 0; i < (size_t)xsizeTest; ++i )
	if ( fabs(xvalsTest[i] - xvalsRef[i]) > 1.e-10 )
	  {
	    differ = true;
	    break;
	  }
    }

  if ( !differ && ysizeTest > 0 && ysizeTest == gridRef->vtable->inqYVals(gridRef, NULL) )
    {
      const double *restrict yvalsRef = gridRef->vtable->inqYValsPtr(gridRef),
        *restrict yvalsTest = gridTest->vtable->inqYValsPtr(gridTest);
      for ( size_t i = 0; i < (size_t)ysizeTest; ++i )
	if ( fabs(yvalsTest[i] - yvalsRef[i]) > 1.e-10 )
	  {
	    differ = true;
	    break;
	  }
    }

  return differ;
}

static
bool compareXYvals2(grid_t *gridRef, grid_t *gridTest)
{
  int gridsize = gridTest->size;
  bool differ = ((gridTest->x.vals == NULL) ^ (gridRef->x.vals == NULL))
             || ((gridTest->y.vals == NULL) ^ (gridRef->y.vals == NULL));

  typedef double (*inqVal)(grid_t *grid, int index);
  inqVal inqXValRef = gridRef->vtable->inqXVal,
         inqYValRef = gridRef->vtable->inqYVal,
         inqXValTest = gridTest->vtable->inqXVal,
         inqYValTest = gridTest->vtable->inqYVal;

  if ( !differ && gridTest->x.vals )
    differ = fabs(inqXValTest(gridTest, 0) - inqXValRef(gridRef, 0)) > 1.e-9
          || fabs(inqXValTest(gridTest, gridsize-1) - inqXValRef(gridRef, gridsize-1)) > 1.e-9;

  if ( !differ && gridTest->y.vals )
    differ = fabs(inqYValTest(gridTest, 0) - inqYValRef(gridRef, 0)) > 1.e-9
          || fabs(inqYValTest(gridTest, gridsize-1) - inqYValRef(gridRef, gridsize-1)) > 1.e-9;

  return differ;
}

static
bool gridCompare(int gridID, const grid_t *grid, bool coord_compare)
{
  bool differ = true;
  const grid_t *gridRef = grid_to_pointer(gridID);

  if ( grid->type == gridRef->type || grid->type == GRID_GENERIC )
    {
      if ( grid->size == gridRef->size )
	{
	  differ = false;
	  if ( grid->type == GRID_LONLAT || grid->type == GRID_PROJECTION )
	    {
	      /*
	      printf("gridID      %d\n", gridID);
	      printf("grid.xdef   %d\n", grid->x.flag);
	      printf("grid.ydef   %d\n", grid->y.flag);
	      printf("grid.xsize  %d\n", grid->x.size);
	      printf("grid.ysize  %d\n", grid->y.size);
	      printf("grid.xfirst %f\n", grid->x.first);
	      printf("grid.yfirst %f\n", grid->y.first);
	      printf("grid.xfirst %f\n", gridInqXval(gridID, 0));
	      printf("grid.yfirst %f\n", gridInqYval(gridID, 0));
	      printf("grid.xinc   %f\n", grid->x.inc);
	      printf("grid.yinc   %f\n", grid->y.inc);
	      printf("grid.xinc   %f\n", gridInqXinc(gridID));
	      printf("grid.yinc   %f\n", gridInqYinc(gridID));
	      */
	      if ( grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size )
		{
		  if ( grid->x.flag == 2 && grid->y.flag == 2 )
		    {
		      if ( ! (IS_EQUAL(grid->x.first, 0) && IS_EQUAL(grid->x.last, 0) && IS_EQUAL(grid->x.inc, 0)) &&
			   ! (IS_EQUAL(grid->y.first, 0) && IS_EQUAL(grid->y.last, 0) && IS_EQUAL(grid->y.inc, 0)) &&
			   IS_NOT_EQUAL(grid->x.first, grid->x.last) && IS_NOT_EQUAL(grid->y.first, grid->y.last) )
			{
			  if ( IS_NOT_EQUAL(grid->x.first, gridInqXval(gridID, 0)) ||
			       IS_NOT_EQUAL(grid->y.first, gridInqYval(gridID, 0)))
			    {
			      differ = true;
			    }
			  if ( !differ && fabs(grid->x.inc) > 0 &&
			       fabs(fabs(grid->x.inc) - fabs(gridRef->x.inc)) > fabs(grid->x.inc/1000))
			    {
			      differ = true;
			    }
			  if ( !differ && fabs(grid->y.inc) > 0 &&
			       fabs(fabs(grid->y.inc) - fabs(gridRef->y.inc)) > fabs(grid->y.inc/1000))
			    {
			      differ = true;
			    }
			}
		    }
		  else if ( grid->x.vals && grid->y.vals )
                    differ = gridRef->vtable->compareXYFull((grid_t *)gridRef, (grid_t *)grid);
		}
	      else
		differ = true;
	    }
	  else if ( grid->type == GRID_GENERIC )
	    {
	      if ( grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size )
		{
		  if ( grid->x.flag == 1 && grid->y.flag == 1
                       && grid->x.vals && grid->y.vals )
                    differ = gridRef->vtable->compareXYFull((grid_t *)gridRef, (grid_t *)grid);
		}
	      else if ( (grid->y.size == 0 || grid->y.size == 1) &&
			grid->x.size == gridRef->x.size*gridRef->y.size )
		{
		}
	      else
		differ = true;
	    }
	  else if ( grid->type == GRID_GAUSSIAN )
	    {
	      if ( grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size )
		{
		  if ( grid->x.flag == 2 && grid->y.flag == 2 )
		    {
		      if ( ! (IS_EQUAL(grid->x.first, 0) && IS_EQUAL(grid->x.last, 0) && IS_EQUAL(grid->x.inc, 0)) &&
			   ! (IS_EQUAL(grid->y.first, 0) && IS_EQUAL(grid->y.last, 0)) )
			if ( fabs(grid->x.first - gridInqXval(gridID, 0)) > 0.0015 ||
			     fabs(grid->y.first - gridInqYval(gridID, 0)) > 0.0015 ||
			     (fabs(grid->x.inc)>0 && fabs(fabs(grid->x.inc) - fabs(gridRef->x.inc)) > fabs(grid->x.inc/1000)) )
			  {
			    differ = true;
			  }
		    }
		  else if ( grid->x.vals && grid->y.vals )
                    differ = gridRef->vtable->compareXYFull((grid_t *)gridRef, (grid_t *)grid);
		}
	      else
		differ = true;
	    }
	  else if ( grid->type == GRID_CURVILINEAR )
	    {
	      /*
	      printf("gridID      %d\n", gridID);
	      printf("grid.xsize  %d\n", grid->x.size);
	      printf("grid.ysize  %d\n", grid->y.size);
	      printf("grid.xfirst %f\n", grid->x.vals[0]);
	      printf("grid.yfirst %f\n", grid->y.vals[0]);
	      printf("grid xfirst %f\n", gridInqXval(gridID, 0));
	      printf("grid yfirst %f\n", gridInqYval(gridID, 0));
	      printf("grid.xlast  %f\n", grid->x.vals[grid->size-1]);
	      printf("grid.ylast  %f\n", grid->y.vals[grid->size-1]);
	      printf("grid xlast  %f\n", gridInqXval(gridID, grid->size-1));
	      printf("grid ylast  %f\n", gridInqYval(gridID, grid->size-1));
	      printf("grid.nv     %d\n", grid->nvertex);
	      printf("grid nv     %d\n", gridInqNvertex(gridID));
	      */
	      if ( grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size )
                differ = gridRef->vtable->compareXYAO((grid_t *)gridRef, (grid_t *)grid);
	    }
	  else if ( grid->type == GRID_UNSTRUCTURED )
	    {
              if ( coord_compare )
                {
                  differ = grid->nvertex != gridRef->nvertex
                    || gridRef->vtable->compareXYAO((grid_t *)gridRef, (grid_t *)grid);
                }
              else
                {
                  /* FIXME: not octet 0 but octet 7 is guaranteed  non-zero for any non-NULL UUID */
                  differ = differ || (gridRef->uuid[0] && grid->uuid[0] && memcmp(gridRef->uuid, grid->uuid, CDI_UUID_SIZE) != 0);

                  if ( !differ &&
                       ((grid->x.vals == NULL) ^ (gridRef->x.vals == NULL)) &&
                       ((grid->y.vals == NULL) ^ (gridRef->y.vals == NULL)) )
                    {
                      int nvertexA, nvertexB, numberA, numberB;
                      differ = ( (nvertexA = grid->nvertex)
                                 && (nvertexB = gridRef->nvertex)
                                 && (nvertexA != nvertexB) )
                        || (numberA = grid->number, numberB = gridRef->number,
                            ( (numberA)
                              && numberB
                              && (numberA != numberB) )
                            || ( (numberA && numberB)
                                 && (grid->position) != (gridRef->position) ) );
                    }
                  else if ( !differ )
                    {
                      differ = grid->nvertex != gridRef->nvertex
                        || grid->number != gridRef->number
                        || (grid->number > 0 && grid->position != gridRef->position)
                        || gridRef->vtable->compareXYAO((grid_t *)gridRef, (grid_t *)grid);
                    }
                }
            }
	}
    }

  if ( (grid->scanningMode != gridInqScanningMode(gridID)) || (grid->uvRelativeToGrid != gridInqUvRelativeToGrid(gridID)) )
    {
      // often grid definition may differ in UV-relativeToGrid
      differ = 1;
#ifdef HIRLAM_EXTENSIONS
      if ( cdiDebugExt>=200 )
        printf("gridCompare(gridID=%d): Differs: grid.scanningMode [%d] != gridInqScanningMode(gridID) [%d] or  grid.uvRelativeToGrid [%ld] != gridInqUvRelativeToGrid(gridID) [%d]\n",
               gridID, grid->scanningMode, gridInqScanningMode(gridID), grid->uvRelativeToGrid, gridInqUvRelativeToGrid(gridID) );
#endif // HIRLAM_EXTENSIONS
    }
  return differ;
}

/*
int gridIsEqual(int gridID1, int gridID2)
{
  const grid_t *grid2 = grid_to_pointer(gridID2);

  int grid_is_equal = gridCompare(gridID1, grid2, true) == false;

  return grid_is_equal;
}
*/

int gridCompareP(void *gridptr1, void *gridptr2)
{
  grid_t *g1 = ( grid_t * ) gridptr1;
  grid_t *g2 = ( grid_t * ) gridptr2;
  enum { equal = 0,
         differ = -1 };
  int i, size;

  xassert ( g1 );
  xassert ( g2 );

  if ( g1->type          != g2->type         ) return differ;
  if ( g1->datatype      != g2->datatype     ) return differ;
  if ( g1->isCyclic      != g2->isCyclic     ) return differ;
  if ( g1->x.flag        != g2->x.flag       ) return differ;
  if ( g1->y.flag        != g2->y.flag       ) return differ;
  if ( g1->gme.nd        != g2->gme.nd       ) return differ;
  if ( g1->gme.ni        != g2->gme.ni       ) return differ;
  if ( g1->gme.ni2       != g2->gme.ni2      ) return differ;
  if ( g1->gme.ni3       != g2->gme.ni3      ) return differ;
  if ( g1->number        != g2->number       ) return differ;
  if ( g1->position      != g2->position     ) return differ;
  if ( g1->trunc         != g2->trunc        ) return differ;
  if ( g1->nvertex       != g2->nvertex      ) return differ;
  if ( g1->nrowlon       != g2->nrowlon      ) return differ;
  if ( g1->size          != g2->size         ) return differ;
  if ( g1->x.size        != g2->x.size       ) return differ;
  if ( g1->y.size        != g2->y.size       ) return differ;
  if ( g1->lcomplex      != g2->lcomplex     ) return differ;

  if ( IS_NOT_EQUAL(g1->x.first       , g2->x.first)       ) return differ;
  if ( IS_NOT_EQUAL(g1->y.first	      , g2->y.first)       ) return differ;
  if ( IS_NOT_EQUAL(g1->x.last        , g2->x.last)        ) return differ;
  if ( IS_NOT_EQUAL(g1->y.last        , g2->y.last)        ) return differ;
  if ( IS_NOT_EQUAL(g1->x.inc	      , g2->x.inc)         ) return differ;
  if ( IS_NOT_EQUAL(g1->y.inc	      , g2->y.inc)         ) return differ;
  if ( IS_NOT_EQUAL(g1->uvRelativeToGrid     , g2->uvRelativeToGrid)     ) return differ;
  if ( IS_NOT_EQUAL(g1->scanningMode         , g2->scanningMode)         ) return differ;

  const double *restrict g1_xvals = g1->vtable->inqXValsPtr(g1),
               *restrict g2_xvals = g2->vtable->inqXValsPtr(g2);
  if ( g1_xvals )
    {
      if ( g1->type == GRID_UNSTRUCTURED || g1->type == GRID_CURVILINEAR )
        size = g1->size;
      else
        size = g1->x.size;
      xassert ( size );

      if ( !g2_xvals ) return differ;

      for ( i = 0; i < size; i++ )
        if ( IS_NOT_EQUAL(g1_xvals[i], g2_xvals[i]) ) return differ;
    }
  else if ( g2_xvals )
    return differ;

  const double *restrict g1_yvals = g1->vtable->inqYValsPtr(g1),
               *restrict g2_yvals = g2->vtable->inqYValsPtr(g2);
  if ( g1_yvals )
    {
      if ( g1->type == GRID_UNSTRUCTURED || g1->type == GRID_CURVILINEAR )
	size = g1->size;
      else
	size = g1->y.size;
      xassert ( size );

      if ( !g2_yvals ) return differ;

      for ( i = 0; i < size; i++ )
        if ( IS_NOT_EQUAL(g1_yvals[i], g2_yvals[i]) ) return differ;
    }
  else if ( g2_yvals )
    return differ;

  const double *restrict g1_area = g1->vtable->inqAreaPtr(g1),
               *restrict g2_area = g2->vtable->inqAreaPtr(g2);
  if ( g1_area )
    {
      xassert ( g1->size );

      if ( !g2_area ) return differ;

      for ( i = 0; i < g1->size; i++ )
	if ( IS_NOT_EQUAL(g1_area[i], g2_area[i]) ) return differ;
    }
  else if ( g2_area )
    return differ;

  {
    const double *restrict g1_xbounds, *restrict g2_xbounds;
    if ( (g1_xbounds = g1->vtable->inqXBoundsPtr(g1)) )
      {
        xassert ( g1->nvertex );
        if ( g1->type == GRID_CURVILINEAR || g1->type == GRID_UNSTRUCTURED )
          size = g1->nvertex * g1->size;
        else
          size = g1->nvertex * g1->x.size;
        xassert ( size );

        if ( !(g2_xbounds = g2->vtable->inqXBoundsPtr(g2)) ) return differ;

        for ( i = 0; i < size; i++ )
          if ( IS_NOT_EQUAL(g1_xbounds[i], g2_xbounds[i]) ) return differ;
      }
    else if ( g2->vtable->inqXBoundsPtr(g2) )
      return differ;
  }

  {
    const double *restrict g1_ybounds, *restrict g2_ybounds;
    if ( (g1_ybounds = g1->vtable->inqYBoundsPtr(g1)) )
      {
        xassert ( g1->nvertex );
        if ( g1->type == GRID_CURVILINEAR || g1->type == GRID_UNSTRUCTURED )
          size = g1->nvertex * g1->size;
        else
          size = g1->nvertex * g1->y.size;
        xassert ( size );

        if ( ! (g2_ybounds = g2->vtable->inqYBoundsPtr(g2)) ) return differ;

        for ( i = 0; i < size; i++ )
          if ( IS_NOT_EQUAL(g1->y.bounds[i], g2->y.bounds[i]) ) return differ;
      }
    else if ( g2->vtable->inqYBoundsPtr(g2) )
      return differ;
  }

  if (strcmp(g1->x.name, g2->x.name)) return differ;
  if (strcmp(g1->y.name, g2->y.name)) return differ;
  if (strcmp(g1->x.longname, g2->x.longname)) return differ;
  if (strcmp(g1->y.longname, g2->y.longname)) return differ;
  if (g1->x.stdname != g2->x.stdname) return differ;
  if (g1->y.stdname != g2->y.stdname) return differ;
  if (strcmp(g1->x.units, g2->x.units)) return differ;
  if (strcmp(g1->y.units, g2->y.units)) return differ;

  if (strcmp(g1->mapping, g2->mapping)) return differ;

  if ( g1->reference )
    {
      if ( !g2->reference ) return differ;
      if ( strcmp(g1->reference, g2->reference) ) return differ;
    }
  else if ( g2->reference )
    return differ;

  if ( g1->mask )
    {
      xassert ( g1->size );
      if ( !g2->mask ) return differ;
      if ( memcmp ( g1->mask, g2->mask, (size_t)g1->size * sizeof(mask_t)) ) return differ;
    }
  else if ( g2->mask )
    return differ;

  if ( g1->mask_gme )
    {
      xassert ( g1->size );
      if ( !g2->mask_gme ) return differ;
      if ( memcmp ( g1->mask_gme, g2->mask_gme, (size_t)g1->size * sizeof(mask_t)) ) return differ;
    }
  else if ( g2->mask_gme )
    return differ;

  if ( memcmp(g1->uuid, g2->uuid, CDI_UUID_SIZE) )
    return differ;

  return equal;
}

static
void gridComplete(grid_t *grid)
{
  int gridID = grid->self;
  gridDefDatatype(gridID, grid->datatype);

  int gridtype = grid->type;
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_UNSTRUCTURED:
    case GRID_CURVILINEAR:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CHARXY:
      {
	if ( grid->x.size > 0 ) gridDefXsize(gridID, grid->x.size);
	if ( grid->y.size > 0 ) gridDefYsize(gridID, grid->y.size);

        if ( gridtype == GRID_GAUSSIAN ) gridDefNP(gridID, grid->np);

	if ( grid->nvertex > 0 )
	  gridDefNvertex(gridID, grid->nvertex);

	if ( grid->x.flag == 2 )
	  {
            assert(gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR);
	    double *xvals = (double *) Malloc((size_t)grid->x.size * sizeof (double));
	    gridGenXvals(grid->x.size, grid->x.first, grid->x.last, grid->x.inc, xvals);
	    grid->x.vals = xvals;
	    // gridDefXinc(gridID, grid->x.inc);
	  }

	if ( grid->y.flag == 2 )
	  {
            assert(gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR);
	    double *yvals = (double *) Malloc((size_t)grid->y.size * sizeof (double));
	    gridGenYvals(gridtype, grid->y.size, grid->y.first, grid->y.last, grid->y.inc, yvals);
	    grid->y.vals = yvals;
	    // gridDefYinc(gridID, grid->y.inc);
	  }

	if ( grid->projtype == CDI_PROJ_RLL )
	  {
	    if ( grid->x.name[0] == 0 || grid->x.name[0] == 'x' ) strcpy(grid->x.name, "rlon");
	    if ( grid->y.name[0] == 0 || grid->y.name[0] == 'y' ) strcpy(grid->y.name, "rlat");
	    if ( grid->x.longname[0] == 0 ) strcpy(grid->x.longname, "longitude in rotated pole grid");
	    if ( grid->y.longname[0] == 0 ) strcpy(grid->y.longname, "latitude in rotated pole grid");
            grid->x.stdname = xystdname_tab[grid_xystdname_grid_latlon][0];
            grid->y.stdname = xystdname_tab[grid_xystdname_grid_latlon][1];
	    if ( grid->x.units[0] == 0 ) strcpy(grid->x.units, "degrees");
	    if ( grid->y.units[0] == 0 ) strcpy(grid->y.units, "degrees");
	  }

        if ( gridtype == GRID_UNSTRUCTURED )
          {
            int number = grid->number;
            int position = grid->position >= 0 ? grid->position : 0;
            if ( number > 0 ) gridDefNumber(gridID, number);
            gridDefPosition(gridID, position);
          }

	break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
	gridDefNP(gridID, grid->np);
	gridDefYsize(gridID, grid->y.size);
        if ( grid->x.flag == 2 )
          {
            double xvals[2] = { grid->x.first, grid->x.last };
            gridDefXvals(gridID, xvals);
          }

        if ( grid->y.flag == 2 )
	  {
	    double *yvals = (double *) Malloc((size_t)grid->y.size * sizeof (double));
	    gridGenYvals(gridtype, grid->y.size, grid->y.first, grid->y.last, grid->y.inc, yvals);
            grid->y.vals = yvals;
	    /*
	    gridDefYinc(gridID, grid->y.inc);
	    */
	  }
	break;
      }
    case GRID_SPECTRAL:
      {
        gridDefTrunc(gridID, grid->trunc);
        if ( grid->lcomplex ) gridDefComplexPacking(gridID, 1);
        break;
      }
    case GRID_FOURIER:
      {
	gridDefTrunc(gridID, grid->trunc);
	break;
      }
    case GRID_GME:
      {
        gridDefParamGME(gridID, grid->gme.nd, grid->gme.ni, grid->gme.ni2, grid->gme.ni3);
        break;
      }
      /*
    case GRID_GENERIC:
      {
        if ( grid->x.size > 0 && grid->y.size > 0 )
          {
            gridDefXsize(gridID, grid->x.size);
            gridDefYsize(gridID, grid->y.size);
            if ( grid->x.vals ) gridDefXvals(gridID, grid->x.vals);
            if ( grid->y.vals ) gridDefYvals(gridID, grid->y.vals);
          }
        break;
      }
      */
    case GRID_TRAJECTORY:
      {
        gridDefXsize(gridID, 1);
        gridDefYsize(gridID, 1);
        break;
      }
    default:
      {
	Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  grid->x.name[CDI_MAX_NAME - 1] = 0;
  grid->x.longname[CDI_MAX_NAME - 1] = 0;
  grid->x.units[CDI_MAX_NAME - 1] = 0;
  grid->y.name[CDI_MAX_NAME - 1] = 0;
  grid->y.longname[CDI_MAX_NAME - 1] = 0;
  grid->y.units[CDI_MAX_NAME - 1] = 0;
}

#define GRID_STR_SERIALIZE(gridP) { gridP->x.dimname, gridP->y.dimname,  \
    gridP->vdimname, gridP->x.name, gridP->y.name,  \
    gridP->x.longname, gridP->y.longname, \
    gridP->x.units, gridP->y.units }

int gridGenerate(const grid_t *grid)
{
  int gridtype = grid->type;
  int gridID = gridCreate(gridtype, grid->size);
  grid_t *restrict gridptr = grid_to_pointer(gridID);
  gridptr->datatype = grid->datatype;
  gridptr->x.size = grid->x.size;
  gridptr->y.size = grid->y.size;
  gridptr->np = grid->np;
  gridptr->nvertex = grid->nvertex;
  gridptr->x.flag = grid->x.flag;
  int valdef_group1 = 0;
  static const int valdef_group1_tab[] = {
    GRID_LONLAT, GRID_GAUSSIAN, GRID_UNSTRUCTURED, GRID_CURVILINEAR,
    GRID_GENERIC, GRID_PROJECTION
  };
  for ( size_t i = 0; i < sizeof (valdef_group1_tab) / sizeof (valdef_group1_tab[0]); ++i)
    valdef_group1 |= (gridtype == valdef_group1_tab[i]);
  if ( valdef_group1 && grid->x.flag == 1 )
    {
      gridDefXvals(gridID, grid->x.vals);
      if ( grid->x.bounds )
        gridDefXbounds(gridID, grid->x.bounds);
    }
  gridptr->x.first = grid->x.first;
  gridptr->x.last = grid->x.last;
  gridptr->x.inc = grid->x.inc;
  gridptr->y.flag = grid->y.flag;
  if ( (valdef_group1 || gridtype == GRID_GAUSSIAN_REDUCED) && grid->y.flag == 1)
    {
      gridDefYvals(gridID, grid->y.vals);
      if ( grid->y.bounds )
        gridDefYbounds(gridID, grid->y.bounds);
    }
  gridptr->y.first = grid->y.first;
  gridptr->y.last = grid->y.last;
  gridptr->y.inc = grid->y.inc;
  if ( valdef_group1 && grid->area)
    gridDefArea(gridID, grid->area);
  gridptr->number = grid->number;
  gridptr->position = grid->position;
  gridptr->uvRelativeToGrid       = grid->uvRelativeToGrid;
  gridptr->scanningMode           = grid->scanningMode;
  gridptr->iScansNegatively       = grid->iScansNegatively;
  gridptr->jScansPositively       = grid->jScansPositively;
  gridptr->jPointsAreConsecutive  = grid->jPointsAreConsecutive;
  memcpy(gridptr->uuid, grid->uuid, CDI_UUID_SIZE);
  if ( gridtype == GRID_UNSTRUCTURED && grid->reference )
    gridDefReference(gridID, grid->reference);
  if ( gridtype == GRID_PROJECTION )
    gridptr->name = strdup(grid->name);
  if ( gridtype == GRID_GAUSSIAN_REDUCED )
    gridDefRowlon(gridID, grid->y.size, grid->rowlon);
  gridptr->trunc = grid->trunc;
  gridptr->lcomplex = grid->lcomplex;
  gridptr->gme.nd = grid->gme.nd;
  gridptr->gme.ni = grid->gme.ni;
  gridptr->gme.ni2 = grid->gme.ni2;
  gridptr->gme.ni3 = grid->gme.ni3;
  const char *grid_str_tab[] = GRID_STR_SERIALIZE(grid);
  char *gridptr_str_tab[] = GRID_STR_SERIALIZE(gridptr);
  for (size_t i = 0; i < sizeof (grid_str_tab) / sizeof (grid_str_tab[0]); ++i)
    if ( grid_str_tab[i][0] )
      memcpy(gridptr_str_tab[i], grid_str_tab[i], CDI_MAX_NAME);
  gridComplete(gridptr);

  return gridID;
}

static void
grid_copy_base_array_fields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  size_t nrowlon = (size_t)gridptrOrig->nrowlon;
  size_t gridsize = (size_t)gridptrOrig->size;
  int gridtype = gridptrOrig->type;
  int irregular = gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED;
  if ( nrowlon )
    {
      gridptrDup->rowlon = (int*) Malloc(nrowlon * sizeof(int));
      memcpy(gridptrDup->rowlon, gridptrOrig->rowlon, nrowlon * sizeof(int));
    }

  if ( gridptrOrig->x.vals != NULL )
    {
      size_t size  = irregular ? gridsize : (size_t)gridptrOrig->x.size;

      gridptrDup->x.vals = (double *)Malloc(size * sizeof (double));
      memcpy(gridptrDup->x.vals, gridptrOrig->x.vals, size * sizeof (double));
    }

  if ( gridptrOrig->y.vals != NULL )
    {
      size_t size  = irregular ? gridsize : (size_t)gridptrOrig->y.size;

      gridptrDup->y.vals = (double *)Malloc(size * sizeof (double));
      memcpy(gridptrDup->y.vals, gridptrOrig->y.vals, size * sizeof (double));
    }

  if ( gridptrOrig->x.bounds != NULL )
    {
      size_t size  = (irregular ? gridsize : (size_t)gridptrOrig->x.size)
        * (size_t)gridptrOrig->nvertex;

      gridptrDup->x.bounds = (double *)Malloc(size * sizeof (double));
      memcpy(gridptrDup->x.bounds, gridptrOrig->x.bounds, size * sizeof (double));
    }

  if ( gridptrOrig->y.bounds != NULL )
    {
      size_t size = (irregular ? gridsize : (size_t)gridptrOrig->y.size)
        * (size_t)gridptrOrig->nvertex;

      gridptrDup->y.bounds = (double *)Malloc(size * sizeof (double));
      memcpy(gridptrDup->y.bounds, gridptrOrig->y.bounds, size * sizeof (double));
    }

  {
    const double *gridptrOrig_area
      = gridptrOrig->vtable->inqAreaPtr(gridptrOrig);
    if ( gridptrOrig_area != NULL )
      {
        size_t size = gridsize;

        gridptrDup->area = (double *)Malloc(size * sizeof (double));
        memcpy(gridptrDup->area, gridptrOrig_area, size * sizeof (double));
      }
  }

  if ( gridptrOrig->mask != NULL )
    {
      size_t size = gridsize;

      gridptrDup->mask = (mask_t *)Malloc(size * sizeof(mask_t));
      memcpy(gridptrDup->mask, gridptrOrig->mask, size * sizeof (mask_t));
    }

  if ( gridptrOrig->mask_gme != NULL )
    {
      size_t size = gridsize;

      gridptrDup->mask_gme = (mask_t *)Malloc(size * sizeof (mask_t));
      memcpy(gridptrDup->mask_gme, gridptrOrig->mask_gme, size * sizeof(mask_t));
    }
}


/*
@Function  gridDuplicate
@Title     Duplicate a horizontal Grid

@Prototype int gridDuplicate(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridDuplicate} duplicates a horizontal Grid.

@Result
@func{gridDuplicate} returns an identifier to the duplicated Grid.

@EndFunction
*/
int gridDuplicate(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  grid_t *gridptrnew = gridptr->vtable->copy(gridptr);
  int gridIDnew = reshPut(gridptrnew, &gridOps);
  gridptrnew->self = gridIDnew;
  return gridIDnew;
}


void gridCompress(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_UNSTRUCTURED )
    {
      if ( gridptr->mask_gme != NULL )
	{
          size_t gridsize = (size_t)gridInqSize(gridID);
	  size_t nv = (size_t)gridptr->nvertex;
          double *restrict area
            = (double *)gridptr->vtable->inqAreaPtr(gridptr),
            *restrict xvals = (double *)gridptr->vtable->inqXValsPtr((grid_t *)gridptr),
            *restrict yvals = (double *)gridptr->vtable->inqYValsPtr((grid_t *)gridptr),
            *restrict xbounds = (double *)gridptr->vtable->inqXBoundsPtr(gridptr),
            *restrict ybounds = (double *)gridptr->vtable->inqYBoundsPtr(gridptr);
          mask_t *restrict mask_gme = gridptr->mask_gme;
          size_t *restrict selection = (size_t *)Malloc(gridsize * sizeof (selection[0]));
          size_t nselect;
          {
            size_t j = 0;
            for (size_t i = 0; i < gridsize; i++ )
              selection[j] = i, j += (mask_gme[i] != 0);
            nselect = j;
          }
          selection = (size_t *)Realloc(selection, nselect * sizeof (selection[0]));
          if (xvals)
            for (size_t i = 0; i < nselect; i++ )
	      xvals[i] = xvals[selection[i]];
          if (yvals)
            for (size_t i = 0; i < nselect; i++ )
              yvals[i] = yvals[selection[i]];
          if (area)
            for (size_t i = 0; i < nselect; i++ )
              area[i] = area[selection[i]];
          if (xbounds)
            for (size_t i = 0; i < nselect; i++ )
              for (size_t iv = 0; iv < nv; iv++)
                xbounds[i * nv + iv] = xbounds[selection[i] * nv + iv];
          if (ybounds)
            for (size_t i = 0; i < nselect; i++ )
              for (size_t iv = 0; iv < nv; iv++)
                ybounds[i * nv + iv] = ybounds[selection[i] * nv + iv];
          Free(selection);

	  /* fprintf(stderr, "grid compress %d %d %d\n", i, j, gridsize); */
	  gridsize = nselect;
	  gridptr->size  = (int)gridsize;
	  gridptr->x.size = (int)gridsize;
	  gridptr->y.size = (int)gridsize;

          double **resizeP[] = { &gridptr->x.vals, &gridptr->y.vals,
                                 &gridptr->area,
                                 &gridptr->x.bounds, &gridptr->y.bounds };
          size_t newSize[] = { gridsize, gridsize, gridsize, nv*gridsize,
                               nv*gridsize };
          for ( size_t i = 0; i < sizeof (resizeP) / sizeof (resizeP[0]); ++i)
            if ( *(resizeP[i]) )
              *(resizeP[i]) = (double *)Realloc(*(resizeP[i]), newSize[i]*sizeof(double));

	  Free(gridptr->mask_gme);
	  gridptr->mask_gme = NULL;
          gridMark4Update(gridID);
	}
    }
  else
    Warning("Unsupported grid type: %s", gridNamePtr(gridtype));
}

static void
gridDefAreaSerial(grid_t *gridptr, const double *area)
{
  size_t size = (size_t)gridptr->size;

  if ( size == 0 )
    Error("size undefined for gridID = %d", gridptr->self);

  if ( gridptr->area == NULL )
    gridptr->area = (double *) Malloc(size*sizeof(double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(gridptr->area, area, size * sizeof(double));
}


void gridDefArea(int gridID, const double *area)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defArea(gridptr, area);
  gridMark4Update(gridID);
}

static void
gridInqAreaSerial(grid_t *gridptr, double *area)
{
  if (gridptr->area)
    memcpy(area, gridptr->area, (size_t)gridptr->size * sizeof (double));
}


void gridInqArea(int gridID, double *area)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->inqArea(gridptr, area);
}

static int
gridHasAreaBase(grid_t *gridptr)
{
  return gridptr->area != NULL;
}

int gridHasArea(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->hasArea(gridptr);
}


static const double *gridInqAreaPtrBase(grid_t *gridptr)
{
  return gridptr->area;
}

const double *gridInqAreaPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqAreaPtr(gridptr);
}


void gridDefNvertex(int gridID, int nvertex)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->nvertex != nvertex)
    {
      gridptr->nvertex = nvertex;
      gridMark4Update(gridID);
    }
}


int gridInqNvertex(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->nvertex;
}

static void
gridDefBoundsGeneric(grid_t *gridptr, const double *bounds, int regularSize,
                     double **field)
{
  int irregular = gridptr->type == GRID_CURVILINEAR
    || gridptr->type == GRID_UNSTRUCTURED;
  size_t nvertex = (size_t)gridptr->nvertex;
  if ( nvertex == 0 )
    {
      Warning("nvertex undefined for gridID = %d. Cannot define bounds!",
              gridptr->self);
      return;
    }
  size_t size = nvertex * (size_t)(irregular ? gridptr->size : regularSize);
  if ( size == 0 )
    Error("size undefined for gridID = %d", gridptr->self);

  if (*field == NULL)
    *field = (double *)Malloc(size * sizeof (double));
  else if ( CDI_Debug )
    Warning("values already defined!");

  memcpy(*field, bounds, size * sizeof (double));
}


static void
gridDefXBoundsSerial(grid_t *gridptr, const double *xbounds)
{
  gridDefBoundsGeneric(gridptr, xbounds, gridptr->x.size, &gridptr->x.bounds);
}

/*
@Function  gridDefXbounds
@Title     Define the bounds of a X-axis

@Prototype void gridDefXbounds(int gridID, const double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xbounds  X-bounds of the grid.

@Description
The function @func{gridDefXbounds} defines all bounds of the X-axis.

@EndFunction
*/
void gridDefXbounds(int gridID, const double *xbounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defXBounds(gridptr, xbounds);
  gridMark4Update(gridID);
}

static int
gridInqXBoundsSerial(grid_t *gridptr, double *xbounds)
{
  size_t nvertex = (size_t)gridptr->nvertex;

  int irregular = gridptr->type == GRID_CURVILINEAR
    || gridptr->type == GRID_UNSTRUCTURED;
  size_t size = nvertex * (size_t)(irregular ? gridptr->size : gridptr->x.size);

  const double *gridptr_xbounds = gridptr->vtable->inqXBoundsPtr(gridptr);
  if ( gridptr_xbounds )
    {
      if ( size && xbounds )
        memcpy(xbounds, gridptr_xbounds, size * sizeof (double));
    }
  else
    size = 0;

  return (int)size;
}

/*
@Function  gridInqXbounds
@Title     Get the bounds of a X-axis

@Prototype int gridInqXbounds(int gridID, double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  xbounds  Pointer to the location into which the X-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXbounds} returns the bounds of the X-axis.

@Result
Upon successful completion @func{gridInqXbounds} returns the number of bounds and
the bounds are stored in @func{xbounds}.
Otherwise, 0 is returned and @func{xbounds} is empty.

@EndFunction
*/
int gridInqXbounds(int gridID, double *xbounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXBounds(gridptr, xbounds);
}

static const double *
gridInqXBoundsPtrSerial(grid_t *gridptr)
{
  return gridptr->x.bounds;
}


const double *gridInqXboundsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXBoundsPtr(gridptr);
}

static void
gridDefYBoundsSerial(grid_t *gridptr, const double *ybounds)
{
  gridDefBoundsGeneric(gridptr, ybounds, gridptr->y.size, &gridptr->y.bounds);
}

/*
@Function  gridDefYbounds
@Title     Define the bounds of a Y-axis

@Prototype void gridDefYbounds(int gridID, const double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ybounds  Y-bounds of the grid.

@Description
The function @func{gridDefYbounds} defines all bounds of the Y-axis.

@EndFunction
*/
void gridDefYbounds(int gridID, const double *ybounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defYBounds(gridptr, ybounds);
  gridMark4Update(gridID);
}

static int
gridInqYBoundsSerial(grid_t *gridptr, double *ybounds)
{
  size_t nvertex = (size_t)gridptr->nvertex;

  int irregular = gridptr->type == GRID_CURVILINEAR
    || gridptr->type == GRID_UNSTRUCTURED;
  size_t size = nvertex * (size_t)(irregular ? gridptr->size : gridptr->y.size);

  const double *gridptr_ybounds = gridptr->vtable->inqYBoundsPtr(gridptr);
  if ( gridptr_ybounds )
    {
      if ( size && ybounds )
        memcpy(ybounds, gridptr_ybounds, size * sizeof (double));
    }
  else
    size = 0;

  return (int)size;
}


/*
@Function  gridInqYbounds
@Title     Get the bounds of a Y-axis

@Prototype int gridInqYbounds(int gridID, double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  ybounds  Pointer to the location into which the Y-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYbounds} returns the bounds of the Y-axis.

@Result
Upon successful completion @func{gridInqYbounds} returns the number of bounds and
the bounds are stored in @func{ybounds}.
Otherwise, 0 is returned and @func{ybounds} is empty.

@EndFunction
*/
int gridInqYbounds(int gridID, double *ybounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYBounds(gridptr, ybounds);
}

static const double *
gridInqYBoundsPtrSerial(grid_t *gridptr)
{
  return gridptr->y.bounds;
}


const double *gridInqYboundsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYBoundsPtr(gridptr);
}

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], size_t nbyte0,
                       size_t n, const double vals[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static void
printIntsPrefixAutoBrk(FILE *fp, const char prefix[], size_t nbyte0,
                       size_t n, const int vals[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%d ", vals[i]);
    }
  fputs("\n", fp);
}

static void
printBounds(FILE *fp, int dig, const char prefix[], size_t nbyte0,
            size_t n, size_t nvertex, const double bounds[])
{
  fputs(prefix, fp);
  if ( n > 0 )
    {
      for ( size_t iv = 0; iv < nvertex; iv++ )
        fprintf(fp, "%.*g ", dig, bounds[iv]);
      for ( size_t i = 1; i < (size_t)n; i++ )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          for ( size_t iv = 0; iv < nvertex; iv++ )
            fprintf(fp, "%.*g ", dig, bounds[i*nvertex+iv]);
        }
      fputs("\n", fp);
    }
}

static void
printMask(FILE *fp, const char prefix[], size_t nbyte0,
          size_t n, const int mask[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%d ", mask[i]);
    }
  fputs("\n", fp);
}

static inline
void *resizeBuffer(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static
void gridPrintAttributes(FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  int atttype, attlen;
  char attname[CDI_MAX_NAME+1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for ( int iatt = 0; iatt < natts; ++iatt )
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      if ( attlen == 0 ) continue;

      if ( atttype == CDI_DATATYPE_TXT )
        {
          size_t attSize = (size_t)(attlen+1)*sizeof(char);
          char *atttxt = (char *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          atttxt[attlen] = 0;
          fprintf(fp, "ATTR_TXT: %s = \"%s\"\n", attname, atttxt);
        }
      else if ( atttype == CDI_DATATYPE_INT8  || atttype == CDI_DATATYPE_UINT8  ||
                atttype == CDI_DATATYPE_INT16 || atttype == CDI_DATATYPE_UINT16 ||
                atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32 )
        {
          size_t attSize = (size_t)attlen*sizeof(int);
          int *attint = (int *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          if ( attlen == 1 )
            fprintf(fp, "ATTR_INT: %s =", attname);
          else
            fprintf(fp, "ATTR_INT_%d: %s =", attlen, attname);
          for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %d", attint[i]);
          fprintf(fp, "\n");
        }
      else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
        {
          size_t attSize = (size_t)attlen * sizeof(double);
          double *attflt = (double *)resizeBuffer(&attBuf, &attBufSize, attSize);
          int dig = (atttype == CDI_DATATYPE_FLT64) ? 15 : 7;
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          if ( attlen == 1 )
            fprintf(fp, "ATTR_FLT: %s =", attname);
          else
            fprintf(fp, "ATTR_FLT_%d: %s =", attlen, attname);
          for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %.*g", dig, attflt[i]);
          fprintf(fp, "\n");
        }
    }

  Free(attBuf);
}

static
void gridPrintKernel(int gridID, int opt, FILE *fp)
{
  size_t xdimLen, ydimLen;
  char attstr[CDI_MAX_NAME];
  char attstr2[CDI_MAX_NAME];
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  const char  **xcvals  = gridInqXCvalsPtr(gridID);
  const char  **ycvals  = gridInqYCvalsPtr(gridID);
  size_t nxvals = (size_t) gridInqXvals(gridID, NULL);
  size_t nyvals = (size_t) gridInqYvals(gridID, NULL);
  size_t nxbounds = (size_t) gridInqXbounds(gridID, NULL);
  size_t nybounds = (size_t) gridInqYbounds(gridID, NULL);

  int type     = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  int xsize    = gridInqXsize(gridID);
  int ysize    = gridInqYsize(gridID);
  int xstrlen  = gridInqXIsc(gridID);
  int ystrlen  = gridInqYIsc(gridID);
  int nvertex  = gridInqNvertex(gridID);
  int datatype = gridInqDatatype(gridID);

  int dig = (datatype == CDI_DATATYPE_FLT64) ? 15 : 7;

  fprintf(fp, "gridtype  = %s\n" "gridsize  = %d\n", gridNamePtr(type), gridsize);

  if ( type != GRID_GME )
    {
      if ( type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER )
        {
          if ( xsize > 0 ) fprintf(fp, "xsize     = %d\n", xsize);
          if ( ysize > 0 ) fprintf(fp, "ysize     = %d\n", ysize);
        }

      if ( nxvals > 0 || xcvals )
        {
          if ( xstrlen )  fprintf(fp, "xstringlen= %d\n", xstrlen);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "xdimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xlongname = %s\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xunits    = %s\n", attstr);
        }

      if ( nyvals > 0 || ycvals )
        {
          if ( ystrlen )  fprintf(fp, "ystringlen= %d\n", ystrlen);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "ydimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "ylongname = %s\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yunits    = %s\n", attstr);
        }

      if ( type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_VDIMNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] ) fprintf(fp, "vdimname  = %s\n", attstr);
        }
      if ( type == GRID_UNSTRUCTURED && nvertex > 0 ) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_CHARXY:
      {
        if ( type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED ) fprintf(fp, "np        = %d\n", gridInqNP(gridID));

	if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED )
	  {
	    xdimLen = (size_t)gridsize;
	    ydimLen = (size_t)gridsize;
	  }
        else if ( type == GRID_GAUSSIAN_REDUCED )
          {
	    xdimLen = 2;
	    ydimLen = (size_t)ysize;
          }
	else
	  {
	    xdimLen = (size_t)xsize;
	    ydimLen = (size_t)ysize;
	  }

	if ( type == GRID_UNSTRUCTURED )
          {
            int number = gridInqNumber(gridID);
            int position = gridInqPosition(gridID);
            // const unsigned char *d;
            if ( number > 0 )
              {
                fprintf(fp, "number    = %d\n", number);
                if ( position >= 0 ) fprintf(fp, "position  = %d\n", position);
              }
            /*
              gridInqUUID(gridID, uuidOfHGrid);
              d = (unsigned char *) &uuidOfHGrid;
              fprintf(fp, "uuid      = %02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x\n",
              d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7],
              d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
            */
            if ( gridInqReference(gridID, NULL) )
              {
                char reference_link[8192];
                gridInqReference(gridID, reference_link);
                fprintf(fp, "uri       = %s\n", reference_link);
              }
          }

	if ( nxvals > 0 )
	  {
	    double xfirst = 0.0, xinc = 0.0;

	    if ( type == GRID_LONLAT     || type == GRID_GAUSSIAN ||
		 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		xfirst = gridInqXval(gridID, 0);
		xinc   = gridInqXinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(xinc, 0) && opt )
	      {
                fprintf(fp, "xfirst    = %.*g\n"
                        "xinc      = %.*g\n", dig, xfirst, dig, xinc);
	      }
	    else
	      {
                double *xvals = (double*) Malloc(nxvals*sizeof(double));
                gridInqXvals(gridID, xvals);
                static const char prefix[] = "xvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, nxvals, xvals);
                Free(xvals);
	      }
	  }

        if ( xcvals )
          {
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, attstr);
            if ( attstr[0] )
              fprintf(fp, "x%ss = %.*s\n", attstr, xstrlen, xcvals[0]);
            else
              fprintf(fp, "xstrings  = %.*s\n", xstrlen, xcvals[0]);
            for ( int i = 1; i < xsize; i++ )
              fprintf(fp, "          = %.*s\n", xstrlen, xcvals[i]);
          }

	if ( nxbounds )
	  {
            double *xbounds = (double*) Malloc(nxbounds*sizeof(double));
            gridInqXbounds(gridID, xbounds);
            static const char prefix[] = "xbounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, xdimLen, (size_t)nvertex, xbounds);
            Free(xbounds);
	  }

	if ( nyvals > 0 )
	  {
	    double yfirst = 0.0, yinc = 0.0;

	    if ( type == GRID_LONLAT || type == GRID_GENERIC ||
                 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		yfirst = gridInqYval(gridID, 0);
		yinc   = gridInqYinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(yinc, 0) && opt )
	      {
	  	fprintf(fp, "yfirst    = %.*g\n"
                        "yinc      = %.*g\n", dig, yfirst, dig, yinc);
	      }
	    else
	      {
                double *yvals = (double*) Malloc(nyvals*sizeof(double));
                gridInqYvals(gridID, yvals);
                static const char prefix[] = "yvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, nyvals, yvals);
                Free(yvals);
	      }
	  }

        if ( ycvals )
          {
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, attstr);
            if ( attstr[0] )
              fprintf(fp, "x%ss = %.*s\n", attstr, ystrlen, ycvals[0]);
            else
              fprintf(fp, "ystrings  = %.*s\n", ystrlen, ycvals[0]);
            for ( int i = 1; i < ysize; i++ )
              fprintf(fp, "          = %.*s\n", ystrlen, ycvals[i]);
          }

	if ( nybounds )
	  {
            double *ybounds = (double*) Malloc(nybounds*sizeof(double));
            gridInqYbounds(gridID, ybounds);
            static const char prefix[] = "ybounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, ydimLen, (size_t)nvertex, ybounds);
            Free(ybounds);
	  }

	if ( gridHasArea(gridID) )
	  {
            double *area = (double*) Malloc((size_t)gridsize*sizeof(double));
            gridInqArea(gridID, area);
            static const char prefix[] = "area      = ";
            printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, (size_t)gridsize, area);
            Free(area);
	  }

        if ( type == GRID_GAUSSIAN_REDUCED )
          {
            static const char prefix[] = "rowlon    = ";
            int *rowlon = (int *)Malloc((size_t)ysize*sizeof(int));
            gridInqRowlon(gridID, rowlon);
            printIntsPrefixAutoBrk(fp, prefix, sizeof(prefix)-1,
                                   (size_t)(ysize > 0 ? ysize : 0), rowlon);
            Free(rowlon);
          }

        if ( type == GRID_PROJECTION ) gridPrintAttributes(fp, gridID);

	break;
      }
    case GRID_SPECTRAL:
      {
        fprintf(fp, "truncation = %d\n"
                "complexpacking = %d\n", gridInqTrunc(gridID), gridInqComplexPacking(gridID) );
        break;
      }
    case GRID_FOURIER:
      {
	fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
	break;
      }
    case GRID_GME:
      {
        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
        fprintf(fp, "ni        = %d\n", ni );
        break;
      }
   default:
      {
	fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
        break;
      }
    }

  gridInqUUID(gridID, uuidOfHGrid);
  if ( !cdiUUIDIsNull(uuidOfHGrid) )
    {
      char uuidOfHGridStr[37];
      cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
      if ( uuidOfHGridStr[0] != 0 && strlen(uuidOfHGridStr) == 36 )
        fprintf(fp, "uuid      = %s\n", uuidOfHGridStr);
    }

  if ( gridInqMask(gridID, NULL) )
    {
      int *mask = (gridsize>0) ? (int*) Malloc((size_t)gridsize*sizeof(int)) : NULL;
      gridInqMask(gridID, mask);
      static const char prefix[] = "mask      = ";
      printMask(fp, prefix, sizeof(prefix)-1,
                (size_t)(gridsize > 0 ? gridsize : 0), mask);
      if ( mask ) Free(mask);
    }
}


void gridPrint(int gridID, int opt)
{
  gridPrintKernel(gridID, opt, stdout);
}


void gridPrintP(void *voidptr, FILE *fp)
{
  grid_t *gridptr = (grid_t *) voidptr;
  int gridID = gridptr->self;

  xassert( gridptr );

  gridPrintKernel(gridID, 0, fp);

  fprintf(fp,
          "datatype  = %d\n"
          "nd        = %d\n"
          "ni        = %d\n"
          "ni2       = %d\n"
          "ni3       = %d\n"
          "number    = %d\n"
          "position  = %d\n"
          "trunc     = %d\n"
          "lcomplex  = %d\n"
          "nrowlon   = %d\n",
          gridptr->datatype, gridptr->gme.nd, gridptr->gme.ni, gridptr->gme.ni2,
          gridptr->gme.ni3, gridptr->number, gridptr->position, gridptr->trunc,
          gridptr->lcomplex, gridptr->nrowlon );

  if ( gridptr->rowlon )
    {
      static const char prefix[] = "rowlon    = ";
      printIntsPrefixAutoBrk(fp, prefix, sizeof(prefix)-1,
                             (size_t)(gridptr->nrowlon > 0
                                      ? gridptr->nrowlon : 0), gridptr->rowlon);
    }

  if ( gridInqMaskGME(gridID, NULL) )
    {
      int gridsize = gridptr->size;
      int *mask = (gridsize>0) ? (int*) Malloc((size_t)gridsize*sizeof(int)) : NULL;
      gridInqMaskGME(gridID, mask);
      static const char prefix[] = "mask_gme  = ";
      printMask(fp, prefix, sizeof(prefix)-1,
                (size_t)(gridptr->size > 0 ? gridptr->size : 0), mask);
      if ( mask ) Free(mask);
    }
}

static const double *gridInqXValsPtrSerial(grid_t *gridptr)
{
  return gridptr->x.vals;
}

static const char **gridInqXCvalsPtrSerial(grid_t *gridptr)
{
  return (const char **) gridptr->x.cvals;
}


const double *gridInqXvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXValsPtr(gridptr);
}


const char **gridInqXCvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXCvalsPtr(gridptr);
}

static const double *gridInqYValsPtrSerial(grid_t *gridptr)
{
  return gridptr->y.vals;
}

static const char **gridInqYCvalsPtrSerial(grid_t *gridptr)
{
  return (const char **) gridptr->y.cvals;
}

const double *gridInqYvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYValsPtr(gridptr);
}

const char **gridInqYCvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYCvalsPtr(gridptr);
}

/*
@Function  gridDefParamLCC
@Title     Define the parameter of a Lambert Conformal Conic grid

@Prototype void gridDefParamLCC(int gridID, double missval, double lon_0, double lat_0, double lat_1, double lat_2, double a, double rf, double xval_0, double yval_0, double x_0, double y_0)
@Parameter
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate}.
    @Item  missval   Missing value
    @Item  lon_0     The East longitude of the meridian which is parallel to the Y-axis.
    @Item  lat_0     Latitude of the projection origin
    @Item  lat_1     First latitude from the pole at which the secant cone cuts the sphere.
    @Item  lat_2     Second latitude at which the secant cone cuts the sphere.
    @Item  a         Earth radius in metres (optional).
    @Item  rf        Inverse flattening (1/f) (optional).
    @Item  xval_0    Longitude of the first grid point in degree (optional).
    @Item  yval_0    Latitude of the first grid point in degree (optional).
    @Item  x_0       False easting (optional).
    @Item  y_0       False northing (optional).

@Description
The function @func{gridDefParamLCC} defines the parameter of a Lambert Conformal Conic grid.

@EndFunction
*/
void gridDefParamLCC(int gridID, double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                     double a, double rf, double xval_0, double yval_0, double x_0, double y_0)
{
  (void)lat_0;
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, CDI_MAX_NAME, "Lambert_Conformal");

  const char *mapname = "lambert_conformal_conic";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapname);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)(strlen(mapname)), mapname);
  int nlats = 0;
  double lats[2];
  lats[nlats++] = lat_1;
  if ( IS_NOT_EQUAL(lat_1, lat_2) ) lats[nlats++] = lat_2;
  cdiDefAttFlt(gridID, CDI_GLOBAL, "standard_parallel", CDI_DATATYPE_FLT64, nlats, lats);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "longitude_of_central_meridian", CDI_DATATYPE_FLT64, 1, &lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lat_2);
  if ( a > 0 ) cdiDefAttFlt(gridID, CDI_GLOBAL, "earth_radius", CDI_DATATYPE_FLT64, 1, &a);
  if ( rf > 0 ) cdiDefAttFlt(gridID, CDI_GLOBAL, "inverse_flattening", CDI_DATATYPE_FLT64, 1, &rf);
  if ( IS_NOT_EQUAL(x_0, missval) ) cdiDefAttFlt(gridID, CDI_GLOBAL, "false_easting", CDI_DATATYPE_FLT64, 1, &x_0);
  if ( IS_NOT_EQUAL(y_0, missval) ) cdiDefAttFlt(gridID, CDI_GLOBAL, "false_northing", CDI_DATATYPE_FLT64, 1, &y_0);
  if ( IS_NOT_EQUAL(xval_0, missval) ) cdiDefAttFlt(gridID, CDI_GLOBAL, "longitudeOfFirstGridPointInDegrees", CDI_DATATYPE_FLT64, 1, &xval_0);
  if ( IS_NOT_EQUAL(yval_0, missval) ) cdiDefAttFlt(gridID, CDI_GLOBAL, "latitudeOfFirstGridPointInDegrees", CDI_DATATYPE_FLT64, 1, &yval_0);

  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->projtype = CDI_PROJ_LCC;

  gridVerifyProj(gridID);
}

/*
@Function  gridInqParamLCC
@Title     Get the parameter of a Lambert Conformal Conic grid

@Prototype void gridInqParamLCC(int gridID, double missval, double *lon_0, double *lat_0, double *lat_1, double *lat_2, double *a, double *rf, double *xval_0, double *yval_0, double *x_0, double *y_0)
@Parameter
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  missval   Missing value
    @Item  lon_0     The East longitude of the meridian which is parallel to the Y-axis.
    @Item  lat_0     Latitude of the projection origin
    @Item  lat_1     First latitude from the pole at which the secant cone cuts the sphere.
    @Item  lat_2     Second latitude at which the secant cone cuts the sphere.
    @Item  a         Earth radius in metres (optional).
    @Item  rf        Inverse flattening (1/f) (optional).
    @Item  xval_0    Longitude of the first grid point in degree (optional).
    @Item  yval_0    Latitude of the first grid point in degree (optional).
    @Item  x_0       False easting (optional).
    @Item  y_0       False northing (optional).

@Description
The function @func{gridInqParamLCC} returns the parameter of a Lambert Conformal Conic grid.

@EndFunction
*/
int gridInqParamLCC(int gridID, double missval, double *lon_0, double *lat_0, double *lat_1, double *lat_2,
                    double *a, double *rf, double *xval_0, double *yval_0, double *x_0, double *y_0)
{
  *a = 0; *rf = 0;
  *lon_0 = missval; *lat_0 = missval, *lat_1 = missval, *lat_2 = missval;
  *xval_0 = missval; *yval_0 = missval; *x_0 = missval, *y_0 = missval;

  int status = -1;
  if ( gridInqType(gridID) != GRID_PROJECTION ) return status;

  status = -2;
  const char *projection = "lambert_conformal_conic";
  char mapname[CDI_MAX_NAME]; mapname[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapname);
  if ( mapname[0] && strcmp(mapname, projection) == 0 )
    {
      int atttype, attlen;
      char attname[CDI_MAX_NAME+1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      if ( natts ) status = 0;

      for ( int iatt = 0; iatt < natts; ++iatt )
        {
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if ( attlen > 2 ) continue;

          double attflt[2];
          if ( cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, attflt) )
            {
              if      ( strcmp(attname, "earth_radius") == 0 )                       *a      = attflt[0];
              else if ( strcmp(attname, "inverse_flattening") == 0 )                 *rf     = attflt[0];
              else if ( strcmp(attname, "longitude_of_central_meridian") == 0 )      *lon_0  = attflt[0];
              else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )      *lat_0  = attflt[0];
              else if ( strcmp(attname, "false_easting")  == 0 )                     *x_0    = attflt[0];
              else if ( strcmp(attname, "false_northing") == 0 )                     *y_0    = attflt[0];
              else if ( strcmp(attname, "longitudeOfFirstGridPointInDegrees") == 0 ) *xval_0 = attflt[0];
              else if ( strcmp(attname, "latitudeOfFirstGridPointInDegrees")  == 0 ) *yval_0 = attflt[0];
              else if ( strcmp(attname, "standard_parallel") == 0 )
                {
                  *lat_1 = attflt[0];
                  *lat_2 = (attlen == 2) ? attflt[1] : attflt[0];
                }
            }
        }
    }

  return status;
}


int gridVerifyGribParamLCC(double missval, double *lon_0, double *lat_0, double *lat_1, double *lat_2,
                           double *a, double *rf, double *xval_0, double *yval_0, double *x_0, double *y_0)
{
  static bool lwarn = true;

  if ( lwarn )
    {
      // lwarn = false;
      const char *projection = "lambert_conformal_conic";
      if ( IS_EQUAL(*lon_0, missval) ) { Warning("%s mapping parameter %s missing!", projection, "longitude_of_central_meridian"); }
      if ( IS_EQUAL(*lat_0, missval) ) { Warning("%s mapping parameter %s missing!", projection, "latitude_of_central_meridian"); }
      if ( IS_EQUAL(*lat_1, missval) ) { Warning("%s mapping parameter %s missing!", projection, "standard_parallel"); }
      if ( IS_NOT_EQUAL(*x_0, missval) && IS_NOT_EQUAL(*y_0, grid_missval) && (IS_EQUAL(*xval_0, missval) || IS_EQUAL(*yval_0, missval)) )
        {
          if ( proj_lcc_to_lonlat_func )
            {
              *xval_0 = -(*x_0); *yval_0 = -(*y_0);
              proj_lcc_to_lonlat_func(missval, *lon_0, *lat_0, *lat_1, *lat_2, *a, *rf, 0.0, 0.0, (size_t)1, xval_0, yval_0);
            }
          if ( IS_EQUAL(*xval_0, missval) || IS_EQUAL(*yval_0, missval) )
            Warning("%s mapping parameter %s missing!", projection, "longitudeOfFirstGridPointInDegrees and latitudeOfFirstGridPointInDegrees");
        }
    }

  return 0;
}


void gridDefComplexPacking(int gridID, int lcomplex)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->lcomplex != lcomplex)
    {
      gridptr->lcomplex = lcomplex != 0;
      gridMark4Update(gridID);
    }
}


int gridInqComplexPacking(int gridID)
{
  grid_t* gridptr = grid_to_pointer(gridID);

  return (int)gridptr->lcomplex;
}


void gridDefHasDims(int gridID, int hasdims)
{
  grid_t* gridptr = grid_to_pointer(gridID);

  if ( gridptr->hasdims != (hasdims != 0) )
    {
      gridptr->hasdims = hasdims != 0;
      gridMark4Update(gridID);
    }
}


int gridInqHasDims(int gridID)
{
  grid_t* gridptr = grid_to_pointer(gridID);

  return (int)gridptr->hasdims;
}

/*
@Function  gridDefNumber
@Title     Define the reference number for an unstructured grid

@Prototype void gridDefNumber(int gridID, const int number)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  number   Reference number for an unstructured grid.

@Description
The function @func{gridDefNumber} defines the reference number for an unstructured grid.

@EndFunction
*/
void gridDefNumber(int gridID, const int number)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->number != number )
    {
      gridptr->number = number;
      gridMark4Update(gridID);
    }
}

/*
@Function  gridInqNumber
@Title     Get the reference number to an unstructured grid

@Prototype int gridInqNumber(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqNumber} returns the reference number to an unstructured grid.

@Result
@func{gridInqNumber} returns the reference number to an unstructured grid.
@EndFunction
*/
int gridInqNumber(int gridID)
{
  grid_t* gridptr = grid_to_pointer(gridID);
  return gridptr->number;
}

/*
@Function  gridDefPosition
@Title     Define the position of grid in the reference file

@Prototype void gridDefPosition(int gridID, const int position)
@Parameter
    @Item  gridID     Grid ID, from a previous call to @fref{gridCreate}.
    @Item  position   Position of grid in the reference file.

@Description
The function @func{gridDefPosition} defines the position of grid in the reference file.

@EndFunction
*/
void gridDefPosition(int gridID, int position)
{
  grid_t* gridptr = grid_to_pointer(gridID);

  if ( gridptr->position != position )
    {
      gridptr->position = position;
      gridMark4Update(gridID);
    }
}

/*
@Function  gridInqPosition
@Title     Get the position of grid in the reference file

@Prototype int gridInqPosition(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqPosition} returns the position of grid in the reference file.

@Result
@func{gridInqPosition} returns the position of grid in the reference file.
@EndFunction
*/
int gridInqPosition(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->position;
}

/*
@Function  gridDefReference
@Title     Define the reference URI for an unstructured grid

@Prototype void gridDefReference(int gridID, const char *reference)
@Parameter
    @Item  gridID      Grid ID, from a previous call to @fref{gridCreate}.
    @Item  reference   Reference URI for an unstructured grid.

@Description
The function @func{gridDefReference} defines the reference URI for an unstructured grid.

@EndFunction
*/
void gridDefReference(int gridID, const char *reference)
{
  grid_t* gridptr = grid_to_pointer(gridID);

  if ( reference )
    {
      if ( gridptr->reference )
        {
          Free(gridptr->reference);
          gridptr->reference = NULL;
        }

      gridptr->reference = strdupx(reference);
      gridMark4Update(gridID);
    }
}

/*
@Function  gridInqReference
@Title     Get the reference URI to an unstructured grid

@Prototype char *gridInqReference(int gridID, char *reference)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqReference} returns the reference URI to an unstructured grid.

@Result
@func{gridInqReference} returns the reference URI to an unstructured grid.
@EndFunction
*/
int gridInqReference(int gridID, char *reference)
{
  size_t len = 0;
  grid_t* gridptr = grid_to_pointer(gridID);

  if ( gridptr->reference )
    {
      len = strlen(gridptr->reference);
      if ( reference )
        strcpy(reference, gridptr->reference);
    }

  return (int)len;
}

const char *gridInqReferencePtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->reference;
}

/*
@Function  gridDefUUID
@Title     Define the UUID for an unstructured grid

@Prototype void gridDefUUID(int gridID, const char *uuid)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  uuid     UUID for an unstructured grid.

@Description
The function @func{gridDefUUID} defines the UUID for an unstructured grid.

@EndFunction
*/
void gridDefUUID(int gridID, const unsigned char uuid[CDI_UUID_SIZE])
{
  grid_t* gridptr = grid_to_pointer(gridID);

  memcpy(gridptr->uuid, uuid, CDI_UUID_SIZE);
  gridMark4Update(gridID);
}

/*
@Function  gridInqUUID
@Title     Get the UUID to an unstructured grid

@Prototype void gridInqUUID(int gridID, char *uuid)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqUUID} returns the UUID to an unstructured grid.

@Result
@func{gridInqUUID} returns the UUID to an unstructured grid to the parameter uuid.
@EndFunction
*/
void gridInqUUID(int gridID, unsigned char uuid[CDI_UUID_SIZE])
{
  grid_t *gridptr = grid_to_pointer(gridID);

  memcpy(uuid, gridptr->uuid, CDI_UUID_SIZE);
}


void gridDefUvRelativeToGrid(int gridID, int uvRelativeToGrid)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->uvRelativeToGrid != uvRelativeToGrid )
    {
      gridMark4Update(gridID);
      gridptr->uvRelativeToGrid = (bool)uvRelativeToGrid;
    }
}


int gridInqUvRelativeToGrid(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->uvRelativeToGrid;
}


void gridDefScanningMode(int gridID, int mode)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if ( gridptr->scanningMode != mode )
    {
      gridMark4Update(gridID);
      gridptr->scanningMode = mode;
    }
}


int gridInqScanningMode(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int scanningModeTMP  = 128 * gridptr->iScansNegatively + 64 * gridptr->jScansPositively + 32 * gridptr->jPointsAreConsecutive;
  if ( scanningModeTMP != gridptr->scanningMode )
    Message("WARNING: scanningMode (%d) ! = (%d) 128 * iScansNegatively(%d) + 64 * jScansPositively(%d) + 32 * jPointsAreConsecutive(%d) ",
            gridptr->scanningMode, scanningModeTMP, gridptr->iScansNegatively,gridptr->jScansPositively,gridptr->jPointsAreConsecutive );

  return gridptr->scanningMode;
}


void cdiGridGetIndexList(unsigned ngrids, int * gridIndexList)
{
  reshGetResHListOfType(ngrids, gridIndexList, &gridOps);
}


static int
gridTxCode ()
{
  return GRID;
}

enum {
  GRID_PACK_INT_IDX_SELF,
  GRID_PACK_INT_IDX_TYPE,
  GRID_PACK_INT_IDX_DATATYPE,
  GRID_PACK_INT_IDX_IS_CYCLIC,
  GRID_PACK_INT_IDX_X_FLAG,
  GRID_PACK_INT_IDX_Y_FLAG,
  GRID_PACK_INT_IDX_GME_ND,
  GRID_PACK_INT_IDX_GME_NI,
  GRID_PACK_INT_IDX_GME_NI2,
  GRID_PACK_INT_IDX_GME_NI3,
  GRID_PACK_INT_IDX_NUMBER,
  GRID_PACK_INT_IDX_POSITION,
  GRID_PACK_INT_IDX_TRUNC,
  GRID_PACK_INT_IDX_NVERTEX,
  GRID_PACK_INT_IDX_NROWLON,
  GRID_PACK_INT_IDX_SIZE,
  GRID_PACK_INT_IDX_X_SIZE,
  GRID_PACK_INT_IDX_Y_SIZE,
  GRID_PACK_INT_IDX_LCOMPLEX,
  GRID_PACK_INT_IDX_MEMBERMASK,
  GRID_PACK_INT_IDX_XTSTDNNAME,
  GRID_PACK_INT_IDX_YTSTDNNAME,
  GRID_PACK_INT_IDX_UVRELATIVETOGRID,
  GRID_PACK_INT_IDX_ISCANSNEGATIVELY,
  GRID_PACK_INT_IDX_JSCANSPOSITIVELY,
  GRID_PACK_INT_IDX_JPOINTSARECONSECUTIVE,
  GRID_PACK_INT_IDX_SCANNINGMODE,
  gridNint
};

enum {
  GRID_PACK_DBL_IDX_X_FIRST,
  GRID_PACK_DBL_IDX_Y_FIRST,
  GRID_PACK_DBL_IDX_X_LAST,
  GRID_PACK_DBL_IDX_Y_LAST,
  GRID_PACK_DBL_IDX_X_INC,
  GRID_PACK_DBL_IDX_Y_INC,
  gridNdouble
};

enum {
       gridHasMaskFlag = 1 << 0,
       gridHasGMEMaskFlag = 1 << 1,
       gridHasXValsFlag = 1 << 2,
       gridHasYValsFlag = 1 << 3,
       gridHasAreaFlag = 1 << 4,
       gridHasXBoundsFlag = 1 << 5,
       gridHasYBoundsFlag = 1 << 6,
       gridHasReferenceFlag = 1 << 7,
       gridHasRowLonFlag = 1 << 8,
       gridHasUUIDFlag = 1 << 9,
};


static int gridGetComponentFlags(const grid_t * gridP)
{
  int flags = (gridHasMaskFlag & (int)((unsigned)(gridP->mask == NULL) - 1U))
    | (gridHasGMEMaskFlag & (int)((unsigned)(gridP->mask_gme == NULL) - 1U))
    | (gridHasXValsFlag
       & (int)((unsigned)(gridP->vtable->inqXValsPtr((grid_t *)gridP) == NULL) - 1U))
    | (gridHasYValsFlag
       & (int)((unsigned)(gridP->vtable->inqYValsPtr((grid_t *)gridP) == NULL) - 1U))
    | (gridHasAreaFlag
       & (int)((unsigned)(gridP->vtable->inqAreaPtr((grid_t *)gridP) == NULL)
               - 1U))
    | (gridHasXBoundsFlag & (int)((unsigned)(gridP->x.bounds == NULL) - 1U))
    | (gridHasYBoundsFlag & (int)((unsigned)(gridP->y.bounds == NULL) - 1U))
    | (gridHasReferenceFlag & (int)((unsigned)(gridP->reference == NULL) - 1U))
    | (gridHasRowLonFlag & (int)((unsigned)(gridP->rowlon == NULL) - 1U))
    | (gridHasUUIDFlag & (int)((unsigned)cdiUUIDIsNull(gridP->uuid) - 1U));
  return flags;
}

static int
gridGetPackSize(void * voidP, void *context)
{
  grid_t * gridP = ( grid_t * ) voidP;
  int packBuffSize = 0, count;

  packBuffSize += serializeGetSize(gridNint, CDI_DATATYPE_INT, context)
    + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (gridP->rowlon)
    {
      xassert(gridP->nrowlon);
      packBuffSize += serializeGetSize(gridP->nrowlon, CDI_DATATYPE_INT, context)
        + serializeGetSize( 1, CDI_DATATYPE_UINT32, context);
    }

  packBuffSize += serializeGetSize(gridNdouble, CDI_DATATYPE_FLT64, context);

  if (gridP->vtable->inqXValsPtr(gridP))
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR)
	count = gridP->size;
      else
	count = gridP->x.size;
      xassert(count);
      packBuffSize += serializeGetSize(count, CDI_DATATYPE_FLT64, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (gridP->vtable->inqYValsPtr(gridP))
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR)
	count = gridP->size;
      else
	count = gridP->y.size;
      xassert(count);
      packBuffSize += serializeGetSize(count, CDI_DATATYPE_FLT64, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (gridP->vtable->inqAreaPtr(gridP))
    {
      xassert(gridP->size);
      packBuffSize +=
        serializeGetSize(gridP->size, CDI_DATATYPE_FLT64, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (gridP->x.bounds)
    {
      xassert(gridP->nvertex);
      if (gridP->type == GRID_CURVILINEAR || gridP->type == GRID_UNSTRUCTURED)
	count = gridP->size;
      else
	count = gridP->x.size;
      xassert(count);
      packBuffSize
        += (serializeGetSize(gridP->nvertex * count, CDI_DATATYPE_FLT64, context)
            + serializeGetSize(1, CDI_DATATYPE_UINT32, context));
    }

  if (gridP->y.bounds)
    {
      xassert(gridP->nvertex);
      if (gridP->type == GRID_CURVILINEAR || gridP->type == GRID_UNSTRUCTURED)
	count = gridP->size;
      else
	count = gridP->y.size;
      xassert(count);
      packBuffSize
        += (serializeGetSize(gridP->nvertex * count, CDI_DATATYPE_FLT64, context)
            + serializeGetSize(1, CDI_DATATYPE_UINT32, context));
    }

  {
    const char *strTab[] = GRID_STR_SERIALIZE(gridP);
    int numStr = (int)(sizeof (strTab) / sizeof (strTab[0]));
    packBuffSize
      += serializeStrTabGetPackSize(strTab, numStr, context);
  }

  if (gridP->reference)
    {
      size_t len = strlen(gridP->reference);
      packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context)
        + serializeGetSize((int)len + 1, CDI_DATATYPE_TXT, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (gridP->mask)
    {
      xassert(gridP->size);
      packBuffSize
        += serializeGetSize(gridP->size, CDI_DATATYPE_UCHAR, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (gridP->mask_gme)
    {
      xassert(gridP->size);
      packBuffSize += serializeGetSize(gridP->size, CDI_DATATYPE_UCHAR, context)
        + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  if (!cdiUUIDIsNull(gridP->uuid))
    packBuffSize += serializeGetSize(CDI_UUID_SIZE, CDI_DATATYPE_UCHAR, context);

  return packBuffSize;
}

void
gridUnpack(char * unpackBuffer, int unpackBufferSize,
           int * unpackBufferPos, int originNamespace, void *context,
           int force_id)
{
  grid_t * gridP;
  uint32_t d;
  int memberMask, size;

  gridInit();

  {
    int intBuffer[gridNint];
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                    intBuffer, gridNint, CDI_DATATYPE_INT, context);
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                    &d, 1, CDI_DATATYPE_UINT32, context);

    xassert(cdiCheckSum(CDI_DATATYPE_INT, gridNint, intBuffer) == d);
    int targetID = namespaceAdaptKey(intBuffer[0], originNamespace);
    gridP = gridNewEntry(force_id?targetID:CDI_UNDEFID);

    xassert(!force_id || targetID == gridP->self);

    gridP->type          =   intBuffer[GRID_PACK_INT_IDX_TYPE];
    gridP->datatype      =   intBuffer[GRID_PACK_INT_IDX_DATATYPE];
    gridP->isCyclic      =   (signed char)intBuffer[GRID_PACK_INT_IDX_IS_CYCLIC];
    gridP->x.flag        =   (short)intBuffer[GRID_PACK_INT_IDX_X_FLAG];
    gridP->y.flag        =   (short)intBuffer[GRID_PACK_INT_IDX_Y_FLAG];
    gridP->gme.nd        =   intBuffer[GRID_PACK_INT_IDX_GME_ND];
    gridP->gme.ni        =   intBuffer[GRID_PACK_INT_IDX_GME_NI];
    gridP->gme.ni2       =   intBuffer[GRID_PACK_INT_IDX_GME_NI2];
    gridP->gme.ni3       =   intBuffer[GRID_PACK_INT_IDX_GME_NI3];
    gridP->number        =   intBuffer[GRID_PACK_INT_IDX_NUMBER];
    gridP->position      =   intBuffer[GRID_PACK_INT_IDX_POSITION];
    gridP->trunc         =   intBuffer[GRID_PACK_INT_IDX_TRUNC];
    gridP->nvertex       =   intBuffer[GRID_PACK_INT_IDX_NVERTEX];
    gridP->nrowlon       =   intBuffer[GRID_PACK_INT_IDX_NROWLON];
    gridP->size          =   intBuffer[GRID_PACK_INT_IDX_SIZE];
    gridP->x.size        =   intBuffer[GRID_PACK_INT_IDX_X_SIZE];
    gridP->y.size        =   intBuffer[GRID_PACK_INT_IDX_Y_SIZE];
    gridP->lcomplex      =   (bool)intBuffer[GRID_PACK_INT_IDX_LCOMPLEX];
    memberMask           =   intBuffer[GRID_PACK_INT_IDX_MEMBERMASK];
    gridP->x.stdname     =
      xystdname_tab[intBuffer[GRID_PACK_INT_IDX_XTSTDNNAME]][0];
    gridP->y.stdname     =
      xystdname_tab[intBuffer[GRID_PACK_INT_IDX_YTSTDNNAME]][1];
    gridP->uvRelativeToGrid         =   intBuffer[GRID_PACK_INT_IDX_UVRELATIVETOGRID];
    gridP->iScansNegatively         =   (bool)intBuffer[GRID_PACK_INT_IDX_ISCANSNEGATIVELY];
    gridP->jScansPositively         =   (bool)intBuffer[GRID_PACK_INT_IDX_JSCANSPOSITIVELY];
    gridP->jPointsAreConsecutive    =   (bool)intBuffer[GRID_PACK_INT_IDX_JPOINTSARECONSECUTIVE];
    gridP->scanningMode             =   intBuffer[GRID_PACK_INT_IDX_SCANNINGMODE];
  }

  if (memberMask & gridHasRowLonFlag)
    {
      xassert(gridP->nrowlon);
      gridP->rowlon = (int *) Malloc((size_t)gridP->nrowlon * sizeof (int));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->rowlon, gridP->nrowlon , CDI_DATATYPE_INT, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_INT, gridP->nrowlon, gridP->rowlon) == d);
    }

  {
    double doubleBuffer[gridNdouble];
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                    doubleBuffer, gridNdouble, CDI_DATATYPE_FLT64, context);
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                    &d, 1, CDI_DATATYPE_UINT32, context);
    xassert(d == cdiCheckSum(CDI_DATATYPE_FLT, gridNdouble, doubleBuffer));

    gridP->x.first = doubleBuffer[GRID_PACK_DBL_IDX_X_FIRST];
    gridP->y.first = doubleBuffer[GRID_PACK_DBL_IDX_Y_FIRST];
    gridP->x.last = doubleBuffer[GRID_PACK_DBL_IDX_X_LAST];
    gridP->y.last = doubleBuffer[GRID_PACK_DBL_IDX_Y_LAST];
    gridP->x.inc = doubleBuffer[GRID_PACK_DBL_IDX_X_INC];
    gridP->y.inc = doubleBuffer[GRID_PACK_DBL_IDX_Y_INC];
  }

  int irregular = gridP->type == GRID_UNSTRUCTURED
    || gridP->type == GRID_CURVILINEAR;
  if (memberMask & gridHasXValsFlag)
    {
      size = irregular ? gridP->size : gridP->x.size;

      gridP->x.vals = (double *) Malloc((size_t)size * sizeof (double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->x.vals, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.vals) == d );
    }

  if (memberMask & gridHasYValsFlag)
    {
      size = irregular ? gridP->size : gridP->y.size;

      gridP->y.vals = (double *) Malloc((size_t)size * sizeof (double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->y.vals, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.vals) == d);
    }

  if (memberMask & gridHasAreaFlag)
    {
      size = gridP->size;
      xassert(size);
      gridP->area = (double *) Malloc((size_t)size * sizeof (double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->area, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->area) == d);
    }

  if (memberMask & gridHasXBoundsFlag)
    {
      size = gridP->nvertex * (irregular ? gridP->size : gridP->x.size);
      xassert(size);

      gridP->x.bounds = (double *) Malloc((size_t)size * sizeof (double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->x.bounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.bounds) == d);
    }

  if (memberMask & gridHasYBoundsFlag)
    {
      size = gridP->nvertex * (irregular ? gridP->size : gridP->y.size);
      xassert(size);

      gridP->y.bounds = (double *) Malloc((size_t)size * sizeof (double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
			  gridP->y.bounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.bounds) == d);
    }

  {
    char *strTab[] = GRID_STR_SERIALIZE(gridP);
    int numStr = sizeof (strTab) / sizeof (strTab[0]);
    serializeStrTabUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                          strTab, numStr, context);
  }

  if (memberMask & gridHasReferenceFlag)
    {
      int referenceSize;
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &referenceSize, 1, CDI_DATATYPE_INT, context);
      gridP->reference = (char *) Malloc((size_t)referenceSize);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->reference, referenceSize, CDI_DATATYPE_TXT, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_TXT, referenceSize, gridP->reference) == d);
    }

  if (memberMask & gridHasMaskFlag)
    {
      xassert((size = gridP->size));
      gridP->mask = (mask_t *) Malloc((size_t)size * sizeof (mask_t));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->mask, gridP->size, CDI_DATATYPE_UCHAR, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_UCHAR, gridP->size, gridP->mask) == d);
    }

  if (memberMask & gridHasGMEMaskFlag)
    {
      xassert((size = gridP->size));
      gridP->mask_gme = (mask_t *) Malloc((size_t)size * sizeof (mask_t));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->mask_gme, gridP->size, CDI_DATATYPE_UCHAR, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_UCHAR, gridP->size, gridP->mask_gme) == d);
    }
  if (memberMask & gridHasUUIDFlag)
    {
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      gridP->uuid, CDI_UUID_SIZE, CDI_DATATYPE_UCHAR, context);
    }

  reshSetStatus(gridP->self, &gridOps,
                reshGetStatus(gridP->self, &gridOps) & ~RESH_SYNC_BIT);
}


static void
gridPack(void * voidP, void * packBuffer, int packBufferSize,
         int * packBufferPos, void *context)
{
  grid_t   * gridP = ( grid_t * )   voidP;
  int size;
  uint32_t d;
  int memberMask;

  {
    int intBuffer[gridNint];

    intBuffer[GRID_PACK_INT_IDX_SELF]         = gridP->self;
    intBuffer[GRID_PACK_INT_IDX_TYPE]         = gridP->type;
    intBuffer[GRID_PACK_INT_IDX_DATATYPE]     = gridP->datatype;
    intBuffer[GRID_PACK_INT_IDX_IS_CYCLIC]    = gridP->isCyclic;
    intBuffer[GRID_PACK_INT_IDX_X_FLAG]       = gridP->x.flag;
    intBuffer[GRID_PACK_INT_IDX_Y_FLAG]       = gridP->y.flag;
    intBuffer[GRID_PACK_INT_IDX_GME_ND]       = gridP->gme.nd;
    intBuffer[GRID_PACK_INT_IDX_GME_NI]       = gridP->gme.ni;
    intBuffer[GRID_PACK_INT_IDX_GME_NI2]      = gridP->gme.ni2;
    intBuffer[GRID_PACK_INT_IDX_GME_NI3]      = gridP->gme.ni3;
    intBuffer[GRID_PACK_INT_IDX_NUMBER]       = gridP->number;
    intBuffer[GRID_PACK_INT_IDX_POSITION]     = gridP->position;
    intBuffer[GRID_PACK_INT_IDX_TRUNC]        = gridP->trunc;
    intBuffer[GRID_PACK_INT_IDX_NVERTEX]      = gridP->nvertex;
    intBuffer[GRID_PACK_INT_IDX_NROWLON]      = gridP->nrowlon;
    intBuffer[GRID_PACK_INT_IDX_SIZE]         = gridP->size;
    intBuffer[GRID_PACK_INT_IDX_X_SIZE]       = gridP->x.size;
    intBuffer[GRID_PACK_INT_IDX_Y_SIZE]       = gridP->y.size;
    intBuffer[GRID_PACK_INT_IDX_LCOMPLEX]     = gridP->lcomplex;
    intBuffer[GRID_PACK_INT_IDX_MEMBERMASK]   = memberMask
                                              = gridGetComponentFlags(gridP);
    intBuffer[GRID_PACK_INT_IDX_XTSTDNNAME]   =
      (int)((const char (*)[2][24])gridP->x.stdname - xystdname_tab);
    intBuffer[GRID_PACK_INT_IDX_YTSTDNNAME]   =
      (int)((const char (*)[2][24])gridP->y.stdname
            - (const char (*)[2][24])xystdname_tab[0][1]);

    intBuffer[GRID_PACK_INT_IDX_UVRELATIVETOGRID] = gridP->uvRelativeToGrid;
    intBuffer[GRID_PACK_INT_IDX_ISCANSNEGATIVELY] = gridP->iScansNegatively;
    intBuffer[GRID_PACK_INT_IDX_JSCANSPOSITIVELY] = gridP->jScansPositively;
    intBuffer[GRID_PACK_INT_IDX_JPOINTSARECONSECUTIVE] = gridP->jPointsAreConsecutive;
    intBuffer[GRID_PACK_INT_IDX_SCANNINGMODE] = gridP->scanningMode;

    serializePack(intBuffer, gridNint, CDI_DATATYPE_INT,
                  packBuffer, packBufferSize, packBufferPos, context);
    d = cdiCheckSum(CDI_DATATYPE_INT, gridNint, intBuffer);
    serializePack(&d, 1, CDI_DATATYPE_UINT32,
                  packBuffer, packBufferSize, packBufferPos, context);
  }

  if (memberMask & gridHasRowLonFlag)
    {
      size = gridP->nrowlon;
      xassert(size > 0);
      serializePack(gridP->rowlon, size, CDI_DATATYPE_INT,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_INT , size, gridP->rowlon);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  {
    double doubleBuffer[gridNdouble];

    doubleBuffer[GRID_PACK_DBL_IDX_X_FIRST]        = gridP->x.first;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_FIRST]        = gridP->y.first;
    doubleBuffer[GRID_PACK_DBL_IDX_X_LAST]         = gridP->x.last;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_LAST]         = gridP->y.last;
    doubleBuffer[GRID_PACK_DBL_IDX_X_INC]          = gridP->x.inc;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_INC]          = gridP->y.inc;

    serializePack(doubleBuffer, gridNdouble, CDI_DATATYPE_FLT64,
                  packBuffer, packBufferSize, packBufferPos, context);
    d = cdiCheckSum(CDI_DATATYPE_FLT, gridNdouble, doubleBuffer);
    serializePack(&d, 1, CDI_DATATYPE_UINT32,
                  packBuffer, packBufferSize, packBufferPos, context);
  }

  if (memberMask & gridHasXValsFlag)
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR)
	size = gridP->size;
      else
	size = gridP->x.size;
      xassert(size);

      const double *gridP_xvals = gridP->vtable->inqXValsPtr(gridP);
      serializePack(gridP_xvals, size, CDI_DATATYPE_FLT64,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP_xvals);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasYValsFlag)
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR )
	size = gridP->size;
      else
	size = gridP->y.size;
      xassert(size);
      const double *gridP_yvals = gridP->vtable->inqYValsPtr(gridP);
      serializePack(gridP_yvals, size, CDI_DATATYPE_FLT64,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP_yvals);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasAreaFlag)
    {
      xassert(gridP->size);

      serializePack(gridP->area, gridP->size, CDI_DATATYPE_FLT64,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, gridP->size, gridP->area);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasXBoundsFlag)
    {
      xassert ( gridP->nvertex );
      if (gridP->type == GRID_CURVILINEAR || gridP->type == GRID_UNSTRUCTURED)
	size = gridP->nvertex * gridP->size;
      else
	size = gridP->nvertex * gridP->x.size;
      xassert ( size );

      serializePack(gridP->x.bounds, size, CDI_DATATYPE_FLT64,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.bounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasYBoundsFlag)
    {
      xassert(gridP->nvertex);
      if (gridP->type == GRID_CURVILINEAR || gridP->type == GRID_UNSTRUCTURED)
	size = gridP->nvertex * gridP->size;
      else
	size = gridP->nvertex * gridP->y.size;
      xassert ( size );

      serializePack(gridP->y.bounds, size, CDI_DATATYPE_FLT64,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.bounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  {
    const char *strTab[] = GRID_STR_SERIALIZE(gridP);
    int numStr = sizeof (strTab) / sizeof (strTab[0]);
    serializeStrTabPack(strTab, numStr,
                        packBuffer, packBufferSize, packBufferPos, context);
  }

  if (memberMask & gridHasReferenceFlag)
    {
      size = (int)strlen(gridP->reference) + 1;
      serializePack(&size, 1, CDI_DATATYPE_INT,
                    packBuffer, packBufferSize, packBufferPos, context);
      serializePack(gridP->reference, size, CDI_DATATYPE_TXT,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_TXT, size, gridP->reference);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasMaskFlag)
    {
      xassert((size = gridP->size));
      serializePack(gridP->mask, size, CDI_DATATYPE_UCHAR,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_UCHAR, size, gridP->mask);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasGMEMaskFlag)
    {
      xassert((size = gridP->size));

      serializePack(gridP->mask_gme, size, CDI_DATATYPE_UCHAR,
                    packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_UCHAR, size, gridP->mask_gme);
      serializePack(&d, 1, CDI_DATATYPE_UINT32,
                    packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasUUIDFlag)
    serializePack(gridP->uuid, CDI_UUID_SIZE, CDI_DATATYPE_UCHAR,
                  packBuffer, packBufferSize, packBufferPos, context);
}

#undef GRID_STR_SERIALIZE


struct gridCompareSearchState
{
  int resIDValue;
  const grid_t *queryKey;
};

static enum cdiApplyRet
gridCompareSearch(int id, void *res, void *data)
{
  struct gridCompareSearchState *state = (struct gridCompareSearchState*)data;
  (void)res;
  if ( gridCompare(id, state->queryKey, true) == false )
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

/* Add grid (which must be Malloc'ed to vlist if not already found) */
struct addIfNewRes cdiVlistAddGridIfNew(int vlistID, grid_t *grid, int mode)
{
  /*
    mode: 0 search in vlist and grid table
          1 search in grid table only
          2 search in grid table only and don't store the grid in vlist
   */
  bool gridglobdefined = false;
  bool griddefined = false;
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  unsigned ngrids = (unsigned)vlistptr->ngrids;

  if ( mode == 0 )
    for ( unsigned index = 0; index < ngrids; index++ )
      {
	if ( (gridID = vlistptr->gridIDs[index]) != CDI_UNDEFID )
          {
            if ( gridCompare(gridID, grid, false) == false )
              {
                griddefined = true;
                break;
              }
          }
        else
          Error("Internal problem: undefined gridID in vlist "
                "%d, position %u!", vlistID, index);
      }

  if ( ! griddefined )
    {
      struct gridCompareSearchState query;
      query.queryKey = grid;// = { .queryKey = grid };
      if ( (gridglobdefined = (cdiResHFilterApply(&gridOps, gridCompareSearch, &query)
              == CDI_APPLY_STOP)) )
        gridID = query.resIDValue;

      if ( mode == 1 && gridglobdefined )
	for ( unsigned index = 0; index < ngrids; index++ )
	  if ( vlistptr->gridIDs[index] == gridID )
	    {
	      gridglobdefined = false;
	      break;
	    }
    }

  if ( ! griddefined )
    {
      if ( ! gridglobdefined )
        {
          grid->self = gridID = reshPut(grid, &gridOps);
          gridComplete(grid);
        }
      if ( mode < 2 )
        {
          vlistptr->gridIDs[ngrids] = gridID;
          vlistptr->ngrids++;
        }
    }

  return (struct addIfNewRes){ .Id = gridID, .isNew = !griddefined && !gridglobdefined };
}


const struct gridVirtTable cdiGridVtable
  = {
  .destroy = gridDestroyKernel,
  .copy = grid_copy_base,
  .copyScalarFields = grid_copy_base_scalar_fields,
  .copyArrayFields = grid_copy_base_array_fields,
  .defXVals = gridDefXValsSerial,
  .defYVals = gridDefYValsSerial,
  .defMask = gridDefMaskSerial,
  .defMaskGME = gridDefMaskGMESerial,
  .defXBounds = gridDefXBoundsSerial,
  .defYBounds = gridDefYBoundsSerial,
  .defArea = gridDefAreaSerial,
  .inqXVal = gridInqXValSerial,
  .inqYVal = gridInqYValSerial,
  .inqXVals = gridInqXValsSerial,
  .inqXCvals = gridInqXCvalsSerial,
  .inqXIsc = gridInqXIscSerial,
  .inqYVals = gridInqYValsSerial,
  .inqYCvals = gridInqYCvalsSerial,
  .inqYIsc = gridInqYIscSerial,
  .inqXValsPtr = gridInqXValsPtrSerial,
  .inqYValsPtr = gridInqYValsPtrSerial,
  .inqXCvalsPtr = gridInqXCvalsPtrSerial,
  .inqYCvalsPtr = gridInqYCvalsPtrSerial,
  .compareXYFull = compareXYvals,
  .compareXYAO = compareXYvals2,
  .inqArea = gridInqAreaSerial,
  .inqAreaPtr = gridInqAreaPtrBase,
  .hasArea = gridHasAreaBase,
  .inqMask = gridInqMaskSerial,
  .inqMaskGME = gridInqMaskGMESerial,
  .inqXBounds = gridInqXBoundsSerial,
  .inqYBounds = gridInqYBoundsSerial,
  .inqXBoundsPtr = gridInqXBoundsPtrSerial,
  .inqYBoundsPtr = gridInqYBoundsPtrSerial,
};

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
