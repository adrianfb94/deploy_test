#ifndef _GRID_H
#define _GRID_H

#include "cdi.h"
#include <stdbool.h>

#include "cdi_att.h"

extern double grid_missval;
extern int (*proj_lonlat_to_lcc_func)();
extern int (*proj_lcc_to_lonlat_func)();

typedef unsigned char mask_t;

typedef struct grid_t grid_t;

struct gridVirtTable
{
  void (*destroy)(grid_t *gridptr);
  grid_t *(*copy)(grid_t *gridptr);
  void (*copyScalarFields)(grid_t *gridptrOrig, grid_t *gridptrDup);
  void (*copyArrayFields)(grid_t *gridptrOrig, grid_t *gridptrDup);
  void (*defXVals)(grid_t *gridptr, const double *xvals);
  void (*defYVals)(grid_t *gridptr, const double *yvals);
  void (*defMask)(grid_t *gridptr, const int *mask);
  void (*defMaskGME)(grid_t *gridptr, const int *mask);
  void (*defXBounds)(grid_t *gridptr, const double *xbounds);
  void (*defYBounds)(grid_t *gridptr, const double *ybounds);
  void (*defArea)(grid_t *gridptr, const double *area);
  double (*inqXVal)(grid_t *gridptr, int index);
  double (*inqYVal)(grid_t *gridptr, int index);
  int (*inqXVals)(grid_t *gridptr, double *xvals);
  int (*inqXCvals)(grid_t *gridptr, char **xcvals);
  int (*inqXIsc)(grid_t *gridptr);
  int (*inqYVals)(grid_t *gridptr, double *yvals);
  int (*inqYCvals)(grid_t *gridptr, char **ycvals);
  int (*inqYIsc)(grid_t *gridptr);
  const double *(*inqXValsPtr)(grid_t *gridptr);
  const char **(*inqXCvalsPtr)(grid_t *gridptr);
  const double *(*inqYValsPtr)(grid_t *gridptr);
  const char **(*inqYCvalsPtr)(grid_t *gridptr);
  /* return if for both grids, all xval and all yval are equal */
  bool (*compareXYFull)(grid_t *gridRef, grid_t *gridTest);
  /* return if for both grids, x[0], y[0], x[size-1] and y[size-1] are
   * respectively equal */
  bool (*compareXYAO)(grid_t *gridRef, grid_t *gridTest);
  void (*inqArea)(grid_t *gridptr, double *area);
  const double *(*inqAreaPtr)(grid_t *gridptr);
  int (*hasArea)(grid_t *gridptr);
  int (*inqMask)(grid_t *gridptr, int *mask);
  int (*inqMaskGME)(grid_t *gridptr, int *mask_gme);
  int (*inqXBounds)(grid_t *gridptr, double *xbounds);
  int (*inqYBounds)(grid_t *gridptr, double *ybounds);
  const double *(*inqXBoundsPtr)(grid_t *gridptr);
  const double *(*inqYBoundsPtr)(grid_t *gridptr);
};

struct gridaxis_t {
  char    name[CDI_MAX_NAME];
  char    longname[CDI_MAX_NAME];
  char    units[CDI_MAX_NAME];
  char    dimname[CDI_MAX_NAME];
  const char *stdname;
  int     size;                  // number of values
  short   flag;                  // 0: undefined 1:vals 2:first+inc
  double  first, last, inc;
  double *vals;
  int clength;
  char  **cvals;
  double *bounds;
};

// GME Grid
struct grid_gme_t {
  int     nd, ni, ni2, ni3;       /* parameter for GRID_GME         */
};

struct grid_t {
  char    vdimname[CDI_MAX_NAME];
  char    mapname[CDI_MAX_NAME];
  char    mapping[CDI_MAX_NAME];
  char   *name;
  int     self;
  int     size;
  int     type;                   /* grid type                      */
  int     datatype;               /* grid data type                 */
  int     proj;                   /* grid projection                */
  int     projtype;               /* grid projection type           */
  mask_t *mask;
  mask_t *mask_gme;
  double *area;
  struct grid_gme_t  gme;
  int     number, position;       /* parameter for GRID_REFERENCE   */
  int     trunc;                  /* parameter for GRID_SPECTEAL    */
  int     nvertex;
  char   *reference;
  unsigned char uuid[CDI_UUID_SIZE]; /* uuid for grid reference        */
  int    *rowlon;
  int     nrowlon;
  int     np;                     /* number of parallels between a pole and the equator */
  signed char isCyclic;           /* three possible states:
                                   * -1 if unknown,
                                   * 0 if found not cyclic, or
                                   * 1 for global cyclic grids
                                   */
  bool    lcomplex;
  bool    hasdims;
  bool uvRelativeToGrid;  /* Some models deliver wind U,V relative to the grid-cell */
  struct gridaxis_t x;
  struct gridaxis_t y;
  const struct gridVirtTable *vtable;
  cdi_atts_t atts;
  int  scanningMode;
  bool iScansNegatively, jScansPositively, jPointsAreConsecutive;
  /* scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
               64  = 128 * 0                + 64 *        1         + 32 * 0
               00  = 128 * 0                + 64 *        0         + 32 * 0
               96  = 128 * 0                + 64 *        1         + 32 * 1
     Default / implicit scanning mode is 64:
                        i and j scan positively, i points are consecutive (row-major)        */
};


void grid_init(grid_t *gridptr);
void cdiGridTypeInit(grid_t *gridptr, int gridtype, int size);
void grid_free(grid_t *gridptr);
grid_t *grid_to_pointer(int gridID);
extern const struct gridVirtTable cdiGridVtable;

unsigned cdiGridCount(void);

void gridVerifyProj(int gridID);

const double *gridInqXvalsPtr(int gridID);
const double *gridInqYvalsPtr(int gridID);

const char **gridInqXCvalsPtr(int gridID);
const char **gridInqYCvalsPtr(int gridID);

const double *gridInqXboundsPtr(int gridID);
const double *gridInqYboundsPtr(int gridID);
const double *gridInqAreaPtr(int gridID);

const char *gridInqReferencePtr(int gridID);

int gridGenerate(const grid_t *grid);

//int gridIsEqual(int gridID1, int gridID2);

void cdiGridGetIndexList(unsigned, int * );

void
gridUnpack(char * unpackBuffer, int unpackBufferSize,
           int * unpackBufferPos, int originNamespace, void *context,
           int force_id);

struct addIfNewRes
{
  int Id;
  int isNew;
};

struct addIfNewRes cdiVlistAddGridIfNew(int vlistID, grid_t *grid, int mode);

int gridVerifyGribParamLCC(double missval, double *lon_0, double *lat_0, double *lat_1, double *lat_2,
                           double *a, double *rf, double *xval_0, double *yval_0, double *x_0, double *y_0);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
