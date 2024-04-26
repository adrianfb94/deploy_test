#ifndef _GRIDDES_H
#define _GRIDDES_H

#include <stdbool.h>

typedef struct {
  int    *mask;
  double *xvals;
  double *yvals;
  double *xbounds;
  double *ybounds;
  double *area;
  double  xfirst, yfirst;
  double  xlast, ylast;
  double  xinc, yinc;
  double  xpole, ypole, angle;    /* rotated north pole             */
  int     scanningMode;
  /* scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
               64  = 128 * 0                + 64 *        1         + 32 * 0  
               00  = 128 * 0                + 64 *        0         + 32 * 0
               96  = 128 * 0                + 64 *        1         + 32 * 1
     Default / implicit scanning mode is 64:
                        i and j scan positively, i points are consecutive (row-major)        */
  bool    uvRelativeToGrid;
  double  a;
  int     datatype;
  int     isRotated;              /* TRUE for rotated grids         */
  int     type;
  int     ntr;
  int    *rowlon;
  bool    genBounds;
  int     nvertex;
  long    size;
  int     xsize;
  int     ysize;
  int     np;
  int     lcomplex;
  bool    def_xfirst;
  bool    def_yfirst;
  bool    def_xlast;
  bool    def_ylast;
  bool    def_xinc;
  bool    def_yinc;
  int     nd, ni, ni2, ni3;
  int     number, position;
  unsigned char uuid[CDI_UUID_SIZE];
  char    path[16384];
  char    xname[CDI_MAX_NAME];
  char    xlongname[CDI_MAX_NAME];
  char    xunits[CDI_MAX_NAME];
  char    xdimname[CDI_MAX_NAME];
  char    yname[CDI_MAX_NAME];
  char    ylongname[CDI_MAX_NAME];
  char    yunits[CDI_MAX_NAME];
  char    ydimname[CDI_MAX_NAME];
  char    vdimname[CDI_MAX_NAME];
}
griddes_t;

void gridInit(griddes_t *grid);
int gridDefine(griddes_t grid);

int gridFromNCfile(const char *gridfile);
int gridFromH5file(const char *gridfile);

#endif  /* _GRIDDES_H */
