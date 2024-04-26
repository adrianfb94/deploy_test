#ifndef _GRID_H
#define _GRID_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <math.h>
#include <stdbool.h>

extern double grid_missval;

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif


#ifndef  RAD2DEG
#define  RAD2DEG  (180./M_PI)   /* conversion for rad to deg */
#endif

#ifndef  DEG2RAD
#define  DEG2RAD  (M_PI/180.)   /* conversion for deg to rad */
#endif


int nfc_to_nlat(int nfc, int ntr);
int nlat_to_ntr(int nlat);
int nlat_to_ntr_linear(int nlat);
int ntr_to_nlat(int ntr);
int ntr_to_nlat_linear(int ntr);
int nlat_to_nlon(int nlat);

void grid_copy_attributes(int gridID1, int gridID2);
void grid_copy_mapping(int gridID1, int gridID2);

void grid_def_param_laea(int gridID, double a, double lon0, double lat0);
void grid_def_param_sinu(int gridID);

bool grid_is_distance_generic(int gridID);

void grid_to_radian(const char *units, size_t nvals, double *restrict values, const char *description);
void grid_to_degree(const char *units, size_t nvals, double *restrict values, const char *description);

void grid_gen_corners(size_t n, const double* restrict vals, double* restrict corners);
void grid_gen_bounds(size_t n, const double *restrict vals, double *restrict bounds);
void grid_check_lat_borders(int n, double *ybounds);

void grid_gen_xbounds2D(size_t nx, size_t ny, const double *restrict xbounds, double *restrict xbounds2D);
void grid_gen_ybounds2D(size_t nx, size_t ny, const double *restrict ybounds, double *restrict ybounds2D);

void grid_cell_center_to_bounds_X2D(const char* xunitstr, size_t xsize, size_t ysize,
				    const double *restrict grid_center_lon, double *restrict grid_corner_lon, double dlon);
void grid_cell_center_to_bounds_Y2D(const char* yunitstr, size_t xsize, size_t ysize,
				    const double *restrict grid_center_lat, double *restrict grid_corner_lat);

int gridWeights(int gridID, double *weights);
int gridGenArea(int gridID, double *area);

int referenceToGrid(int gridID);
int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToUnstructured(int gridID, int lbounds);
int gridToUnstructuredSelecton(int gridID1, size_t selectionSize, int *selectionIndexList, int nocoords ,int nobounds);
int gridToCurvilinear(int gridID, int lbounds);
int gridCurvilinearToRegular(int gridID);
int gridToRegular(int gridID);
void field2regular(int gridID1, int gridID2, double missval, double *array, int nmiss, int lnearest);

/* GME grid */
struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

void correct_sinxvals(int xsize, int ysize, double *xvals);

struct cart gc2cc(struct geo *position);
void factorni(int kni, int *kni2, int *kni3);
void gme_grid_restore(double *p, int ni, int nd);
void gme_grid(int lbounds, int gridsize, double *rlon, double *rlat,
	      double *blon, double *blat, int *imask,
              int ni, int nd, int ni2, int ni3);

/* Rotated grid */
double lamrot_to_lam(double phis, double rlas, double polphi, double pollam, double polgam);
double phirot_to_phi(double phis, double rlas, double polphi, double polgam);
void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v);

void cdo_print_grid(int gridID, int opt);

bool grid_has_proj4param(int gridID);

// Define a de-staggered grid for U and V
int cdo_define_destagered_grid(int gridID_u_stag, int gridID_v_stag, double *destagGridOffsets);

// Define a sampled grid of another grid
int cdo_define_sample_grid(int gridID, int sampleFactor);

// Define a sub-grid of another grid
int cdo_define_subgrid_grid(int gridSrcID, int subI0, int subI1, int subJ0, int subJ1);

#endif  /* _GRID_H */
