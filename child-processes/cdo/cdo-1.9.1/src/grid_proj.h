#ifndef _GRID_PROJ_H
#define _GRID_PROJ_H


int cdo_lonlat_to_lcc(int gridID, size_t nvals, double *xvals, double *yvals);
int cdo_lcc_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_sinu_to_lonlat(size_t nvals, double *xvals, double *yvals);
void cdo_laea_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_proj_to_lonlat(char *proj4param, size_t nvals, double *xvals, double *yvals);

#ifdef __cplusplus
extern "C" {
#endif
int proj_lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, size_t nvals, double *xvals, double *yvals);
int proj_lcc_to_lonlat(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, double x_0, double y_0, size_t nvals, double *xvals, double *yvals);
#if defined (__cplusplus)
}
#endif

#endif  /* _GRID_PROJ_H */
