#ifndef _ARRAY_H
#define _ARRAY_H

const char *fpe_errstr(int fpeRaised);

int array_minmaxsum_val(size_t len, const double *array, double *rmin, double *rmax, double *rsum);
int array_minmaxmean_val(size_t len, const double *array, double *rmin, double *rmax, double *rmean);

int array_mean_val(size_t len, const double *restrict array, double *rmean);
int array_mean_val_weighted(size_t len, const double *restrict array, const double *restrict w, double missval, double *rmean);

int array_add_array(size_t len, double *restrict array1, const double *restrict array2);

#endif // _ARRAY_H

