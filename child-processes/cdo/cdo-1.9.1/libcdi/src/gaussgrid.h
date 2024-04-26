#ifndef _GAUSSGRID_H
#define _GAUSSGRID_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

void gaussaw(double *restrict pa, double *restrict pw, size_t nlat);
bool isGaussGrid(size_t ysize, double yinc, const double *yvals);

#if defined (__cplusplus)
}
#endif

#endif  /* _GAUSSGRID_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
