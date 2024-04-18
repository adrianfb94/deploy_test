#ifndef _HETAETA_H
#define _HETAETA_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdbool.h>

void hetaeta(bool ltq, int ngp, const int *imiss,
	     int nlev1, const double *ah1, const double *bh1,
             const double *fis1, const double *ps1, 
             const double *t1, const double *q1,
             int nlev2, const double *ah2, const double *bh2, 
             const double *fis2, double *ps2, 
             double *t2, double *q2,
	     int nvars, double **vars1, double **vars2,
	     double *tscor, double *pscor, double *secor);

#endif  /* _HETAETA_H */
