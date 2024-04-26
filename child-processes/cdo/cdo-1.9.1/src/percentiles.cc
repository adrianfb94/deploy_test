#include <math.h>
#include <string.h>
#include "util.h"
#include "percentiles.h"
#include "nth_element.h"

enum percentile_methods {NRANK=1, NIST, NUMPY};
enum interpolation_methods {LINEAR=1, LOWER, HIGHER, NEAREST};

static int percentile_method = NRANK;
static int interpolation_method = LINEAR;

static
double percentile_nrank(double *array, size_t len, double pn)
{
  size_t irank = (size_t)ceil(len*(pn/100.0));
  if ( irank <   1 ) irank = 1;
  if ( irank > len ) irank = len;
  return nth_element(array, len, irank-1);
}

static
double percentile_nist(double *array, size_t len, double pn)
{
  double rank = (len+1)*(pn/100.0);
  size_t k = (size_t) rank;
  double d = rank - k;
  double percentil = 0;
  if      ( k ==   0 ) percentil = nth_element(array, len, 0);
  else if ( k >= len ) percentil = nth_element(array, len, len-1);
  else
    {
      double vk1 = nth_element(array, len, k);
      double vk  = nth_element(array, len, k-1);;
      percentil = vk + d*(vk1 - vk);
    }

  return percentil;
}

static
double percentile_numpy(double *array, size_t len, double pn)
{
  double rank = (len-1)*(pn/100.0) + 1;
  size_t k = (size_t) rank;
  double d = rank - k;
  double percentil = 0;
  if      ( k ==   0 ) percentil = nth_element(array, len, 0);
  else if ( k >= len ) percentil = nth_element(array, len, len-1);
  else
    {
      if ( interpolation_method == LINEAR )
        {
          double vk1 = nth_element(array, len, k);
          double vk  = nth_element(array, len, k-1);;
          percentil = vk + d*(vk1 - vk);
        }
      else
        {
          size_t irank = 0;
          if      ( interpolation_method == LOWER   ) irank = (size_t) floor(rank);
          else if ( interpolation_method == HIGHER  ) irank = (size_t) ceil(rank);
          else if ( interpolation_method == NEAREST ) irank = (size_t) lround(rank);
          // numpy is using around(), with rounds to the nearest even value

          if ( irank <   1 ) irank = 1;
          if ( irank > len ) irank = len;

          percentil = nth_element(array, len, irank-1);
        }
    }

  return percentil;
}


double percentile(double *array, size_t len, double pn)
{
  double percentil = 0;

  if      ( percentile_method == NRANK ) percentil = percentile_nrank(array, len, pn);
  else if ( percentile_method == NIST  ) percentil = percentile_nist(array, len, pn);
  else if ( percentile_method == NUMPY ) percentil = percentile_numpy(array, len, pn);
  else cdoAbort("Internal error: percentile method %d not implemented!", percentile_method);

  return percentil;
}


void percentile_set_method(const char *methodstr)
{
  char *methodname = strdup(methodstr);
  strtolower(methodname);

  if      ( strcmp("nrank", methodname) == 0 ) percentile_method = NRANK;
  else if ( strcmp("nist",  methodname) == 0 ) percentile_method = NIST;
  else if ( strcmp("numpy", methodname) == 0 ) percentile_method = NUMPY;
  else if ( strcmp("numpy_linear",  methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = LINEAR;}
  else if ( strcmp("numpy_lower",   methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = LOWER;}
  else if ( strcmp("numpy_higher",  methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = HIGHER;}
  else if ( strcmp("numpy_nearest", methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = NEAREST;}
  else cdoAbort("Percentile method %s not available!", methodstr);
}


void percentile_check_number(double pn)
{
  if ( pn < 0 || pn > 100 )
    cdoAbort("Percentile number %g out of range! Percentiles must be in the range [0,100].", pn);
}

/*
  CDO check
#/bin/sh
CDO=cdo
#
cdo -f nc input,r5x1 testfile <<EOF
 15 20 35 40 50
EOF
cdo -f nc input,r6x1 testfile <<EOF
 15 20 35 40 50 55
EOF
#
PERS="30 40 50 75 100"
METS="nrank nist numpy numpy_lower numpy_higher numpy_nearest"
#
for MET in $METS; do
    for PER in $PERS; do
        echo "$MET: $PER"
        $CDO -s --percentile $MET output -fldpctl,$PER testfile 
    done
done
*/
/*
  numpy check
#python with numpy 1.9.0
import numpy as np
np.version.version
a=np.array([15, 20, 35, 40, 50, 55])
for p in [30, 40, 50, 75, 100] : print np.percentile(a, p, interpolation='linear')
*/
