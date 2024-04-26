#ifndef __CONVERT_UNITS_H_
#define __CONVERT_UNITS_H_

#if  defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBUDUNITS2) && (defined(HAVE_UDUNITS2_H) || defined(HAVE_UDUNITS2_UDUNITS2_H))
#define HAVE_UDUNITS2
#endif

#if defined(HAVE_UDUNITS2)
#if defined(HAVE_UDUNITS2_UDUNITS2_H)
#  include <udunits2/udunits2.h>
#else
#  include <udunits2.h>
#endif

void cdoConvertFree(void *ut_converter);
void cdoConvertDestroy();
#endif

#include "cdo_int.h"

void cdoConvertUnits(void **ut_converter, bool *changeunits, char *units, char *units_old, const char *name);

#endif // __CONVERT_UNITS_H_
