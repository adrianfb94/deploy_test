#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdbool.h>

#include "error.h"
#include "cdi.h"
#include "basetime.h"


void basetimeInit(basetime_t *basetime)
{
  if ( basetime == NULL )
    Error("Internal problem! Basetime not allocated.");

  basetime->ncvarid       = CDI_UNDEFID;
  basetime->ncdimid       = CDI_UNDEFID;
  basetime->ncvarboundsid = CDI_UNDEFID;
  basetime->leadtimeid    = CDI_UNDEFID;
  basetime->lwrf          = false;
  basetime->timevar_cache = NULL;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
