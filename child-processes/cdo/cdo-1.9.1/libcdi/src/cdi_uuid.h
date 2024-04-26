#ifndef CDI_UUID_H
#define CDI_UUID_H

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "cdi.h"


#ifdef __cplusplus
extern "C" {
#endif

static inline int cdiUUIDIsNull(const unsigned char uuid[])
{
  int isNull = 1;
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i)
    isNull &= (uuid[i] == 0);
  return isNull;
}

void cdiCreateUUID(unsigned char uuid[CDI_UUID_SIZE]);

void cdiUUID2Str(const unsigned char uuid[], char uuidstr[]);
int cdiStr2UUID(const char *uuidstr, unsigned char uuid[]);

#if defined (__cplusplus)
}
#endif

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
