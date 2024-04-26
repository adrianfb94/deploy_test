#ifndef CDI_ATT_H
#define CDI_ATT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _CDI_LIMITS_H
#include "cdi_limits.h"
#endif

/*
 * CDI attribute
 */
typedef struct {
  size_t    xsz;	  /* amount of space at xvalue                      */
  size_t    namesz;       /* size of name                                   */
  char     *name;         /* attribute name                                 */
  int       indtype;	  /* internal data type of xvalue (INT, FLT or TXT) */
  int       exdtype;      /* external data type                             */
                          /* indtype    exdtype                             */
                          /* TXT        TXT                                 */
                          /* INT        INT16, INT32                        */
                          /* FLT        FLT32, FLT64                        */
  size_t    nelems;    	  /* number of elements                             */
  void     *xvalue;       /* the actual data                                */
} cdi_att_t;


typedef struct {
  size_t     nalloc;		/* number allocated >= nelems */
  size_t     nelems;		/* length of the array */
  cdi_att_t  value[MAX_ATTRIBUTES];
} cdi_atts_t;


int cdiAttsGetSize(void *p, int varID, void *context);

void cdiAttsPack(void *p, int varID, void *buf, int size, int *position, void *context);

void cdiAttsUnpack(int cdiID, int varID, void *buf, int size, int *position, void *context);

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
