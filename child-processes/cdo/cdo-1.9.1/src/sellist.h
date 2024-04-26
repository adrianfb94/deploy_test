#ifndef _SELLIST_H
#define _SELLIST_H

#include "pmlist.h"

typedef union {
  int ival;
  double dval;
  const char *cval;
} cvalues_t;

typedef struct {
  int nvalues;
  char *key;
  char **values;
  bool *flag;
  int type;
  char *txt;
  void *cvalues;
} selentry_t;


typedef struct {
  int size;
  selentry_t *entry;
} sellist_t;


#define  SELLIST_INT         1
#define  SELLIST_FLT         2
#define  SELLIST_WORD        3

#define  SELLIST_DEF_INT(name)                  int name = 0
#define  SELLIST_DEF_FLT(name)                  double name = 0
#define  SELLIST_DEF_WORD(name)                 const char *name = 0
#define  SELLIST_ADD_INT(name, txt)             SELLIST_DEF_INT(name);  int idx_##name = sellist_add(sellist, txt, #name, SELLIST_INT)
#define  SELLIST_ADD_FLT(name, txt)             SELLIST_DEF_FLT(name);  int idx_##name = sellist_add(sellist, txt, #name, SELLIST_FLT)
#define  SELLIST_ADD_WORD(name, txt)            SELLIST_DEF_WORD(name); int idx_##name = sellist_add(sellist, txt, #name, SELLIST_WORD)
#define  SELLIST_NVAL(name)                     sellist_nvalues(sellist, idx_##name)
#define  SELLIST_CHECK_FLAG(name)               sellist_check_flag(sellist, idx_##name)
#define  SELLIST_CHECK(name)                    sellist_check(sellist, idx_##name, &name)
#define  SELLIST_CHECK_DATE(name)               sellist_check_date(sellist, idx_##name, name)
#define  SELLIST_CHECK_SEASON(name, month)      sellist_check_season(sellist, idx_##name, month)
#define  SELLIST_DEF_FLAG(name, vindex, flag)   sellist_def_flag(sellist, idx_##name, vindex, flag)
#define  SELLIST_GET_VAL(name, vindex, val)     sellist_get_val(sellist, idx_##name, vindex, val)
#define  SELLIST_DEF_VAL(name, vindex, val)     sellist_def_val(sellist, idx_##name, vindex, val)


sellist_t *sellist_create(list_t *kvlist);
void sellist_destroy(sellist_t *sellist);
void sellist_verify(sellist_t *sellist);
int sellist_add(sellist_t *sellist, const char *txt, const char *name, int type);
int sellist_nvalues(sellist_t *sellist, int idx);
void sellist_check_flag(sellist_t *sellist, int idx);
bool sellist_check(sellist_t *sellist, int idx, void *par);
bool sellist_check_date(sellist_t *sellist, int idx, const char *par);
bool sellist_check_season(sellist_t *sellist, int idx, int month);
void sellist_def_flag(sellist_t *sellist, int idx, int vindex, bool flag);
void sellist_get_val(sellist_t *sellist, int idx, int vindex, void *val);
void sellist_def_val(sellist_t *sellist, int idx, int vindex, void *val);
void sellist_print(sellist_t *sellist);

#endif
