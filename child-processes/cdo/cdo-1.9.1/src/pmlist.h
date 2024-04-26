#ifndef _PMLIST_H
#define _PMLIST_H

#include "list.h"

typedef struct {
  int nvalues;
  char *key;
  char **values;
} keyValues_t;

keyValues_t *kvlist_search(list_t *kvlist, const char *key);
list_t *pmlist_search_kvl(list_t *pmlist, const char *key, const char *value);

bool kvlist_print_iter(void *data);
bool pmlist_print_iter(void *data);

void kvlist_print(list_t *kvlist);

void free_keyval(void *data);
void free_kvlist(void *data);

list_t *kvlist_new(const char *name);
void kvlist_destroy(list_t *list);
void kvlist_append(list_t *kvlist, const char *key, const char **values, int nvalues);
int kvlist_parse_cmdline(list_t *kvlist, int nparams, char **params);

list_t *pmlist_search_kvlist_ventry(list_t *pmlist, const char *key, const char *value, int nentry, const char **entry);
list_t *pmlist_get_kvlist_ventry(list_t *pmlist, int nentry, const char **entry);

#endif
