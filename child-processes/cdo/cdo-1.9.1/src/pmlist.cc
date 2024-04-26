#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pmlist.h"

keyValues_t *kvlist_search(list_t *kvlist, const char *key)
{
  if ( key )
    {
      listNode_t *node = kvlist->head;
      while ( node )
        {
          keyValues_t *kv = *(keyValues_t **)node->data;
          if ( kv->key && *(kv->key) == *key && strcmp(kv->key, key) == 0 ) return kv;
          node = node->next;
        }
    }

  return NULL;
}


list_t *pmlist_search_kvlist(list_t *pmlist, const char *key, const char *value)
{
  if ( pmlist && key && value )
    {
      listNode_t *node = pmlist->head;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvlist = *(list_t **)node->data;
              keyValues_t *kv = kvlist_search(kvlist, key);
              if ( kv && kv->nvalues > 0 && *(kv->values[0]) == *value && strcmp(kv->values[0], value) == 0 ) return kvlist;
            }
          node = node->next;
        }
    }

  return NULL;
}


bool kvlist_print_iter(void *data)
{
  keyValues_t *keyval = *(keyValues_t **)data;
  char *key = keyval->key;
  char **values = keyval->values;
  int nvalues = keyval->nvalues;
  printf("  %s =", key);
  for ( int i = 0; i < nvalues; ++i ) printf(" '%s'", values[i]);
  printf("\n");

  return true;
}


void kvlist_print(list_t *kvlist)
{
  printf("Key/Value list %s:\n", list_name(kvlist));
  list_for_each(kvlist, kvlist_print_iter);
}


bool pmlist_print_iter(void *data)
{
  list_t *kvlist = *(list_t **)data;
  const char *listname = list_name(kvlist);
  printf("\nFound %s list with %d key/values: \n", listname?listname:"", list_size(kvlist));
  list_for_each(kvlist, kvlist_print_iter);
  return true;
}


void free_keyval(void *data)
{
  keyValues_t *keyval = *(keyValues_t **)data;
  if ( keyval->key ) free(keyval->key);
  int nvalues = keyval->nvalues;
  for ( int i = 0; i < nvalues; ++i )
    if ( keyval->values[i] ) free(keyval->values[i]);
  free(keyval->values);
  free(keyval);
}


void free_kvlist(void *data)
{
  list_t *kvlist = *(list_t **)data;
  //int n = list_size(kvlist);
  list_destroy(kvlist);
  //printf("Successfully freed %d keyvalues...\n", n);
}


list_t *kvlist_new(const char *name)
{
  return list_new(sizeof(keyValues_t *), free_keyval, name);
}


void kvlist_destroy(list_t *list)
{
  list_destroy(list);
}


void kvlist_append(list_t *kvlist, const char *key, const char **values, int nvalues)
{
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i ) keyval->values[i] = strdup(values[i]);
  list_append(kvlist, &keyval);
}


int kvlist_parse_cmdline(list_t *kvlist, int nparams, char **params)
{
  /* Assume key = value pairs. That is, if params[i] contains no '='
   * then treat it as if it belongs to the values of params[i-1]. */
  char key[256];
  int i = 0;
  while ( i < nparams )
    {
      char *end = strchr(params[i], '=');
      if ( end == NULL )
        {
          fprintf(stderr, "Missing '=' in key/value string: >%s<\n", params[i]);
          return -1;
        }

      snprintf(key, sizeof(key), "%.*s", (int)(end-params[i]), params[i]);
      key[sizeof(key)-1] = 0;

      int j = 1;
      while ( i + j < nparams && strchr(params[i + j], '=') == NULL ) j++;

      int nvalues = j;

      const char **values = nvalues ? (const char**) malloc(nvalues*sizeof(char*)) : NULL;

      values[0] = end + 1;
      if ( *values[0] == 0 ) nvalues = 0;

      for ( j = 1; j < nvalues; ++j ) values[j] = params[i + j];
      kvlist_append(kvlist, key, values, nvalues);

      if ( values ) free(values);
      
      i += j;
    }

  return 0;
}


list_t *pmlist_search_kvlist_ventry(list_t *pmlist, const char *key, const char *value, int nentry, const char **entry)
{
  if ( pmlist && key && value )
    {
      listNode_t *node = pmlist->head;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvlist = *(list_t **)node->data;
              const char *listname = list_name(kvlist);
              for ( int i = 0; i < nentry; ++i )
                if ( strcmp(listname, entry[i]) == 0 )
                  {
                    keyValues_t *kv = kvlist_search(kvlist, key);
                    if ( kv && kv->nvalues > 0 && *(kv->values[0]) == *value && strcmp(kv->values[0], value) == 0 ) return kvlist;
                  }
            }
          node = node->next;
        }
    }

  return NULL;
}


list_t *pmlist_get_kvlist_ventry(list_t *pmlist, int nentry, const char **entry)
{
  if ( pmlist )
    {
      listNode_t *node = pmlist->head;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvlist = *(list_t **)node->data;
              const char *listname = list_name(kvlist);
              for ( int i = 0; i < nentry; ++i )
                if ( strcmp(listname, entry[i]) == 0 ) return kvlist;
            }
          node = node->next;
        }
    }

  return NULL;
}

/*
int main(void)
{
  int numLists = 0;
  printf("Generating list with lists of keyValues...\n");

  list_t *pmlist = list_new(sizeof(list_t *), free_kvlist, "parameter");

  {
    const char *k1 = "longname", *k1vals[] = {"surface temperature"};
    const char *k2 = "name",     *k2vals[] = {"temperature"};
    const char *k3 = "values",   *k3vals[] = {"273.15", "292.5", "301.4"};

    list_t *kvlist = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvlist, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvlist, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvlist, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));

    list_append(pml, &kvlist);
    numLists++;
  }
  {
    const char *k1 = "longname", *k1vals[] = {"surface pressure"};
    const char *k2 = "name",     *k2vals[] = {"pressure"};
    const char *k3 = "values",   *k3vals[] = {"1000", "850", "500"};
    const char *k4 = "units",    *k4vals[] = {"hPa"};

    list_t *kvlist = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvlist, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvlist, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvlist, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvlist, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));

    list_append(pmlist, &kvlist);
    numLists++;
  }
  {
    const char *k1 = "longname", *k1vals[] = {"Air Temperature"};
    const char *k2 = "name",     *k2vals[] = {"ta"};
    const char *k3 = "valid_max",*k3vals[] = {"336"};
    const char *k4 = "units",    *k4vals[] = {"K"};

    list_t *kvlist = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvlist, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvlist, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvlist, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvlist, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));

    list_append(pmlist, &kvlist);
    numLists++;
  }
  
 
  printf("pmlist=%s n=%d:\n", list_name(pmlist), list_size(pmlist));
  list_for_each(pml, pmlist_print_iter);

  list_t *kvlist = pmlist_search_kvl(pml, "name", "pressure");
  printf("kvlist of name=pressure:\n");
  if ( kvlist ) list_for_each(kvlist, kvlist_print_iter);

  printf("\n");
  int n = list_size(kvlist);
  for ( int i = 0; i < n; ++i )
    {
      keyValues_t *kv = (keyValues_t *) list_entry(kvlist, i);
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
    }

  printf("\n");
  int i = 0;
  for ( listNode_t *node = kvlist->head; node; node = node->next )
    {
      keyValues_t *kv = *(keyValues_t **)node->data;
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
      ++i;
    }

  list_destroy(pmlist);
  printf("Successfully freed %d kvlists...\n", numLists);

  return 0;
}
*/
