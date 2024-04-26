#include <cdo_int.h>
#include "sellist.h"

//#define SELDEBUG 1

sellist_t *sellist_create(list_t *kvlist)
{
  sellist_t *sellist = (sellist_t *) Malloc(sizeof(sellist_t));
  sellist->size = list_size(kvlist);
  sellist->entry = (selentry_t *) Malloc(sellist->size*sizeof(selentry_t));

  int i = 0;
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      selentry_t *e = &(sellist->entry[i]);
      e->key = kv->key;
      e->values = kv->values;
      e->nvalues = kv->nvalues;
#ifdef SELDEBUG
      printf("%s =", e->key);
      for ( int ii = 0; ii < e->nvalues; ++ii ) printf(" '%s'", e->values[ii]);
      printf("\n");
#endif
      ++i;
    }

  for ( int i = 0; i < sellist->size; ++i )
    {
      selentry_t *e = &(sellist->entry[i]);
      e->flag = NULL;
      e->cvalues = NULL;
#ifdef SELDEBUG
      printf("%s =", e->key);
      for ( int ii = 0; ii < e->nvalues; ++ii ) printf(" '%s'", e->values[ii]);
      printf("\n");
#endif
    }

  return sellist;
}


void sellist_destroy(sellist_t *sellist)
{
  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          selentry_t *e = &(sellist->entry[i]);
          if ( e->txt ) Free(e->txt);
          if ( e->flag ) Free(e->flag);
          if ( e->cvalues ) Free(e->cvalues);
        }

      Free(sellist);
    }
}


void sellist_verify(sellist_t *sellist)
{
  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          selentry_t *e = &(sellist->entry[i]);
          if ( e->type == 0 ) cdoAbort("Unsupported selection keyword: '%s'!", e->key);
        }
    }
}

void split_intstring(const char *intstr, int *first, int *last, int *inc);

int sellist_add(sellist_t *sellist, const char *txt, const char *name, int type)
{
  int idx = -1;

  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          const char *key = sellist->entry[i].key;
          if ( strcmp(key, name) == 0 )
            {
              idx = i;
              break;
            }
        }

      if ( idx >= 0 && idx < sellist->size )
        {
          selentry_t *e = &(sellist->entry[idx]);
          e->type = type;
          e->txt = strdup(txt);
          if ( e->nvalues && e->cvalues == NULL )
            {
              switch (type)
                {
                case SELLIST_INT:  e->cvalues = Malloc(e->nvalues*sizeof(int)); break;
                case SELLIST_FLT:  e->cvalues = Malloc(e->nvalues*sizeof(double)); break;
                case SELLIST_WORD: e->cvalues = Malloc(e->nvalues*sizeof(char*)); break;
                }
            }

          int j = 0;
          int nvalues = e->nvalues;
          for ( int i = 0; i < nvalues; ++i )
            switch (type)
              {
              case SELLIST_INT:
                {
                  int first, last, inc;
                  split_intstring(e->values[i], &first, &last, &inc);

                  if ( first == last )
                    {
                      ((int*)e->cvalues)[j++] = parameter2int(e->values[i]);
                    }
                  else
                    {
                      int k = 0;
                      if ( inc >= 0 )
                        for ( int ival = first; ival <= last; ival += inc ) k++;
                      else
                        for ( int ival = first; ival >= last; ival += inc ) k++;

                      e->nvalues += k-1;
                      e->cvalues = Realloc(e->cvalues, e->nvalues*sizeof(int));

                      if ( inc >= 0 )
                        {
                          for ( int ival = first; ival <= last; ival += inc )
                            ((int*)e->cvalues)[j++] = ival;
                        }
                      else
                        {
                          for ( int ival = first; ival >= last; ival += inc )
                            ((int*)e->cvalues)[j++] = ival;
                        }
                    }

                  break;
                }
              case SELLIST_FLT:  ((double*)e->cvalues)[i] = parameter2double(e->values[i]); break;
              case SELLIST_WORD: ((const char**)e->cvalues)[i] = parameter2word(e->values[i]); break;
              }

          if ( e->nvalues ) e->flag = (bool*) Calloc(e->nvalues, sizeof(bool));      
#ifdef SELDEBUG          
          printf("%s =", e->key);
          for ( int i = 0; i < e->nvalues; ++i )
            switch (type)
              {
              case SELLIST_INT:  printf(" %d", ((int*)e->cvalues)[i]); break;
              case SELLIST_FLT:  printf(" %g", ((double*)e->cvalues)[i]); break;
              case SELLIST_WORD: printf(" %s", ((char**)e->cvalues)[i]); break;
              }
          printf("\n");
#endif
        }
    }

  return idx;
}


int sellist_nvalues(sellist_t *sellist, int idx)
{
  int nvalues = 0;

  if ( sellist && idx >= 0 && idx < sellist->size ) nvalues = sellist->entry[idx].nvalues;

  return nvalues;
}


void sellist_check_flag(sellist_t *sellist, int idx)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      for ( int i = 0; i < nvalues; ++i )
        if ( e->flag[i] == false )
          switch (e->type)
            {
            case SELLIST_INT:  cdoWarning("%s >%d< not found!", e->txt, ((int*)e->cvalues)[i]); break;
            case SELLIST_FLT:  cdoWarning("%s >%g< not found!", e->txt, ((double*)e->cvalues)[i]); break;
            case SELLIST_WORD: cdoWarning("%s >%s< not found!", e->txt, ((char**)e->cvalues)[i]); break;
            }
    }
}


bool sellist_check(sellist_t *sellist, int idx, void *par)
{
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      for ( int i = 0; i < nvalues; ++i )
        {
          switch (e->type)
            {
            case SELLIST_INT:  if ( *(int*)par == ((int*)e->cvalues)[i] )                       { found = true; e->flag[i] = true; } break;
            case SELLIST_FLT:  if ( fabs(*(double*)par - ((double*)e->cvalues)[i]) < 1.e-4 )    { found = true; e->flag[i] = true; } break;
            case SELLIST_WORD: if ( wildcardmatch(((char**)e->cvalues)[i], *(char**)par) == 0 ) { found = true; e->flag[i] = true; } break;
            }
        }
    }

  return found;
}


bool sellist_check_date(sellist_t *sellist, int idx, const char *par)
{
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      char wcdate[512];
      selentry_t *e = &(sellist->entry[idx]);

      if ( *par == ' ' ) ++par;

      for ( int i = 0; i < nvalues; ++i )
        {
          strcpy(wcdate, e->values[i]);
          strcat(wcdate, "*");
          if ( wildcardmatch(wcdate, par) == 0 ) { found = true; e->flag[i] = true; }
        }
    }

  return found;
}

void season_to_months(const char *season, int *imonths);

bool sellist_check_season(sellist_t *sellist, int idx, int month)
{
  assert(month>=1&&month<=12);
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      int imon[13]; /* 1-12 ! */
      selentry_t *e = &(sellist->entry[idx]);

      for ( int i = 0; i < nvalues; ++i )
        {
          for ( int m = 0; m < 13; ++m ) imon[m] = 0;
          season_to_months(e->values[i], imon);
          if ( imon[month] ) { found = true; e->flag[i] = true; }
        }
    }

  return found;
}


void sellist_def_flag(sellist_t *sellist, int idx, int vindex, bool flag)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      if ( vindex >= 0 && vindex < nvalues ) e->flag[vindex] = flag;
    }  
}


void sellist_get_val(sellist_t *sellist, int idx, int vindex, void *val)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      if ( vindex >= 0 && vindex < nvalues )
        {
          switch (e->type)
            {
            case SELLIST_INT:  *(int*)val = ((int*)e->cvalues)[vindex]; break;
            case SELLIST_FLT:  *(double*)val = ((double*)e->cvalues)[vindex]; break;
            case SELLIST_WORD: *(const char**)val = ((const char**)e->cvalues)[vindex]; break;
            }
        }
    }
}


void sellist_def_val(sellist_t *sellist, int idx, int vindex, void *val)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      if ( vindex >= 0 && vindex < nvalues )
        {
          switch (e->type)
            {
            case SELLIST_INT:  ((int*)e->cvalues)[vindex] = *(int*)val; break;
            case SELLIST_FLT:  ((double*)e->cvalues)[vindex] = *(double*)val; break;
            case SELLIST_WORD: ((const char**)e->cvalues)[vindex] = *(const char**)val; break;
            }
        }
    }
}

static
void sellist_print_val(int type, cvalues_t *cvalues, int i)
{
  switch (type)
    {
    case SELLIST_INT:  printf(" %d", ((int*)cvalues)[i]); break;
    case SELLIST_FLT:  printf(" %g", ((double*)cvalues)[i]); break;
    case SELLIST_WORD: printf(" %s", ((char**)cvalues)[i]); break;
    }
}


void sellist_print(sellist_t *sellist)
{
  if ( sellist )
    {
      // printf("Parameter list: %s\n", sellist->
      printf("Num  Name             Type  Size  Entries\n");
      for ( int idx = 0; idx < sellist->size; ++idx )
        {
          selentry_t *e = &(sellist->entry[idx]);
          printf("%3d  %-16s %4d  %4d ", idx+1, e->key, e->type, e->nvalues);
          int nvalues = e->nvalues;
          if ( nvalues > 12 ) nvalues = 11;
          for ( int i = 0; i < nvalues; ++i ) sellist_print_val(e->type, (cvalues_t *)e->cvalues, i);
          if ( nvalues < e->nvalues ) printf(" ...");
          sellist_print_val(e->type, (cvalues_t *)e->cvalues, e->nvalues-1);
          printf("\n");
        }
    }
}

