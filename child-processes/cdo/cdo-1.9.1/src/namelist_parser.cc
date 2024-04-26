#include "cdo_int.h"
#include "namelist.h"


static
void kvlist_append_namelist(list_t *kvlist, const char *key, const char *buffer, namelisttok_t *t, int nvalues)
{
  char vbuf[4096];
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = NULL;

  if ( nvalues > 0 ) keyval->values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i )
    {
      size_t len = t[i].end - t[i].start;
      char *value = (char*) malloc((len+1)*sizeof(char));
      //printf(" value >%.*s<\n", len, buffer+t[i].start);
      const char *pval = buffer+t[i].start;
      if ( len < sizeof(vbuf) ) // snprintf seems to call strlen(pval)
        {
          memcpy(vbuf, buffer+t[i].start, len);
          vbuf[len] = 0;
          pval = vbuf;
        }
      snprintf(value, len+1, "%.*s", (int)len, pval);
      value[len] = 0;
      keyval->values[i] = value;
    }
  list_append(kvlist, &keyval);
}

static
int get_number_of_values(int ntok, namelisttok_t *tokens)
{
  int it;
  
  for ( it = 0; it < ntok; ++it )
    {
      namelisttok_t *t = &tokens[it];
      if ( t->type != NAMELIST_WORD && t->type != NAMELIST_STRING ) break;
    }
  
  return it;
}

static
int namelist_to_pml(list_t *pmlist, namelist_parser *parser, char *buf)
{
  char name[4096];
  list_t *kvlist = NULL;
  namelisttok_t *t;
  namelisttok_t *tokens = parser->tokens;
  unsigned int ntok = parser->toknext;
  // printf("Number of tokens %d\n", ntok);

  for ( unsigned int it = 0; it < ntok; ++it )
    {
      t = &tokens[it];
      // printf("Token %u", it+1);
      if ( t->type == NAMELIST_OBJECT )
        {
          name[0] = 0;
          if ( it+1 < ntok && tokens[it+1].type == NAMELIST_WORD )
            {
              it++;
              t = &tokens[it];
              snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf+t->start);
              name[sizeof(name)-1] = 0;
            }
          kvlist = kvlist_new(name);
          list_append(pmlist, &kvlist);
        }
      else if ( t->type == NAMELIST_KEY )
        {
          if ( kvlist == NULL )
            {
              kvlist = kvlist_new(NULL);
              list_append(pmlist, &kvlist);
            }
          // printf(" key >%.*s<\n", t->end - t->start, buf+t->start);
          snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf+t->start);
          name[sizeof(name)-1] = 0;
          int nvalues = get_number_of_values(ntok-it-1, &tokens[it+1]);
          kvlist_append_namelist(kvlist, name, buf, &tokens[it+1], nvalues);
          it += nvalues;
        }
      else
        {
          // printf(" token >%.*s<\n", t->end - t->start, buf+t->start);
          break;
        }
    }

  return 0;
}


list_t *namelistbuf_to_pmlist(listbuf_t *listbuf)
{
  const char *name = listbuf->name;
  namelist_parser *p = namelist_new();

  int status = namelist_parse(p, listbuf->buffer, listbuf->size);
  if ( status )
    {
      switch (status)
        {
        case NAMELIST_ERROR_INVAL: fprintf(stderr, "Namelist error: Invalid character in %s (line=%d character='%c')!\n", name, p->lineno, listbuf->buffer[p->pos]); break;
        case NAMELIST_ERROR_PART:  fprintf(stderr, "Namelist error: End of string not found in %s (line=%d)!\n", name, p->lineno); break;
        case NAMELIST_ERROR_INKEY: fprintf(stderr, "Namelist error: Invalid key word in %s (line=%d)!\n", name, p->lineno); break;
        case NAMELIST_ERROR_INTYP: fprintf(stderr, "Namelist error: Invalid key word type in %s (line=%d)!\n", name, p->lineno); break;
        case NAMELIST_ERROR_INOBJ: fprintf(stderr, "Namelist error: Invalid object in %s (line=%d)!\n", name, p->lineno); break;
        case NAMELIST_ERROR_EMKEY: fprintf(stderr, "Namelist error: Emtry key name in %s (line=%d)!\n", name, p->lineno); break;
        default:                   fprintf(stderr, "Namelist error in %s (line=%d)!\n", name, p->lineno); break;
        }
      cdoAbort("Namelist error!");
    }

  // namelist_dump(p, listbuf->buffer);
  status = namelist_verify(p, listbuf->buffer);
  if ( status )
    {
      fprintf(stderr, "Namelist error: Invalid contents in %s!\n", name);
      cdoAbort("Namelist error!");      
    }

  list_t *pmlist = list_new(sizeof(list_t *), free_kvlist, listbuf->name);

  namelist_to_pml(pmlist, p, listbuf->buffer);
  namelist_destroy(p);
  
  return pmlist;
}


list_t *namelist_to_pmlist(FILE *fp, const char *name)
{
  listbuf_t *listbuf = listbuf_new();
  if ( listbuf_read(listbuf, fp, name) ) cdoAbort("Read error on namelist %s!", name);
  
  list_t *pmlist = namelistbuf_to_pmlist(listbuf);
  if ( pmlist == NULL ) cdoAbort("Namelist not found!");

  listbuf_destroy(listbuf);

  return pmlist;
}
