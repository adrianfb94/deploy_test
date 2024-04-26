#include <errno.h>
#include "cdo_int.h"
#include "pmlist.h"
#include "json/jsmn.h"


static
char *readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  int ichar;
  size_t ipos = 0;

  while ( *buffersize )
    {
      ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if ( ichar == '\r' )
        {
          if ( *buffersize )
            {
              ichar = *buffer;
              if ( ichar == '\n' )
                {
                  (*buffersize)--;
                  buffer++;
                }
            }
          break;
        }
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %ld)!\n", len);
          break;
        }
    }

  line[ipos] = 0;

  if ( *buffersize == 0 && ipos == 0 ) buffer = NULL;

  return buffer;
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return pline;
}

static
char *getElementName(char *pline, char *name)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  size_t pos = 0;
  while ( pos < len && !isspace((int) *(pline+pos)) && *(pline+pos) != '=' && *(pline+pos) != ':' ) pos++;

  strncpy(name, pline, pos);
  name[pos] = 0;

  pline += pos;
  return pline;
}

static
char *getElementValue(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  if ( *pline != '"' && *pline != '\'' )
    for ( size_t i = 1; i < len; ++i )
      if ( pline[i] == '!' ) { pline[i] = 0; len = i; break; }
  while ( isspace((int) *(pline+len-1)) && len ) { *(pline+len-1) = 0; len--; }

  return pline;
}


void cmortablebuf_to_pmlist(list_t *pmlist, size_t buffersize, char *buffer)
{
  char line[4096];
  char name[256];
  char *pline;
  const char *listentry[] = {"axis_entry", "variable_entry"};
  int nentry = sizeof(listentry)/sizeof(listentry[0]);
  int linenumber = 0;
  list_t *kvlist = NULL;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);

      int ientry = -1;
      for ( ientry = 0; ientry < nentry; ++ientry )
        if ( strncmp(pline, listentry[ientry], strlen(listentry[ientry])) == 0 ) break;
      
      if ( ientry < nentry )
	{
	  pline += strlen(listentry[ientry]);

          kvlist = kvlist_new(listentry[ientry]);
          list_append(pmlist, &kvlist);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvlist, "name", (const char **)&pline, 1);
	}
      else
	{
	  pline = getElementName(pline, name);
	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( kvlist == NULL )
            {
              kvlist = kvlist_new("Header");
              list_append(pmlist, &kvlist);
            }

	  if ( *pline ) kvlist_append(kvlist, name, (const char **)&pline, 1);
	}
    }
}

// not used
int dump_json(const char *js, jsmntok_t *t, size_t count, int level)
{
  int i, j, k;
  if (count == 0)  return 0;

  if (t->type == JSMN_PRIMITIVE)
    {
      printf("%.*s", t->end - t->start, js+t->start);
      return 1;
    }
  else if (t->type == JSMN_STRING)
    {
      printf("'%.*s'", t->end - t->start, js+t->start);
      return 1;
    }
  else if (t->type == JSMN_OBJECT)
    {
      printf("\n");
      //  printf("Object: size %d\n", t->size);
      printf("Object: size %d count %d level %d\n", t->size, (int)count, level);
      j = 0;
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < level; k++) printf("  ");
          j += dump_json(js, t+1+j, count-j, level+1);
          printf(": ");
          j += dump_json(js, t+1+j, count-j, level+1);
          printf("\n");
        }
      return j+1;
    }
  else if (t->type == JSMN_ARRAY)
    {
      j = 0;
      printf("\n");
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < level-1; k++) printf("  ");
          printf("   - ");
          j += dump_json(js, t+1+j, count-j, level+1);
          printf("\n");
        }
      return j+1;
    }
  return 0;
}

static
void kvlist_append_json(list_t *kvlist, const char *key, const char *js, jsmntok_t *t, int nvalues)
{
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i )
    {
      size_t len = t[i].end - t[i].start;
      char *value = (char*) malloc((len+1)*sizeof(char));
      snprintf(value, len+1, "%.*s", (int)len, js+t[i].start);
      value[len] = 0;
      // printf("set %s: '%s'\n", key, value);
      keyval->values[i] = value;
    }
  list_append(kvlist, &keyval);
}

static
int json_to_pmlist(list_t *pmlist, const char *js, jsmntok_t *t, int count)
{
  bool debug = false;
  char name[4096];
  int i = 0;
  int nobj = t[0].size;
  if ( t[0].type == JSMN_OBJECT )
    while ( nobj-- )
      {
        ++i;
        int pmlname = i;
        if ( debug ) printf("  object: %.*s\n", t[i].end - t[i].start, js+t[i].start);
        ++i;
        if ( t[i].type == JSMN_OBJECT )
          {
            int ic = 0;
          NEXT:
            snprintf(name, sizeof(name), "%.*s", t[pmlname].end - t[pmlname].start, js+t[pmlname].start);
            name[sizeof(name)-1] = 0;
            // printf("new object: %s\n", name);
            list_t *kvlist = kvlist_new(name);
            list_append(pmlist, &kvlist);
                
            if ( t[i+2].type == JSMN_OBJECT )
              {
                if ( ic == 0 ) ic = t[i].size;
                else           ic--;
                
                ++i;
                kvlist_append_json(kvlist, "name", js, &t[i], 1);
                if ( debug ) printf("    name: '%.*s'\n", t[i].end - t[i].start, js+t[i].start);
                ++i;
              }
            int n = t[i].size;
            while ( n-- )
              {
                ++i;
                snprintf(name, sizeof(name), "%.*s", t[i].end - t[i].start, js+t[i].start);
                name[sizeof(name)-1] = 0;
                if ( debug ) printf("    %.*s:", t[i].end - t[i].start, js+t[i].start);
                ++i;
                if ( t[i].type == JSMN_ARRAY )
                  {
                    int nae = t[i].size;
                    kvlist_append_json(kvlist, name, js, &t[i+1], nae);
                    while ( nae-- )
                      {
                        ++i;
                        if ( debug ) printf(" '%.*s'", t[i].end - t[i].start, js+t[i].start);
                      }
                  }
                else
                  {
                    kvlist_append_json(kvlist, name, js, &t[i], 1);
                    if ( debug ) printf(" '%.*s'", t[i].end - t[i].start, js+t[i].start);
                  }
                if ( debug ) printf("\n");
              }
            if ( ic > 1 ) goto NEXT;
          }
      }

  if ( debug ) printf("Processed %d of %d tokens!\n", i, count-1);

  return 0;
}


void cmortablebuf_to_pmlist_json(list_t *pmlist, size_t buffersize, char *buffer, const char *filename)
{
  /* Prepare parser */
  jsmn_parser *p = jsmn_new();

  int status = jsmn_parse(p, buffer, buffersize);
  if ( status != 0 )
    {
      switch (status)
        {
        case JSMN_ERROR_INVAL: fprintf(stderr, "JSON error: Invalid character in %s (line=%d character='%c')!\n", filename, p->lineno, buffer[p->pos]); break;
        case JSMN_ERROR_PART:  fprintf(stderr, "JSON error: End of string not found in %s (line=%d)!\n", filename, p->lineno); break;
        default:               fprintf(stderr, "JSON error in %s (line=%d)\n", filename, p->lineno); break;
        }
    }

  json_to_pmlist(pmlist, buffer, p->tokens, (int)p->toknext);
  jsmn_destroy(p);
}


list_t *cmortable_to_pmlist(FILE *fp, const char *name)
{
  listbuf_t *listbuf = listbuf_new();
  if ( listbuf_read(listbuf, fp, name) ) cdoAbort("Read error on CMOR table %s!", name);
  
  list_t *pmlist = NULL;

  if ( listbuf->buffer[0] == '{' )
    {
      pmlist = list_new(sizeof(list_t *), free_kvlist, name);
      cmortablebuf_to_pmlist_json(pmlist, listbuf->size, listbuf->buffer, name);
    }
  else if ( strncmp(listbuf->buffer, "table_id:", 9) == 0 )
    {
      pmlist = list_new(sizeof(list_t *), free_kvlist, name);
      cmortablebuf_to_pmlist(pmlist, listbuf->size, listbuf->buffer);
    }
  else if ( listbuf->buffer[0] == '&' || listbuf->buffer[0] == '#' )
    {
      pmlist = namelistbuf_to_pmlist(listbuf);
    }
  else
    cdoAbort("Invalid CMOR table (file: %s)!", name);

  listbuf_destroy(listbuf);

  return pmlist;
}
