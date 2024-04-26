#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "listbuf.h"

listbuf_t *listbuf_new(void)
{
  listbuf_t *listbuf = (listbuf_t *) malloc(sizeof(listbuf_t));
  if ( listbuf )
    {
      listbuf->size = 0;
      listbuf->buffer = NULL;
      listbuf->name = NULL;
    }

  return listbuf;
}


int listbuf_read(listbuf_t *listbuf, FILE *fp, const char *name)
{
  int filedes = fileno(fp);
  struct stat buf;
  size_t filesize = 0;
  if ( fstat(filedes, &buf) == 0 ) filesize = (size_t) buf.st_size;

  if ( filesize > MAX_LISTBUF_SIZE )
    {
      fprintf(stderr, "%s: max buffer size of %d exceeded (%s)!\n", __func__, (int)MAX_LISTBUF_SIZE, name);
      return -1;
    }
  else if ( filesize == 0 )
    {
      fprintf(stderr, "%s: empty stream: %s\n", __func__, name);
      return -1;
    }

  char *buffer = (char*) malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  if ( nitems != filesize )
    {
      free(buffer);
      fprintf(stderr, "%s: read failed on %s!\n", __func__, name);
      return -1;
    }

  listbuf->size = filesize;
  listbuf->buffer = buffer;

  if ( name ) listbuf->name = strdup(name);

  return 0;
}


void listbuf_destroy(listbuf_t *listbuf)
{
  if ( listbuf )
    {
      if ( listbuf->buffer ) free(listbuf->buffer);
      if ( listbuf->name ) free(listbuf->name);
      listbuf->size = 0;
      listbuf->buffer = NULL;
      listbuf->name = NULL;
      free(listbuf);
    }
}


