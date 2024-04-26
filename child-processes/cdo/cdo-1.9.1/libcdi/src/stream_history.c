#ifdef  HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"

#ifdef HAVE_LIBNETCDF
static inline bool
filetypeIsNetCDF(int filetype)
{
  return filetype == CDI_FILETYPE_NC
    ||   filetype == CDI_FILETYPE_NC2
    ||   filetype == CDI_FILETYPE_NC5
    ||   filetype == CDI_FILETYPE_NC4
    ||   filetype == CDI_FILETYPE_NC4C;
}
#endif

void streamDefHistory(int streamID, int length, const char *history)
{
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( filetypeIsNetCDF(streamptr->filetype) )
    {
      char *histstring;
      size_t len;
      if ( history )
	{
	  len = strlen(history);
	  if ( len )
	    {
              /* FIXME: what's the point of strdupx? Why not use
               * history argument directly? */
	      histstring = strdupx(history);
	      cdfDefHistory(streamptr, length, histstring);
	      Free(histstring);
	    }
	}
    }
#else
  (void)streamID; (void)length; (void)history;
#endif
}


int streamInqHistorySize(int streamID)
{
  int size = 0;
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( filetypeIsNetCDF(streamptr->filetype) )
    {
      size = cdfInqHistorySize(streamptr);
    }
#else
  (void)streamID;
#endif
  return (size);
}


void streamInqHistoryString(int streamID, char *history)
{
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( filetypeIsNetCDF(streamptr->filetype) )
    {
      cdfInqHistoryString(streamptr, history);
    }
#else
  (void)streamID; (void)history;
#endif
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
