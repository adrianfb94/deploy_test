#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "dmemory.h"
#include "namespace.h"


/* the single image implementation */
static
int cdiStreamReadVar(int streamID, int varID, int memtype, void *data, int *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);
  int filetype = streamptr->filetype;

  *nmiss = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      {
        grb_read_var(streamptr, varID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        srvReadVarDP(streamptr, varID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        extReadVarDP(streamptr, varID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        iegReadVarDP(streamptr, varID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      {
        cdf_read_var(streamptr, varID, memtype, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }

  return status;
}

/*
@Function  streamReadVar
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void streamReadVar(int streamID, int varID, double *data, int *nmiss)
{
  cdiStreamReadVar(streamID, varID, MEMTYPE_DOUBLE, data, nmiss);
}

/*
@Function  streamReadVarF
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, float *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarF(int streamID, int varID, float *data, int *nmiss)
{
  if ( cdiStreamReadVar(streamID, varID, MEMTYPE_FLOAT, data, nmiss) )
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      elementCount *= (size_t) zaxisInqSize(vlistInqVarZaxis(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      streamReadVar(streamID, varID, conversionBuffer, nmiss);
      for ( size_t i = elementCount; i--; ) data[i] = (float) conversionBuffer[i];
      Free(conversionBuffer);
    }
}


static
int cdiStreamReadVarSlice(int streamID, int varID, int levelID, int memtype, void *data, int *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);
  int filetype = streamptr->filetype;

  *nmiss = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      {
        grb_read_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        srvReadVarSliceDP(streamptr, varID, levelID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        extReadVarSliceDP(streamptr, varID, levelID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        iegReadVarSliceDP(streamptr, varID, levelID, (double *)data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      {
        cdf_read_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
        break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
        status = 2;
	break;
      }
    }

  return status;
}

/*
@Function  streamReadVarSlice
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSlice reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
{
  if ( cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, data, nmiss) )
    {
      Warning("Unexpected error returned from cdiStreamReadVarSlice()!");
      size_t elementCount = (size_t)gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      memset(data, 0, elementCount * sizeof(*data));
    }
}

/*
@Function  streamReadVarSliceF
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSliceF reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, int *nmiss)
{
  if ( cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, data, nmiss) )
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      streamReadVarSlice(streamID, varID, levelID, conversionBuffer, nmiss);
      for ( size_t i = elementCount; i--; ) data[i] = (float) conversionBuffer[i];
      Free(conversionBuffer);
    }
}

static
int stream_read_record(int streamID, int memtype, void *data, int *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);

  *nmiss = 0;

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      grb_read_record(streamptr, memtype, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      srvReadRecord(streamptr, (double *)data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      extReadRecord(streamptr, (double *)data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      iegReadRecord(streamptr, (double *)data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      cdf_read_record(streamptr, memtype, data, nmiss);
      break;
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(streamptr->filetype));
	break;
      }
    }

  return status;
}


void streamReadRecord(int streamID, double *data, int *nmiss)
{
  stream_read_record(streamID, MEMTYPE_DOUBLE, (void *) data, nmiss);
}


void streamReadRecordF(int streamID, float *data, int *nmiss)
{
  if ( stream_read_record(streamID, MEMTYPE_FLOAT, (void *) data, nmiss) )
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      stream_t *streamptr = stream_to_pointer(streamID);
      int tsID   = streamptr->curTsID;
      int vrecID = streamptr->tsteps[tsID].curRecID;
      int recID  = streamptr->tsteps[tsID].recIDs[vrecID];
      int varID  = streamptr->tsteps[tsID].records[recID].varID;
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      streamReadRecord(streamID, conversionBuffer, nmiss);
      for ( size_t i = elementCount; i--; ) data[i] = (float) conversionBuffer[i];
      Free(conversionBuffer);
    }
}
