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
int cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data, int nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision writing.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0)
    Error("Writing of non-trivial subtypes not yet implemented!");

  // check taxis
  if ( streamptr->curTsID == CDI_UNDEFID ) streamDefTimestep(streamID, 0);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      {
        grb_write_var(streamptr, varID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        srvWriteVarDP(streamptr, varID, (double *)data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        extWriteVarDP(streamptr, varID, (double *)data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        iegWriteVarDP(streamptr, varID, (double *)data);
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
        cdf_write_var(streamptr, varID, memtype, data, nmiss);
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
@Function  streamWriteVar
@Title     Write a variable

@Prototype void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVar writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
{
  void (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype, const void *data, int nmiss)
    = (void (*)(int, int, int, const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;

  myCdiStreamWriteVar_(streamID, varID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

/*
@Function  streamWriteVarF
@Title     Write a variable

@Prototype void streamWriteVarF(int streamID, int varID, const float *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarF writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVarF(int streamID, int varID, const float *data, int nmiss)
{
  int (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype, const void *data, int nmiss)
    = (int (*)(int, int, int, const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;

  if ( myCdiStreamWriteVar_(streamID, varID, MEMTYPE_FLOAT, (const void *) data, nmiss) )
    {
      // In case the file format does not support single precision writing,
      // we fall back to double precision writing, converting the data
      // on the fly.
      int vlistID = streamInqVlist(streamID);
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(vlistID, varID));
      elementCount *= (size_t) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      for ( size_t i = elementCount; i--; ) conversionBuffer[i] = (double) data[i];
      myCdiStreamWriteVar_(streamID, varID, MEMTYPE_DOUBLE, (const void *) conversionBuffer, nmiss);
      Free(conversionBuffer);
    }
}

static
int cdiStreamWriteVarSlice(int streamID, int varID, int levelID, int memtype, const void *data, int nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision writing.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0)
    Error("Writing of non-trivial subtypes not yet implemented!");

  // check taxis
  if ( streamptr->curTsID == CDI_UNDEFID ) streamDefTimestep(streamID, 0);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      {
        grb_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        srvWriteVarSliceDP(streamptr, varID, levelID, (double *)data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        extWriteVarSliceDP(streamptr, varID, levelID, (double *)data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        iegWriteVarSliceDP(streamptr, varID, levelID, (double *)data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
      break;
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
@Function  streamWriteVarSlice
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSlice writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
{
  cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

/*
@Function  streamWriteVarSliceF
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSliceF writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, int nmiss)
{
  if ( cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, (const void *) data, nmiss) )
    {
      // In case the file format does not support single precision writing,
      // we fall back to double precision writing, converting the data on the fly.
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      for ( size_t i = elementCount; i--; ) conversionBuffer[i] = (double) data[i];
      streamWriteVarSlice(streamID, varID, levelID, conversionBuffer, nmiss);
      Free(conversionBuffer);
    }
}


void streamWriteVarChunk(int streamID, int varID,
                         const int rect[][2], const double *data, int nmiss)
{
  void (*myCdiStreamWriteVarChunk_)(int streamID, int varID, int memtype,
                                    const int rect[][2], const void *data, int nmiss)
    = (void (*)(int, int, int, const int [][2], const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_).func;
  myCdiStreamWriteVarChunk_(streamID, varID, MEMTYPE_DOUBLE, rect, data, nmiss);
}

/* single image implementation */
void cdiStreamWriteVarChunk_(int streamID, int varID, int memtype,
                             const int rect[][2], const void *data, int nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  stream_t *streamptr = stream_to_pointer(streamID);

  // streamDefineTaxis(streamID);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if defined (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
#endif
#if defined (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
#endif
#if defined (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
#endif
#if defined (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
#endif
#if  defined (HAVE_LIBGRIB) || defined (HAVE_LIBSERVICE)      \
  || defined (HAVE_LIBEXTRA) || defined (HAVE_LIBIEG)
      xabort("streamWriteVarChunk not implemented for filetype %s!",
             strfiletype(filetype));
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      cdf_write_var_chunk(streamptr, varID, memtype, rect, data, nmiss);
      break;
#endif
    default:
      Error("%s support not compiled in!", strfiletype(filetype));
      break;
    }
}

static
int stream_write_record(int streamID, int memtype, const void *data, int nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision writing.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      grb_write_record(streamptr, memtype, data, nmiss);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      srvWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      extWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      if ( memtype == MEMTYPE_FLOAT ) return 1;
      iegWriteRecord(streamptr, (const double *)data);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      {
	cdf_write_record(streamptr, memtype, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(streamptr->filetype));
	break;
      }
    }

  return status;
}

/*
@Function  streamWriteRecord
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteRecord(int streamID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteRecord writes the values of a horizontal slice (record) of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteRecord(int streamID, const double *data, int nmiss)
{
  stream_write_record(streamID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}


void streamWriteRecordF(int streamID, const float *data, int nmiss)
{
  if ( stream_write_record(streamID, MEMTYPE_FLOAT, (const void *) data, nmiss) )
    {
      // In case the file format does not support single precision writing,
      // we fall back to double precision writing, converting the data on the fly.
      stream_t *streamptr = stream_to_pointer(streamID);
      int varID = streamptr->record->varID;
      size_t elementCount = (size_t) gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount*sizeof(*conversionBuffer));
      for ( size_t i = elementCount; i--; ) conversionBuffer[i] = (double) data[i];
      streamWriteRecord(streamID, conversionBuffer, nmiss);
      Free(conversionBuffer);
    }
}

