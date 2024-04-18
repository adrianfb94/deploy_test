#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBGRIB

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cgribex.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "file.h"
#include "cgribex.h"  /* gribZip gribGetZip gribGinfo */
#include "namespace.h"


static
int grbDecode(int filetype, int memtype, void *gribbuffer, int gribsize, void *data, size_t datasize,
	      int unreduced, int *nmiss, double missval, int vlistID, int varID)
{
  int status = 0;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == CDI_FILETYPE_GRB )
    {
#if  defined  (HAVE_LIBGRIB_API)
      extern int cdiNAdditionalGRIBKeys;
      if ( cdiNAdditionalGRIBKeys > 0 )
	Error("CGRIBEX decode does not support reading of additional GRIB keys!");
#endif
      status = cgribexDecode(memtype, gribbuffer, gribsize, data, (long) datasize, unreduced, nmiss, missval);
    }
  else
#endif
#ifdef HAVE_LIBGRIB_API
    {
      void *datap = data;
      if ( memtype == MEMTYPE_FLOAT )
        datap = Malloc(datasize*sizeof(double));

      status = gribapiDecode(gribbuffer, gribsize, datap, (long) datasize, unreduced, nmiss, missval, vlistID, varID);

      if ( memtype == MEMTYPE_FLOAT )
        {
          float *dataf = (float*) data;
          double *datad = (double*) datap;
          for ( size_t i = 0; i < datasize; ++i ) dataf[i] = (float) datad[i];
          Free((void*)datap);
        }
    }
#else
    {
      (void)vlistID; (void)varID;
      Error("GRIB_API support not compiled in!");
    }
#endif

  return status;
}

static
int grbUnzipRecord(void *gribbuffer, size_t *gribsize)
{
  int zip = 0;
  int izip;
  size_t unzipsize;

  size_t igribsize = *gribsize;
  size_t ogribsize = *gribsize;

  if ( (izip = gribGetZip(igribsize, (unsigned char *)gribbuffer, &unzipsize)) > 0 )
    {
      zip = izip;
      if ( izip == 128 ) /* szip */
	{
	  unsigned char *itmpbuffer = NULL;
	  size_t itmpbuffersize = 0;

	  if ( unzipsize < igribsize )
	    {
	      fprintf(stderr, "Decompressed size smaller than compressed size (in %zu; out %zu)!\n", igribsize, unzipsize);
	      return 0;
	    }

	  if ( itmpbuffersize < igribsize )
	    {
	      itmpbuffersize = igribsize;
	      itmpbuffer = (unsigned char *) Realloc(itmpbuffer, itmpbuffersize);
	    }

	  memcpy(itmpbuffer, gribbuffer, itmpbuffersize);

	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */

	  ogribsize = (size_t)gribUnzip((unsigned char *)gribbuffer, (long)unzipsize, itmpbuffer, (long)igribsize);

	  Free(itmpbuffer);

	  if ( ogribsize <= 0 ) Error("Decompression problem!");
	}
      else
	{
	  Error("Decompression for %d not implemented!", izip);
	}
    }

  *gribsize = ogribsize;

  return zip;
}


void grb_read_record(stream_t * streamptr, int memtype, void *data, int *nmiss)
{
  int filetype = streamptr->filetype;

  void *gribbuffer = streamptr->record->buffer;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int tsID    = streamptr->curTsID;
  int vrecID  = streamptr->tsteps[tsID].curRecID;
  int recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  off_t recpos  = streamptr->tsteps[tsID].records[recID].position;
  size_t recsize = streamptr->tsteps[tsID].records[recID].size;
  int varID   = streamptr->tsteps[tsID].records[recID].varID;

  int gridID   = vlistInqVarGrid(vlistID, varID);
  int gridsize = gridInqSize(gridID);

  streamptr->numvals += gridsize;

  fileSetPos(fileID, recpos, SEEK_SET);

  if (fileRead(fileID, gribbuffer, recsize) != recsize)
    Error("Failed to read GRIB record");

  double missval = vlistInqVarMissval(vlistID, varID);

  streamptr->tsteps[tsID].records[recID].zip = grbUnzipRecord(gribbuffer, &recsize);

  grbDecode(filetype, memtype, gribbuffer, (int)recsize, data, (size_t)gridsize, streamptr->unreduced, nmiss, missval, vlistID, varID);
}


void grb_read_var(stream_t * streamptr, int varID, int memtype, void *data, int *nmiss)
{
  int filetype = streamptr->filetype;

  void *gribbuffer = streamptr->record->buffer;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int tsID    = streamptr->curTsID;

  int gridID   = vlistInqVarGrid(vlistID, varID);
  int gridsize = gridInqSize(gridID);

  off_t currentfilepos = fileGetPos(fileID);

  int isub     = subtypeInqActiveIndex(streamptr->vars[varID].subtypeID);
  int nlevs    = streamptr->vars[varID].recordTable[0].nlevs;
  if ( CDI_Debug )
    Message("nlevs = %d gridID = %d gridsize = %d", nlevs, gridID, gridsize);
  *nmiss = 0;
  for (int levelID = 0; levelID < nlevs; levelID++ )
    {
      int    recID   = streamptr->vars[varID].recordTable[isub].recordID[levelID];
      off_t  recpos  = streamptr->tsteps[tsID].records[recID].position;
      size_t recsize = streamptr->tsteps[tsID].records[recID].size;

      fileSetPos(fileID, recpos, SEEK_SET);

      fileRead(fileID, gribbuffer, recsize);

      double missval = vlistInqVarMissval(vlistID, varID);

      int imiss;

      streamptr->tsteps[tsID].records[recID].zip = grbUnzipRecord(gribbuffer, &recsize);

      void *datap = NULL;
      if ( memtype == MEMTYPE_FLOAT )
        datap = (float*)data + levelID*gridsize;
      else
        datap = (double*)data + levelID*gridsize;

      grbDecode(filetype, memtype, gribbuffer, (int)recsize, datap, (size_t)gridsize,
                streamptr->unreduced, &imiss, missval, vlistID, varID);

      *nmiss += imiss;
    }

  fileSetPos(fileID, currentfilepos, SEEK_SET);
}


void grb_read_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, void *data, int *nmiss)
{
  int filetype = streamptr->filetype;

  void *gribbuffer = streamptr->record->buffer;

  int vlistID = streamptr->vlistID;
  int gridID   = vlistInqVarGrid(vlistID, varID);
  int gridsize = gridInqSize(gridID);
  int tsID = streamptr->curTsID;

  if ( CDI_Debug )
    Message("gridID = %d gridsize = %d", gridID, gridsize);

  int fileID = streamptr->fileID;

  off_t currentfilepos = fileGetPos(fileID);

  int    isub    = subtypeInqActiveIndex(streamptr->vars[varID].subtypeID);

  int    recID   = streamptr->vars[varID].recordTable[isub].recordID[levelID];
  off_t  recpos  = streamptr->tsteps[tsID].records[recID].position;
  size_t recsize = streamptr->tsteps[tsID].records[recID].size;

  if ( recsize == 0 )
    Error("Internal problem! Recordsize is zero for record %d at timestep %d",
	  recID+1, tsID+1);

  fileSetPos(fileID, recpos, SEEK_SET);
  fileRead(fileID, gribbuffer, recsize);

  streamptr->tsteps[tsID].records[recID].zip = grbUnzipRecord(gribbuffer, &recsize);

  double missval = vlistInqVarMissval(vlistID, varID);
  grbDecode(filetype, memtype, gribbuffer, (int)recsize, data, (size_t)gridsize, streamptr->unreduced, nmiss, missval, vlistID, varID);

  fileSetPos(fileID, currentfilepos, SEEK_SET);
}

#endif
