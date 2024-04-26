#ifdef  HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdf_int.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"


void recordInitEntry(record_t *record)
{
  record->position = CDI_UNDEFID;
  record->size     = 0;
  record->param    = 0;
  record->ilevel   = CDI_UNDEFID;
  record->used     = false;
  record->varID    = CDI_UNDEFID;
  record->levelID  = CDI_UNDEFID;
  memset(record->varname, 0, sizeof(record->varname));
  memset(&record->tiles, 0, sizeof(record->tiles));
}


int recordNewEntry(stream_t *streamptr, int tsID)
{
  size_t recordID = 0;
  size_t recordSize = (size_t)streamptr->tsteps[tsID].recordSize;
  record_t *records = streamptr->tsteps[tsID].records;
  /*
    Look for a free slot in record.
    (Create the table the first time through).
  */
  if ( ! recordSize )
    {
      recordSize = 1;   /*  <<<<----  */
      records = (record_t *) Malloc(recordSize * sizeof (record_t));

      for ( size_t i = 0; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }
  else
    {
      while ( recordID < recordSize
              && records[recordID].used != CDI_UNDEFID )
        ++recordID;
    }
  /*
    If the table overflows, double its size.
  */
  if ( recordID == recordSize )
    {
      if (recordSize <= INT_MAX / 2)
        recordSize *= 2;
      else if (recordSize < INT_MAX)
        recordSize = INT_MAX;
      else
        Error("Cannot handle this many records!\n");
      records = (record_t *) Realloc(records,
                                     recordSize * sizeof (record_t));

      for ( size_t i = recordID; i < recordSize; i++ )
	records[i].used = CDI_UNDEFID;
    }

  recordInitEntry(&records[recordID]);

  records[recordID].used = 1;

  streamptr->tsteps[tsID].recordSize = (int)recordSize;
  streamptr->tsteps[tsID].records    = records;

  return (int)recordID;
}

static
void cdiInitRecord(stream_t *streamptr)
{
  Record *record = (Record *) Malloc(sizeof(Record));
  streamptr->record = record;

  record->param      = 0;
  record->level      = 0;
  record->date       = 0;
  record->time       = 0;
  record->gridID     = 0;
  record->buffer     = NULL;
  record->buffersize = 0;
  record->position   = 0;
  record->varID      = 0;
  record->levelID    = CDI_UNDEFID;
}


void streamInqRecord(int streamID, int *varID, int *levelID)
{
  check_parg(varID);
  check_parg(levelID);

  stream_t *streamptr = stream_to_pointer(streamID);

  cdiDefAccesstype(streamID, TYPE_REC);

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  int tsID   = streamptr->curTsID;
  int rindex = streamptr->tsteps[tsID].curRecID + 1;

  if ( rindex >= streamptr->tsteps[tsID].nrecs )
    Error("record %d not available at timestep %d", rindex+1, tsID+1);

  int recID  = streamptr->tsteps[tsID].recIDs[rindex];

  if ( recID == -1 || recID >= streamptr->tsteps[tsID].nallrecs )
    Error("Internal problem! tsID = %d recID = %d", tsID, recID);

  *varID   = streamptr->tsteps[tsID].records[recID].varID;
  int lindex = streamptr->tsteps[tsID].records[recID].levelID;

  int isub = subtypeInqActiveIndex(streamptr->vars[*varID].subtypeID);
  *levelID = streamptr->vars[*varID].recordTable[isub].lindex[lindex];

  if ( CDI_Debug )
    Message("tsID = %d, recID = %d, varID = %d, levelID = %d", tsID, recID, *varID, *levelID);

  streamptr->curTsID = tsID;
  streamptr->tsteps[tsID].curRecID = rindex;
}

/*
@Function  streamDefRecord
@Title     Define the next record

@Prototype void streamDefRecord(int streamID, int varID, int levelID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.

@Description
The function streamDefRecord defines the meta-data of the next record.
@EndFunction
*/
void streamDefRecord(int streamID, int varID, int levelID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  int tsID = streamptr->curTsID;

  if ( tsID == CDI_UNDEFID )
    {
      tsID++;
      streamDefTimestep(streamID, tsID);
    }

  if ( ! streamptr->record ) cdiInitRecord(streamptr);

  int vlistID = streamptr->vlistID;
  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int param   = vlistInqVarParam(vlistID, varID);
  int level   = (int)(zaxisInqLevel(zaxisID, levelID));

  Record *record = streamptr->record;
  record->varID    = varID;
  record->levelID  = levelID;
  record->param    = param;
  record->level    = level;
  record->date     = streamptr->tsteps[tsID].taxis.vdate;
  record->time     = streamptr->tsteps[tsID].taxis.vtime;
  record->gridID   = gridID;
  record->prec     = vlistInqVarDatatype(vlistID, varID);

  switch (streamptr->filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      grbDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      srvDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      extDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      iegDefRecord(streamptr);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
      cdfDefRecord(streamptr);
      break;
#endif
    default:
      Error("%s support not compiled in!", strfiletype(streamptr->filetype));
      break;
    }
}


void streamCopyRecord(int streamID2, int streamID1)
{
  stream_t *streamptr1 = stream_to_pointer(streamID1),
    *streamptr2 = stream_to_pointer(streamID2);
  int filetype1 = streamptr1->filetype,
    filetype2 = streamptr2->filetype,
    filetype  = CDI_FILETYPE_UNDEF;

  if ( filetype1 == filetype2 ) filetype = filetype2;
  else
    {
      switch (filetype1)
        {
        case CDI_FILETYPE_NC:
        case CDI_FILETYPE_NC2:
        case CDI_FILETYPE_NC4:
        case CDI_FILETYPE_NC4C:
        case CDI_FILETYPE_NC5:
          switch (filetype2)
            {
            case CDI_FILETYPE_NC:
            case CDI_FILETYPE_NC2:
            case CDI_FILETYPE_NC4:
            case CDI_FILETYPE_NC4C:
            case CDI_FILETYPE_NC5:
              // Warning("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));
              filetype = filetype2;
              break;
            }
          break;
        }
    }

  if ( filetype == CDI_FILETYPE_UNDEF )
    Error("Streams have different file types (%s -> %s)!", strfiletype(filetype1), strfiletype(filetype2));

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      grbCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
      srvCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
      extCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
      iegCopyRecord(streamptr2, streamptr1);
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NC5:
      cdfCopyRecord(streamptr2, streamptr1);
      break;
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}


void cdi_create_records(stream_t *streamptr, int tsID)
{
  unsigned nrecords, maxrecords;

  tsteps_t *sourceTstep = streamptr->tsteps;
  tsteps_t *destTstep = sourceTstep + tsID;

  if ( destTstep->records ) return;

  int vlistID = streamptr->vlistID;

  if ( tsID == 0 )
    {
      maxrecords = 0;
      int nvars = streamptr->nvars;
      for ( int varID = 0; varID < nvars; varID++)
        for (int isub=0; isub<streamptr->vars[varID].subtypeSize; isub++)
          maxrecords += (unsigned)streamptr->vars[varID].recordTable[isub].nlevs;
    }
  else
    {
      maxrecords = (unsigned)sourceTstep->recordSize;
    }

  if ( tsID == 0 )
    {
      nrecords = maxrecords;
    }
  else if ( tsID == 1 )
    {
      nrecords = 0;
      maxrecords = (unsigned)sourceTstep->recordSize;
      for ( unsigned recID = 0; recID < maxrecords; recID++ )
	{
	  int varID = sourceTstep->records[recID].varID;
	  nrecords += (varID == CDI_UNDEFID /* varID = CDI_UNDEFID for write mode !!! */
                       || vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT);
          //    printf("varID nrecords %d %d %d \n", varID, nrecords, vlistInqVarTsteptype(vlistID, varID));
	}
    }
  else
    {
      nrecords = (unsigned)streamptr->tsteps[1].nallrecs;
    }
  //  printf("tsID, nrecords %d %d\n", tsID, nrecords);

  record_t *records = NULL;
  if ( maxrecords > 0 ) records = (record_t *) Malloc(maxrecords*sizeof(record_t));

  destTstep->records    = records;
  destTstep->recordSize = (int)maxrecords;
  destTstep->nallrecs   = (int)nrecords;

  if ( tsID == 0 )
    {
      for ( unsigned recID = 0; recID < maxrecords; recID++ )
        recordInitEntry(&destTstep->records[recID]);
    }
  else
    {
      memcpy(destTstep->records, sourceTstep->records, (size_t)maxrecords*sizeof(record_t));

      for ( unsigned recID = 0; recID < maxrecords; recID++ )
	{
          record_t *curRecord = &sourceTstep->records[recID];
          destTstep->records[recID].used = curRecord->used;
          if ( curRecord->used != CDI_UNDEFID && curRecord->varID != -1 ) /* curRecord->varID = -1 for write mode !!! */
            {
              if ( vlistInqVarTimetype(vlistID, curRecord->varID) != TIME_CONSTANT )
                {
                  destTstep->records[recID].position = CDI_UNDEFID;
                  destTstep->records[recID].size     = 0;
                  destTstep->records[recID].used     = false;
                }
            }
	}
    }
}

#include "file.h"

void streamFCopyRecord(stream_t *streamptr2, stream_t *streamptr1, const char *container_name)
{
  int fileID1 = streamptr1->fileID;
  int fileID2 = streamptr2->fileID;

  int tsID    = streamptr1->curTsID;
  int vrecID  = streamptr1->tsteps[tsID].curRecID;
  int recID   = streamptr1->tsteps[tsID].recIDs[vrecID];
  off_t recpos  = streamptr1->tsteps[tsID].records[recID].position;
  size_t recsize = streamptr1->tsteps[tsID].records[recID].size;

  if (fileSetPos(fileID1, recpos, SEEK_SET) != 0)
    Error("Cannot seek input file for %s record copy!", container_name);

  char *buffer = (char *) Malloc(recsize);

  if (fileRead(fileID1, buffer, recsize) != recsize)
    Error("Failed to read record from %s file for copying!", container_name);

  if (fileWrite(fileID2, buffer, recsize) != recsize)
    Error("Failed to write record to %s file when copying!", container_name);

  Free(buffer);
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
