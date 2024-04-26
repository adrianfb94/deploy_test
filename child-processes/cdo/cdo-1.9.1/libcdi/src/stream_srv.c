#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "varscan.h"
#include "datetime.h"
#include "service.h"
#include "stream_srv.h"
#include "vlist.h"
#include "exse.h"


#if defined (HAVE_LIBSERVICE)

typedef struct {
  int param;
  int level;
} srvcompvar_t;


static
int srvInqDatatype(int prec)
{
  return (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
}

static
int srvDefDatatype(int datatype)
{
  if ( datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64 )
    Error("CDI/SERVICE library does not support complex numbers!");

  if ( datatype != CDI_DATATYPE_FLT32 && datatype != CDI_DATATYPE_FLT64 )
    datatype = CDI_DATATYPE_FLT32;

  return (datatype == CDI_DATATYPE_FLT64) ? EXSE_DOUBLE_PRECISION : EXSE_SINGLE_PRECISION;
}

/* not used
int srvInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int header[8];
  int vlistID;
  void *srvp = streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  *varID   = -1;
  *levelID = -1;

  status = srvRead(fileID, srvp);
  if ( status != 0 ) return (0);

  srvInqHeader(srvp, header);

  icode  = header[0];
  ilevel = header[1];

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == CDI_UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);

  return 1;
}
*/

void srvReadRecord(stream_t *streamptr, double *data, int *nmiss)
{
  int vlistID  = streamptr->vlistID;
  int fileID   = streamptr->fileID;
  int tsID     = streamptr->curTsID;
  int vrecID   = streamptr->tsteps[tsID].curRecID;
  int recID    = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID    = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *srvp = streamptr->record->exsep;
  int status = srvRead(fileID, srvp);
  if ( status != 0 )
    Error("Failed to read record from SRV file");

  int header[8];
  srvInqHeader(srvp, header);
  srvInqDataDP(srvp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  int gridID  = vlistInqVarGrid(vlistID, varID);
  int size    = gridInqSize(gridID);

  streamptr->numvals += size;

  *nmiss = 0;
  for ( int i = 0; i < size; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void srvCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "SRV");
}


void srvDefRecord(stream_t *streamptr)
{
  Record *record = streamptr->record;
  srvrec_t *srvp = (srvrec_t*) record->exsep;

  int pdis, pcat, pnum;
  cdiDecodeParam(record->param, &pnum, &pcat, &pdis);

  int header[8];
  header[0] = pnum;
  header[1] = record->level;
  header[2] = record->date;
  header[3] = record->time;

  int gridID = record->gridID;
  int xsize = gridInqXsize(gridID),
      ysize = gridInqYsize(gridID);
  if ( xsize == 0 || ysize == 0 )
    {
      xsize = gridInqSize(gridID);
      ysize = 1;
    }
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ysize = 1;
  if ( gridInqSize(gridID) != xsize*ysize )
    Error("Internal problem with gridsize!");

  header[4] = xsize;
  header[5] = ysize;
  header[6] = 0;
  header[7] = 0;

  int datatype = record->prec;
  srvp->dprec = srvDefDatatype(datatype);

  srvDefHeader(srvp, header);
}


void srvWriteRecord(stream_t *streamptr, const double *data)
{
  int fileID = streamptr->fileID;
  void *srvp = streamptr->record->exsep;

  srvDefDataDP(srvp, data);
  srvWrite(fileID, srvp);
}

static
void srv_add_record(stream_t *streamptr, int param, int level, int xsize, int ysize,
                    size_t recsize, off_t position, int prec)
{
  int vlistID = streamptr->vlistID;
  int tsID    = streamptr->curTsID;
  int recID   = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  record->size     = recsize;
  record->position = position;
  record->param    = param;
  record->ilevel   = level;

  grid_t *grid = (grid_t*) Malloc(sizeof(*grid));
  grid_init(grid);
  cdiGridTypeInit(grid, GRID_GENERIC, xsize*ysize);
  grid->x.size = xsize;
  grid->y.size = ysize;
  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 0);
  int gridID = gridAdded.Id;
  if (!gridAdded.isNew) Free(grid);
  /*
  if ( level == 0 ) leveltype = ZAXIS_SURFACE;
  else              leveltype = ZAXIS_GENERIC;
  */
  int leveltype = ZAXIS_GENERIC;

  int datatype = srvInqDatatype(prec);

  int levelID = 0;
  int varID;
  varAddRecord(recID, param, gridID, leveltype, 0, level, 0, 0, 0,
	       datatype, &varID, &levelID, TSTEP_INSTANT, 0, 0, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);

  xassert(varID <= SHRT_MAX && levelID <= SHRT_MAX);
  record->varID   = (short)varID;
  record->levelID = (short)levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d gridID = %d levelID = %d", varID, gridID, levelID);
}

static
void srvScanTimestep1(stream_t *streamptr)
{
  DateTime datetime0 = { LONG_MIN, LONG_MIN };
  off_t recpos;
  taxis_t *taxis;
  srvrec_t *srvp = (srvrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 0;

  {
    int tsID  = tstepsNewEntry(streamptr);
    if ( tsID != 0 )
      Error("Internal problem! tstepsNewEntry returns %d", tsID);
    taxis = &streamptr->tsteps[tsID].taxis;
  }

  int fileID = streamptr->fileID;

  int nrecs = 0;
  while ( true )
    {
      int header[8];
      recpos = fileGetPos(fileID);
      int status = srvRead(fileID, srvp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

      srvInqHeader(srvp, header);

      int prec   = srvp->dprec;
      int rcode  = header[0];
      int rlevel = header[1];
      int vdate  = header[2];
      int vtime  = header[3];
      int rxsize = header[4];
      int rysize = header[5];

      int param = cdiEncodeParam(rcode, 255, 255);

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}
      else
	{
	  for ( int recID = 0; recID < nrecs; recID++ )
            if (    streamptr->tsteps[0].records[recID].param  == param
                 && streamptr->tsteps[0].records[recID].ilevel == rlevel )
              goto tstepScanLoopFinished;
	  DateTime datetime = { .date = vdate, .time = vtime };
	  if ( datetimeCmp(datetime, datetime0) )
	    Warning("Inconsistent verification time for code %d level %d", rcode, rlevel);
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", nrecs, (int)recpos, rcode, rlevel, vdate, vtime);

      srv_add_record(streamptr, param, rlevel, rxsize, rysize, recsize, recpos, prec);
    }

  tstepScanLoopFinished:
  streamptr->rtsteps = 1;

  cdi_generate_vars(streamptr);

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type  = TAXIS_ABSOLUTE;
  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  vlist_check_contents(vlistID);

  int nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
	(record_t *) Realloc(streamptr->tsteps[0].records,
                             (size_t)nrecords * sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) Malloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( int recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  if ( streamptr->ntsteps == -1 )
    {
      int tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = true;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( int varID = 0; varID < streamptr->nvars; varID++ )
            vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
    }
}

static
int srvScanTimestep2(stream_t *streamptr)
{
  int header[8];
  off_t recpos = 0;
  srvcompvar_t compVar, compVar0;
  void *srvp = streamptr->record->exsep;

  streamptr->curTsID = 1;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  int nrecords = streamptr->tsteps[0].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) Malloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[1].nrecs = 0;
  for ( int recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( int recID = 0; recID < nrecords; recID++ )
    {
      streamptr->tsteps[tsID].records[recID].position =
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =
	streamptr->tsteps[0].records[recID].size;
    }

  for ( int rindex = 0; rindex <= nrecords; rindex++ )
    {
      recpos = fileGetPos(fileID);
      int status = srvRead(fileID, srvp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

      srvInqHeader(srvp, header);

      int rcode  = header[0];
      int rlevel = header[1];
      int vdate  = header[2];
      int vtime  = header[3];

      int param = cdiEncodeParam(rcode, 255, 255);

      if ( rindex == 0 )
	{
	  taxis->type  = TAXIS_ABSOLUTE;
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;
	}

      compVar.param = param;
      compVar.level = rlevel;
      bool nextstep = false;
      int recID;
      for ( recID = 0; recID < nrecords; recID++ )
	{
	  compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(srvcompvar_t)) == 0 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  nextstep = true;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = true;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	      break;
	    }
	}
      if ( recID == nrecords )
	{
	  Warning("Code %d level %d not found at timestep %d", rcode, rlevel, tsID+1);
	  return CDI_EUFSTRUCT;
	}

      if ( nextstep ) break;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", rindex+1, (int)recpos, rcode, rlevel, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      compVar0.param  = streamptr->tsteps[tsID].records[recID].param;
      compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

      if ( memcmp(&compVar0, &compVar, sizeof(srvcompvar_t)) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	  return CDI_EUFSTRUCT;
	}

      streamptr->tsteps[1].records[recID].position = recpos;
    }

  int nrecs = 0;
  for ( int recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  int varID = streamptr->tsteps[tsID].records[recID].varID;
          vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = true;
      streamptr->tsteps[tsID].position = recpos;
    }

  return 0;
}


int srvInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;

  srvScanTimestep1(streamptr);

  int status = 0;
  if ( streamptr->ntsteps == -1 ) status = srvScanTimestep2(streamptr);

  int fileID = streamptr->fileID;
  fileSetPos(fileID, 0, SEEK_SET);

  return status;
}

static
long srvScanTimestep(stream_t *streamptr)
{
  int header[8];
  /* int rxsize = 0, rysize = 0; */
  off_t recpos = 0;
  int recID;
  int nrecs = 0;
  void *srvp = streamptr->record->exsep;
  /*
  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }
  */

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) Malloc((size_t)nrecs * sizeof (int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for ( int rindex = 0; rindex <= nrecs; rindex++ )
	{
	  recpos = fileGetPos(fileID);
	  int status = srvRead(fileID, srvp);
	  if ( status != 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

	  srvInqHeader(srvp, header);

	  int rcode  = header[0];
	  int rlevel = header[1];
	  int vdate  = header[2];
	  int vtime  = header[3];
          /* rxsize = header[4]; */
          /* rysize = header[5]; */

	  int param = cdiEncodeParam(rcode, 255, 255);

	  // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
	  if ( rindex == nrecs ) continue;
	  recID = streamptr->tsteps[tsID].recIDs[rindex];

	  if ( rindex == 0 )
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;
	    }

          if (    param  != streamptr->tsteps[tsID].records[recID].param
               || rlevel != streamptr->tsteps[tsID].records[recID].ilevel )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, rlevel);
	      Error("Invalid, unsupported or inconsistent record structure!");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d%8d%4d%8d%8d%6d", rindex, (int)recpos, rcode, rlevel, vdate, vtime);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = true;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return streamptr->ntsteps;
}


int srvInqTimestep(stream_t *streamptr, int tsID)
{
  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  long ntsteps = CDI_UNDEFID;
  while ( ( tsID + 1 ) > streamptr->rtsteps && ntsteps == CDI_UNDEFID )
    ntsteps = srvScanTimestep(streamptr);

  int nrecs = 0;
  if ( !(tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID) )
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}


void srvReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  void *srvp = streamptr->record->exsep;

  int vlistID  = streamptr->vlistID;
  int fileID   = streamptr->fileID;
  /* NOTE: tiles are not supported here! */
  double missval = vlistInqVarMissval(vlistID, varID);
  int gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  int tsid     = streamptr->curTsID;

  off_t currentfilepos = fileGetPos(fileID);

  /* NOTE: tiles are not supported here! */
  int recID = streamptr->vars[varID].recordTable[0].recordID[levID];
  off_t recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  if ( srvRead(fileID, srvp) < 0 ) abort();
  int header[8];
  srvInqHeader(srvp, header);
  srvInqDataDP(srvp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  for ( int i = 0; i < gridsize; i++ )
    if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
      {
	data[i] = missval;
	(*nmiss)++;
      }
}


void srvReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = (size_t) gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs    = (size_t) streamptr->vars[varID].recordTable[0].nlevs;

  for ( size_t levID = 0; levID < nlevs; levID++)
    srvReadVarSliceDP(streamptr, varID, (int)levID, &data[levID*gridsize], nmiss);
}


void srvWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  int vlistID  = streamptr->vlistID;
  int fileID   = streamptr->fileID;
  int tsID     = streamptr->curTsID;
  int gridID   = vlistInqVarGrid(vlistID, varID);

  int pdis, pcat, pnum;
  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  int header[8];
  header[0] = pnum;
  header[1] = (int)(zaxisInqLevel(vlistInqVarZaxis(vlistID, varID), levID));
  header[2] = streamptr->tsteps[tsID].taxis.vdate;
  header[3] = streamptr->tsteps[tsID].taxis.vtime;

  int xsize = gridInqXsize(gridID);
  int ysize = gridInqYsize(gridID);
  if ( xsize == 0 || ysize == 0 )
    {
      xsize = gridInqSize(gridID);
      ysize = 1;
    }
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ysize = 1;
  if ( gridInqSize(gridID) != xsize*ysize )
    Error("Internal problem with gridsize!");

  header[4] = xsize;
  header[5] = ysize;
  header[6] = 0;
  header[7] = 0;

  int datatype = vlistInqVarDatatype(vlistID, varID);

  srvrec_t *srvp = (srvrec_t*) streamptr->record->exsep;
  srvp->dprec = srvDefDatatype(datatype);

  srvDefHeader(srvp, header);
  srvDefDataDP(srvp, data);
  srvWrite(fileID, srvp);
}


void srvWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = (size_t) gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs    = (size_t) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));

  for ( size_t levID = 0; levID < nlevs; levID++ )
    srvWriteVarSliceDP(streamptr, varID, (int)levID, &data[levID*gridsize]);
}

#endif /* HAVE_LIBSERVICE */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
