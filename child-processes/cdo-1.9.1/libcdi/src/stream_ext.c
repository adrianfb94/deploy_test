#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_ext.h"
#include "varscan.h"
#include "datetime.h"
#include "extra.h"
#include "vlist.h"
#include "exse.h"


#if defined (HAVE_LIBEXTRA)

typedef struct {
  int param;
  int level;
} extcompvar_t;

static
int extInqDatatype(int prec, int number)
{
  int datatype;

  if ( number == 2 )
    datatype = (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
  else
    datatype = (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;

  return datatype;
}

static
void extDefDatatype(int datatype, int *prec, int *number)
{

  if ( datatype != CDI_DATATYPE_FLT32 && datatype != CDI_DATATYPE_FLT64 &&
       datatype != CDI_DATATYPE_CPX32 && datatype != CDI_DATATYPE_CPX64 )
    datatype = CDI_DATATYPE_FLT32;

  *number = (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) ? 2 : 1;

  *prec = (datatype == CDI_DATATYPE_FLT64 || datatype == CDI_DATATYPE_CPX64) ? EXSE_DOUBLE_PRECISION : EXSE_SINGLE_PRECISION;
}

/* not used
int extInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int header[4];
  int vlistID;
  void *extp = streamptr->record->exsep;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  *varID   = -1;
  *levelID = -1;

  status = extRead(fileID, extp);
  if ( status != 0 ) return 0;

  extInqHeader(extp, header);

  icode  = header[1];
  ilevel = header[2];

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == CDI_UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);

  return 1;
}
*/

void extReadRecord(stream_t *streamptr, double *data, int *nmiss)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int tsID    = streamptr->curTsID;
  int vrecID  = streamptr->tsteps[tsID].curRecID;
  int recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID   = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *extp = streamptr->record->exsep;
  int status = extRead(fileID, extp);
  if ( status != 0 ) Error("Failed to read EXTRA record");

  int header[4];
  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  int gridID  = vlistInqVarGrid(vlistID, varID);
  int size    = gridInqSize(gridID);

  streamptr->numvals += size;

  *nmiss = 0;
  if ( vlistInqVarNumber(vlistID, varID) == CDI_REAL )
    {
      for ( int i = 0; i < size; i++ )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
  else
    {
      for ( int i = 0; i < 2*size; i+=2 )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
}


void extCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "EXTRA");
}


void extDefRecord(stream_t *streamptr)
{
  int pdis, pcat, pnum;
  cdiDecodeParam(streamptr->record->param, &pnum, &pcat, &pdis);

  int header[4];
  header[0] = streamptr->record->date;
  header[1] = pnum;
  header[2] = streamptr->record->level;
  int gridID = streamptr->record->gridID;
  header[3] = gridInqSize(gridID);

  extrec_t *extp = (extrec_t*) streamptr->record->exsep;
  extDefDatatype(streamptr->record->prec, &extp->prec, &extp->number);
  extDefHeader(extp, header);
}


void extWriteRecord(stream_t *streamptr, const double *data)
{
  int fileID = streamptr->fileID;
  void *extp = streamptr->record->exsep;

  extDefDataDP(extp, data);
  extWrite(fileID, extp);
}

static
void extAddRecord(stream_t *streamptr, int param, int level, int xysize,
		  size_t recsize, off_t position, int prec, int number)
{
  int vlistID = streamptr->vlistID;
  int tsID    = streamptr->curTsID;
  int recID   = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  record->size     = recsize;
  record->position = position;
  record->param    = param;
  record->ilevel   = level;

  grid_t *grid = (grid_t *)Malloc(sizeof (*grid));
  grid_init(grid);
  cdiGridTypeInit(grid, GRID_GENERIC, xysize);
  grid->x.size = xysize;
  grid->y.size = 0;
  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 0);
  int gridID = gridAdded.Id;
  if (!gridAdded.isNew) Free(grid);
  /*
  if ( level == 0 ) leveltype = ZAXIS_SURFACE;
  else              leveltype = ZAXIS_GENERIC;
  */
  int leveltype = ZAXIS_GENERIC;

  int varID, levelID = 0;
  varAddRecord(recID, param, gridID, leveltype, 0, level, 0, 0, 0,
               extInqDatatype(prec, number), &varID, &levelID, TSTEP_INSTANT, 0, 0, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);
  record->varID   = (short)varID;
  record->levelID = (short)levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d gridID = %d levelID = %d", varID, gridID, levelID);
}

static
void extScanTimestep1(stream_t *streamptr)
{
  int header[4];
  DateTime datetime0 = { LONG_MIN, LONG_MIN };
  int varID;
  off_t recpos = 0;
  long recsize;
  int recID;
  extcompvar_t compVar, compVar0;
  extrec_t *extp = (extrec_t*) streamptr->record->exsep;

  streamptr->curTsID = 0;

  int tsID  = tstepsNewEntry(streamptr);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  int fileID = streamptr->fileID;

  int nrecs = 0;
  while ( true )
    {
      recpos = fileGetPos(fileID);
      int status = extRead(fileID, extp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 1;
	  break;
	}
      recsize = fileGetPos(fileID) - recpos;

      extInqHeader(extp, header);

      int vdate   = header[0];
      int vtime   = 0;
      int rcode   = header[1];
      int rlevel  = header[2];
      int rxysize = header[3];

      int param = cdiEncodeParam(rcode, 255, 255);

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}
      else
	{
	  compVar.param = param;
          compVar.level = rlevel;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      compVar0.param  = streamptr->tsteps[0].records[recID].param;
	      compVar0.level = streamptr->tsteps[0].records[recID].ilevel;

	      if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) == 0 ) break;
	    }
	  if ( recID < nrecs ) break;
	  DateTime datetime = { .date = vdate, .time = vtime};
	  if ( datetimeCmp(datetime, datetime0) )
	    Warning("Inconsistent verification time for code %d level %d", rcode, rlevel);
	}

      nrecs++;

      if ( CDI_Debug )
	Message("%4d%8d%4d%8d%8d%6d", nrecs, (int)recpos, rcode, rlevel, vdate, vtime);

      extAddRecord(streamptr, param, rlevel, rxysize, (size_t)recsize, recpos, extp->prec, extp->number);
    }

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
        (record_t *) Realloc(streamptr->tsteps[0].records, (size_t)nrecords * sizeof (record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) Malloc((size_t)nrecords * sizeof (int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
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
	  for ( varID = 0; varID < streamptr->nvars; varID++ )
            vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
    }
}

static
int extScanTimestep2(stream_t *streamptr)
{
  int header[4];
  int varID;
  off_t recpos = 0;
  extcompvar_t compVar, compVar0;
  void *extp = streamptr->record->exsep;

  streamptr->curTsID = 1;

  int fileID  = streamptr->fileID;
  int vlistID = streamptr->vlistID;

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
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position =
	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =
	streamptr->tsteps[0].records[recID].size;
    }

  for ( int rindex = 0; rindex <= nrecords; rindex++ )
    {
      recpos = fileGetPos(fileID);
      int status = extRead(fileID, extp);
      if ( status != 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

      extInqHeader(extp, header);

      int vdate  = header[0];
      int vtime  = 0;
      int rcode  = header[1];
      int rlevel = header[2];
      // rxysize = header[3];

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
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) == 0 )
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

      if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) != 0 )
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
	  varID = streamptr->tsteps[tsID].records[recID].varID;
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


int extInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;

  extScanTimestep1(streamptr);

  int status = 0;
  if ( streamptr->ntsteps == -1 ) status = extScanTimestep2(streamptr);

  int fileID = streamptr->fileID;
  fileSetPos(fileID, 0, SEEK_SET);

  return status;
}

static
long extScanTimestep(stream_t *streamptr)
{
  int header[4];
  off_t recpos = 0;
  int recID;
  int nrecs = 0;
  extcompvar_t compVar, compVar0;
  void *extp = streamptr->record->exsep;
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
	  int status = extRead(fileID, extp);
	  if ( status != 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  size_t recsize = (size_t)(fileGetPos(fileID) - recpos);

	  extInqHeader(extp, header);

	  int vdate  = header[0];
	  int vtime  = 0;
	  int rcode  = header[1];
	  int rlevel = header[2];
	  // rxysize = header[3];

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

	  compVar.param  = param;
          compVar.level  = rlevel;
	  compVar0.param = streamptr->tsteps[tsID].records[recID].param;
	  compVar0.level = streamptr->tsteps[tsID].records[recID].ilevel;

	  if ( memcmp(&compVar0, &compVar, sizeof(extcompvar_t)) != 0 )
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


int extInqTimestep(stream_t *streamptr, int tsID)
{
  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  long ntsteps = CDI_UNDEFID;
  while ( ( tsID + 1 ) > streamptr->rtsteps && ntsteps == CDI_UNDEFID )
    ntsteps = extScanTimestep(streamptr);

  int nrecs = 0;
  if ( !(tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID) )
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}


void extReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  void *extp = streamptr->record->exsep;

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
  extRead(fileID, extp);
  int header[4];
  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *nmiss = 0;
  if ( vlistInqVarNumber(vlistID, varID) == CDI_REAL )
    {
      for ( int i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
  else
    {
      for ( int i = 0; i < 2*gridsize; i+=2 )
	if ( DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float)missval) )
	  {
	    data[i] = missval;
	    (*nmiss)++;
	  }
    }
}


void extReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = (size_t) gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs    = (size_t) streamptr->vars[varID].recordTable[0].nlevs;

  for ( size_t levID = 0; levID < nlevs; levID++)
    extReadVarSliceDP(streamptr, varID, (int)levID, &data[levID*gridsize], nmiss);
}


void extWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  int vlistID  = streamptr->vlistID;
  int fileID   = streamptr->fileID;
  int tsID     = streamptr->curTsID;

  int pdis, pcat, pnum;
  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  int header[4];
  header[0] = streamptr->tsteps[tsID].taxis.vdate;
  header[1] = pnum;
  header[2] = (int)(zaxisInqLevel(vlistInqVarZaxis(vlistID, varID), levID));
  header[3] = gridInqSize(vlistInqVarGrid(vlistID, varID));

  extrec_t *extp = (extrec_t*) streamptr->record->exsep;
  extDefDatatype(vlistInqVarDatatype(vlistID, varID), &extp->prec, &extp->number);
  extDefHeader(extp, header);

  extDefDataDP(extp, data);
  extWrite(fileID, extp);
}


void extWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = (size_t) gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs    = (size_t) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));

  for ( size_t levID = 0;  levID < nlevs; levID++ )
    extWriteVarSliceDP(streamptr, varID, (int)levID, &data[levID*gridsize]);
}

#endif /* HAVE_LIBEXTRA */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
