#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include "dmemory.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "stream_cdf.h"
#include "cdf_int.h"
#include "varscan.h"
#include "vlist.h"
#include "zaxis.h"


#define  POSITIVE_UP    1
#define  POSITIVE_DOWN  2


static const char bndsName[] = "bnds";


void cdfCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  int vlistID1 = streamptr1->vlistID;
  int tsID     = streamptr1->curTsID;
  int vrecID   = streamptr1->tsteps[tsID].curRecID;
  int recID    = streamptr1->tsteps[tsID].recIDs[vrecID];
  int ivarID   = streamptr1->tsteps[tsID].records[recID].varID;
  int gridID   = vlistInqVarGrid(vlistID1, ivarID);
  int datasize = gridInqSize(gridID);
  int datatype = vlistInqVarDatatype(vlistID1, ivarID);
  int memtype  = datatype != CDI_DATATYPE_FLT32 ? MEMTYPE_DOUBLE : MEMTYPE_FLOAT;

  void *data = Malloc((size_t)datasize
             * (memtype == MEMTYPE_DOUBLE ? sizeof(double) : sizeof(float)));

  int nmiss;
  cdf_read_record(streamptr1, memtype, data, &nmiss);
  cdf_write_record(streamptr2, memtype, data, nmiss);

  Free(data);
}


void cdfDefRecord(stream_t *streamptr)
{
  (void)streamptr;
}

static
void cdfDefTimeValue(stream_t *streamptr, int tsID)
{
  int fileID = streamptr->fileID;

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d", streamptr->self, fileID);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  double timevalue = cdiEncodeTimeval(taxis->vdate, taxis->vtime, &streamptr->tsteps[0].taxis);
  if ( CDI_Debug ) Message("tsID = %d  timevalue = %f", tsID, timevalue);

  int ncvarid = streamptr->basetime.ncvarid;
  size_t index = (size_t)tsID;
  cdf_put_var1_double(fileID, ncvarid, &index, &timevalue);

  if ( taxis->has_bounds )
    {
      ncvarid = streamptr->basetime.ncvarboundsid;
      if ( ncvarid == CDI_UNDEFID ) Error("Call to taxisWithBounds() missing!");

      timevalue = cdiEncodeTimeval(taxis->vdate_lb, taxis->vtime_lb, &streamptr->tsteps[0].taxis);
      size_t start[2], count[2];
      start[0] = (size_t)tsID; count[0] = 1; start[1] = 0; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);

      timevalue = cdiEncodeTimeval(taxis->vdate_ub, taxis->vtime_ub, &streamptr->tsteps[0].taxis);
      start[0] = (size_t)tsID; count[0] = 1; start[1] = 1; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);
    }

  ncvarid = streamptr->basetime.leadtimeid;
  if ( taxis->type == TAXIS_FORECAST && ncvarid != CDI_UNDEFID )
    {
      timevalue = taxis->fc_period;
      cdf_put_var1_double(fileID, ncvarid, &index, &timevalue);
    }
}

void cdfDefTimestep(stream_t *streamptr, int tsID)
{
  cdfDefTimeValue(streamptr, tsID);
}

static
void cdfDefComplex(stream_t *streamptr, int gridID, int gridindex)
{
  int dimID;
  ncgrid_t *ncgrid = streamptr->ncgrid;

  for ( int index = 0; index < gridindex; ++index )
    {
      if ( ncgrid[index].ncIDs[CDF_DIMID_X] != CDI_UNDEFID )
        {
          int gridID0 = ncgrid[index].gridID;
          int gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_SPECTRAL || gridtype0 == GRID_FOURIER )
            {
              dimID = ncgrid[index].ncIDs[CDF_DIMID_X];
              goto dimIDEstablished;
            }
        }
    }

  {
    static const char axisname[] = "nc2";
    size_t dimlen = 2;
    int fileID  = streamptr->fileID;

    if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
    cdf_def_dim(fileID, axisname, dimlen, &dimID);
    cdf_enddef(fileID);
    streamptr->ncmode = 2;
  }
  dimIDEstablished:
  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_X] = dimID;
}

struct idSearch
{
  int numNonMatching, foundID;
  size_t foundIdx;
};

static inline struct idSearch
cdfSearchIDBySize(size_t startIdx, size_t numIDs, const ncgrid_t ncgrid[numIDs],
                  int ncIDType, int searchType, int searchSize,
                  int (*typeInq)(int id), int (*sizeInq)(int id))
{
  int numNonMatching = 0,
    foundID = CDI_UNDEFID;
  size_t foundIdx = SIZE_MAX;
  for ( size_t index = startIdx; index < numIDs; index++ )
    {
      if ( ncgrid[index].ncIDs[ncIDType] != CDI_UNDEFID )
        {
          int id0 = ncgrid[index].gridID,
            id0Type = typeInq(id0);
          if ( id0Type == searchType )
            {
              int size0 = sizeInq(id0);
              if ( searchSize == size0 )
                {
                  foundID = ncgrid[index].ncIDs[ncIDType];
                  foundIdx = index;
                  break;
                }
              numNonMatching++;
            }
        }
    }
  return (struct idSearch){ .numNonMatching = numNonMatching,
      .foundID = foundID, .foundIdx = foundIdx };
}

static int
cdfGridInqHalfSize(int gridID)
{
  return gridInqSize(gridID)/2;
}


static void
cdfDefSPorFC(stream_t *streamptr, int gridID, int gridindex,
             char *restrict axisname, int gridRefType)
{
  ncgrid_t *ncgrid = streamptr->ncgrid;

  size_t dimlen = (size_t)(gridInqSize(gridID))/2;

  int iz;
  int dimID;
  {
    struct idSearch search
      = cdfSearchIDBySize(0, (size_t)gridindex, ncgrid, CDF_DIMID_Y,
                          gridRefType, (int)dimlen,
                          gridInqType, cdfGridInqHalfSize);
    dimID = search.foundID;
    iz = search.numNonMatching;
  }

  if ( dimID == CDI_UNDEFID )
    {
      int fileID  = streamptr->fileID;
      if ( iz == 0 ) axisname[3] = '\0';
      else           sprintf(&axisname[3], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_Y] = dimID;
}

static
void cdfDefSP(stream_t *streamptr, int gridID, int gridindex)
{
  /*
  char longname[] = "Spherical harmonic coefficient";
  */
  char axisname[5] = "nspX";
  cdfDefSPorFC(streamptr, gridID, gridindex, axisname, GRID_SPECTRAL);
}


static
void cdfDefFC(stream_t *streamptr, int gridID, int gridindex)
{
  char axisname[5] = "nfcX";
  cdfDefSPorFC(streamptr, gridID, gridindex, axisname, GRID_FOURIER);
}

static const struct cdfDefGridAxisInqs {
  int (*axisSize)(int gridID);
  int (*axisDimname)(int cdiID, int key, int size, char *mesg);
  int (*axisName)(int cdiID, int key, int size, char *mesg);
  int (*axisLongname)(int cdiID, int key, int size, char *mesg);
  int (*axisUnits)(int cdiID, int key, int size, char *mesg);
  void (*axisStdname)(int cdiID, char *dimstdname);
  double (*axisVal)(int gridID, int index);
  const double *(*axisValsPtr)(int gridID);
  const double *(*axisBoundsPtr)(int gridID);
} gridInqsX = {
  .axisSize = gridInqXsize,
  .axisDimname = cdiGridInqKeyStr,
  .axisName = cdiGridInqKeyStr,
  .axisLongname = cdiGridInqKeyStr,
  .axisUnits = cdiGridInqKeyStr,
  .axisStdname = gridInqXstdname,
  .axisVal = gridInqXval,
  .axisValsPtr = gridInqXvalsPtr,
  .axisBoundsPtr = gridInqXboundsPtr,
}, gridInqsY = {
  .axisSize = gridInqYsize,
  .axisDimname = cdiGridInqKeyStr,
  .axisName = cdiGridInqKeyStr,
  .axisLongname = cdiGridInqKeyStr,
  .axisUnits = cdiGridInqKeyStr,
  .axisStdname = gridInqYstdname,
  .axisVal = gridInqYval,
  .axisValsPtr = gridInqYvalsPtr,
  .axisBoundsPtr = gridInqYboundsPtr,
}, gridInqsZ = {
  .axisLongname = cdiZaxisInqKeyStr,
  .axisUnits = cdiZaxisInqKeyStr,
  .axisStdname = zaxisInqStdname,
};

static
void cdfPutGridStdAtts(int fileID, int ncvarid, int gridID, int dimtype, const struct cdfDefGridAxisInqs *inqs)
{
  size_t len;

  char stdname[CDI_MAX_NAME];
  inqs->axisStdname(gridID, stdname);
  if ( (len = strlen(stdname)) )
    cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);

  char longname[CDI_MAX_NAME]; longname[0] = 0;
  int keyname = (dimtype == 'Z') ? CDI_KEY_LONGNAME : (dimtype == 'X') ? CDI_KEY_XLONGNAME : CDI_KEY_YLONGNAME;
  inqs->axisLongname(gridID, keyname, CDI_MAX_NAME, longname);
  if ( longname[0] && (len = strlen(longname)) )
    cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);

  char units[CDI_MAX_NAME]; units[0] = 0;
  keyname = (dimtype == 'Z') ? CDI_KEY_UNITS : (dimtype == 'X') ? CDI_KEY_XUNITS : CDI_KEY_YUNITS;
  inqs->axisUnits(gridID, keyname, CDI_MAX_NAME, units);
  if ( units[0] && (len = strlen(units)) )
    cdf_put_att_text(fileID, ncvarid, "units", len, units);
}

static void
cdfDefTrajLatLon(stream_t *streamptr, int gridID, int gridindex,
                 const struct cdfDefGridAxisInqs *inqs, int dimtype)
{
  nc_type xtype = (gridInqDatatype(gridID) == CDI_DATATYPE_FLT32) ? NC_FLOAT : NC_DOUBLE;
  ncgrid_t *ncgrid = streamptr->ncgrid;

  int dimlen = inqs->axisSize(gridID);
  if ( dimlen != 1 )
    Error("%c size isn't 1 for %s grid!", dimtype, gridNamePtr(gridInqType(gridID)));

  int ncvarid = ncgrid[gridindex].ncIDs[dimtype == 'X' ? CDF_DIMID_X : CDF_DIMID_Y];

  if ( ncvarid == CDI_UNDEFID )
    {
      int dimNcID = streamptr->basetime.ncvarid;
      int fileID  = streamptr->fileID;
      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      char axisname[CDI_MAX_NAME]; axisname[0] = 0;
      int keyname = (dimtype == 'X') ? CDI_KEY_XNAME : CDI_KEY_YNAME;
      inqs->axisName(gridID, keyname, CDI_MAX_NAME, axisname);
      cdf_def_var(fileID, axisname, xtype, 1, &dimNcID, &ncvarid);
      cdfPutGridStdAtts(fileID, ncvarid, gridID, dimtype, inqs);
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  ncgrid[gridindex].gridID = gridID;
  /* var ID for trajectory !!! */
  ncgrid[gridindex].ncIDs[dimtype == 'X' ? CDF_DIMID_X : CDF_DIMID_Y] = ncvarid;
}

static
void cdfDefTrajLon(stream_t *streamptr, int gridID, int gridindex)
{
  cdfDefTrajLatLon(streamptr, gridID, gridindex, &gridInqsX, 'X');
}


static
void cdfDefTrajLat(stream_t *streamptr, int gridID, int gridindex)
{
  cdfDefTrajLatLon(streamptr, gridID, gridindex, &gridInqsY, 'Y');
}

static
int checkDimName(int fileID, size_t dimlen, char *dimname)
{
  /* check whether the dimenion name is already defined with the same length */
  unsigned iz = 0;
  int dimid = CDI_UNDEFID;
  char name[CDI_MAX_NAME];

  size_t len = strlen(dimname);
  memcpy(name, dimname, len + 1);

  do
    {
      if ( iz ) sprintf(name + len, "_%u", iz+1);

      int dimid0, status = nc_inq_dimid(fileID, name, &dimid0);
      if ( status != NC_NOERR )
        break;
      size_t dimlen0;
      cdf_inq_dimlen(fileID, dimid0, &dimlen0);
      if ( dimlen0 == dimlen )
        {
          dimid = dimid0;
          break;
        }
      iz++;
    }
  while ( iz <= 99 );


  if ( iz ) sprintf(dimname + len, "_%u", iz+1);

  return dimid;
}

static
void checkGridName(char *axisname, int fileID)
{
  int ncdimid;
  char axisname2[CDI_MAX_NAME];

  /* check that the name is not already defined */
  unsigned iz = 0;

  size_t axisnameLen = strlen(axisname);
  memcpy(axisname2, axisname, axisnameLen + 1);
  do
    {
      if ( iz ) sprintf(axisname2 + axisnameLen, "_%u", iz+1);

      int status = nc_inq_varid(fileID, axisname2, &ncdimid);
      if ( status != NC_NOERR ) break;

      ++iz;
    }
  while ( iz <= 99 );

  if ( iz ) sprintf(axisname + axisnameLen, "_%u", iz+1);
}

static
int checkZaxisName(char *axisname, int fileID, int vlistID, int zaxisID, int nzaxis)
{
  char axisname2[CDI_MAX_NAME];

  /* check that the name is not already defined */
  unsigned iz = 0;

  size_t axisnameLen = strlen(axisname);
  memcpy(axisname2, axisname, axisnameLen + 1);
  do
    {
      if ( iz ) sprintf(axisname2 + axisnameLen, "_%u", iz+1);

      int ncdimid, status = nc_inq_varid(fileID, axisname2, &ncdimid);

      if ( status != NC_NOERR )
        {
          if ( iz )
            {
              /* check that the name does not exist for other zaxes */
              for ( int index = 0; index < nzaxis; index++ )
                {
                  int zaxisID0 = vlistZaxis(vlistID, index);
                  if ( zaxisID != zaxisID0 )
                    {
                      const char *axisname0 = zaxisInqNamePtr(zaxisID0);
                      if ( strcmp(axisname0, axisname2) == 0 ) goto nextSuffix;
                    }
                }
            }
          break;
        }
      nextSuffix:
      ++iz;
    }
  while (iz <= 99);


  if ( iz ) sprintf(axisname + axisnameLen, "_%u", iz+1);

  return (int)iz;
}

static void
cdfDefAxisCommon(stream_t *streamptr, int gridID, int gridindex, int ndims,
                 const struct cdfDefGridAxisInqs *gridAxisInq, int dimKey, char axisLetter,
                 void (*finishCyclicBounds)(double *pbounds, size_t dimlen, const double *pvals))
{
  int dimID = CDI_UNDEFID;
  int ncvarid = CDI_UNDEFID, ncbvarid = CDI_UNDEFID;
  int nvdimID = CDI_UNDEFID;
  int fileID  = streamptr->fileID;
  size_t dimlen = (size_t)gridAxisInq->axisSize(gridID);
  nc_type xtype = (nc_type)cdfDefDatatype(gridInqDatatype(gridID), streamptr->filetype);

  ncgrid_t *ncgrid = streamptr->ncgrid;

  const double *pvals = gridAxisInq->axisValsPtr(gridID);
  char dimname[CDI_MAX_NAME+3]; dimname[0] = 0;
  if ( ndims && pvals == NULL ) cdiGridInqKeyStr(gridID, dimKey, CDI_MAX_NAME, dimname);

  for ( int index = 0; index < gridindex; ++index )
    {
      int gridID0 = ncgrid[index].gridID;
      assert(gridID0 != CDI_UNDEFID);
      int gridtype0 = gridInqType(gridID0);
      if ( gridtype0 == GRID_GAUSSIAN    ||
           gridtype0 == GRID_LONLAT      ||
           gridtype0 == GRID_PROJECTION  ||
           gridtype0 == GRID_CURVILINEAR ||
           gridtype0 == GRID_GENERIC )
        {
          size_t dimlen0 = (size_t)gridAxisInq->axisSize(gridID0);
          char dimname0[CDI_MAX_NAME]; dimname0[0] = 0;
          if ( dimname[0] ) cdiGridInqKeyStr(gridID0, dimKey, CDI_MAX_NAME, dimname0);
          bool lname = dimname0[0] ? strcmp(dimname, dimname0) == 0 : true;
          if ( dimlen == dimlen0 && lname )
            {
              double (*inqVal)(int gridID, int index) = gridAxisInq->axisVal;
              if ( IS_EQUAL(inqVal(gridID0, 0), inqVal(gridID, 0)) &&
                   IS_EQUAL(inqVal(gridID0, (int)dimlen-1), inqVal(gridID, (int)dimlen-1)) )
                {
                  dimID = ncgrid[index].ncIDs[dimKey == CDI_KEY_XDIMNAME
                                              ? CDF_DIMID_X : CDF_DIMID_Y];
                  break;
                }
            }
        }
    }

  if ( dimID == CDI_UNDEFID )
    {
      char axisname[CDI_MAX_NAME]; axisname[0] = 0;
      int keyname = (axisLetter == 'X') ? CDI_KEY_XNAME : CDI_KEY_YNAME;
      gridAxisInq->axisName(gridID, keyname, CDI_MAX_NAME, axisname);
      if ( axisname[0] == 0 ) Error("axis name undefined!");
      size_t axisnameLen = strlen(axisname);

      /* enough to append _ plus up to 100 decimal and trailing \0 */
      char extendedAxisname[axisnameLen + 4 + 1];
      memcpy(extendedAxisname, axisname, axisnameLen + 1);
      checkGridName(extendedAxisname, fileID);
      size_t extendedAxisnameLen = axisnameLen + strlen(extendedAxisname + axisnameLen);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      if ( ndims )
        {
          if ( dimname[0] == 0 ) strcpy(dimname, extendedAxisname);
          dimID = checkDimName(fileID, dimlen, dimname);

          if ( dimID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, dimlen, &dimID);
        }

      bool gen_bounds = false;
      bool grid_is_cyclic = gridIsCircular(gridID) > 0;
      double *pbounds = NULL;
      if ( pvals )
        {
          cdf_def_var(fileID, extendedAxisname, xtype, ndims, &dimID, &ncvarid);

          cdfPutGridStdAtts(fileID, ncvarid, gridID, axisLetter, gridAxisInq);
          {
            char axisStr[2] = { axisLetter, '\0' };
            cdf_put_att_text(fileID, ncvarid, "axis", 1, axisStr);
          }

          size_t nvertex = gridInqNvertex(gridID);
          pbounds = (double *)gridAxisInq->axisBoundsPtr(gridID);

          if ( CDI_cmor_mode && grid_is_cyclic && !pbounds )
            {
              gen_bounds = true;
              nvertex = 2;
              pbounds = (double*) Malloc(2*dimlen*sizeof(double));
              for ( size_t i = 0; i < dimlen-1; ++i )
                {
                  pbounds[i*2+1]   = (pvals[i] + pvals[i+1])/2;
                  pbounds[(i+1)*2] = (pvals[i] + pvals[i+1])/2;
                }
              finishCyclicBounds(pbounds, dimlen, pvals);
            }
          if ( pbounds )
            {
              if ( nc_inq_dimid(fileID, bndsName, &nvdimID) != NC_NOERR )
                cdf_def_dim(fileID, bndsName, nvertex, &nvdimID);
            }
          if ( pbounds && nvdimID != CDI_UNDEFID )
            {
              char boundsname[extendedAxisnameLen + 1 + sizeof(bndsName)];
              memcpy(boundsname, axisname, extendedAxisnameLen);
              boundsname[extendedAxisnameLen] = '_';
              memcpy(boundsname + extendedAxisnameLen + 1, bndsName, sizeof bndsName);
              int dimIDs[2] = { dimID, nvdimID };
              cdf_def_var(fileID, boundsname, xtype, 2, dimIDs, &ncbvarid);
              cdf_put_att_text(fileID, ncvarid, "bounds", extendedAxisnameLen + sizeof (bndsName), boundsname);
            }
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncvarid  != CDI_UNDEFID ) cdf_put_var_double(fileID, ncvarid, pvals);
      if ( ncbvarid != CDI_UNDEFID ) cdf_put_var_double(fileID, ncbvarid, pbounds);
      if ( gen_bounds ) Free(pbounds);

      if ( ndims == 0 )
        ncgrid[gridindex].ncIDs[dimKey == CDI_KEY_XDIMNAME
                                ? CDF_VARID_X : CDF_VARID_Y] = ncvarid;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[dimKey == CDI_KEY_XDIMNAME
                          ? CDF_DIMID_X : CDF_DIMID_Y] = dimID;
}

static
void finishCyclicXBounds(double *pbounds, size_t dimlen, const double *pvals)
{
  pbounds[0] = (pvals[0] + pvals[dimlen-1]-360)*0.5;
  pbounds[2*dimlen-1] = (pvals[dimlen-1] + pvals[0]+360)*0.5;
}

static
void cdfDefXaxis(stream_t *streamptr, int gridID, int gridindex, int ndims)
{
  cdfDefAxisCommon(streamptr, gridID, gridindex, ndims, &gridInqsX,
                   CDI_KEY_XDIMNAME, 'X', finishCyclicXBounds);
}

static
void finishCyclicYBounds(double *pbounds, size_t dimlen, const double *pvals)
{
  pbounds[0] = copysign(90.0, pvals[0]);
  pbounds[2*dimlen-1] = copysign(90.0, pvals[dimlen-1]);
}

static
void cdfDefYaxis(stream_t *streamptr, int gridID, int gridindex, int ndims)
{
  cdfDefAxisCommon(streamptr, gridID, gridindex, ndims, &gridInqsY,
                   CDI_KEY_YDIMNAME, 'Y', finishCyclicYBounds);
}

static
void cdfGridCompress(int fileID, int ncvarid, int gridsize, int filetype, int comptype)
{
#if  defined  (HAVE_NETCDF4)
  if ( gridsize > 1 && comptype == CDI_COMPRESS_ZIP && (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C) )
    {
      cdf_def_var_chunking(fileID, ncvarid, NC_CHUNKED, NULL);
      cdfDefVarDeflate(fileID, ncvarid, 1);
    }
#endif
}

static
void cdfDefGridReference(stream_t *streamptr, int gridID)
{
  int fileID  = streamptr->fileID;
  int number = gridInqNumber(gridID);

  if ( number > 0 )
    {
      cdf_put_att_int(fileID, NC_GLOBAL, "number_of_grid_used", NC_INT, 1, &number);
    }

  const char *gridfile = gridInqReferencePtr(gridID);
  if ( gridfile && gridfile[0] != 0 )
    cdf_put_att_text(fileID, NC_GLOBAL, "grid_file_uri", strlen(gridfile), gridfile);
}

static
void cdfDefGridUUID(stream_t *streamptr, int gridID)
{
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];

  gridInqUUID(gridID, uuidOfHGrid);
  if ( !cdiUUIDIsNull(uuidOfHGrid) )
    {
      char uuidOfHGridStr[37];
      cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
      if ( uuidOfHGridStr[0] != 0 && strlen(uuidOfHGridStr) == 36 )
        {
          int fileID  = streamptr->fileID;
          //if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
          cdf_put_att_text(fileID, NC_GLOBAL, "uuidOfHGrid", 36, uuidOfHGridStr);
          //if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
        }
    }
}

struct cdfDefIrregularGridCommonIDs
{
  int xdimID, ydimID, ncxvarid, ncyvarid, ncavarid;
};

static struct cdfDefIrregularGridCommonIDs
cdfDefIrregularGridCommon(stream_t *streamptr, int gridID,
                          size_t xdimlen, size_t ydimlen,
                          int ndims, const char *xdimname_default,
                          size_t nvertex, const char *vdimname_default,
                          bool setVdimname)
{
  nc_type xtype = (nc_type)cdfDefDatatype(gridInqDatatype(gridID), streamptr->filetype);
  int xdimID = CDI_UNDEFID;
  int ydimID = CDI_UNDEFID;
  int ncxvarid = CDI_UNDEFID, ncyvarid = CDI_UNDEFID, ncavarid = CDI_UNDEFID;
  int ncbxvarid = CDI_UNDEFID, ncbyvarid = CDI_UNDEFID;
  int fileID  = streamptr->fileID;
  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  {
    char xdimname[CDI_MAX_NAME+3];
    xdimname[0] = 0;
    cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, xdimname);
    if ( xdimname[0] == 0 ) strcpy(xdimname, xdimname_default);
    xdimID = checkDimName(fileID, xdimlen, xdimname);
    if ( xdimID == CDI_UNDEFID ) cdf_def_dim(fileID, xdimname, xdimlen, &xdimID);
  }

  if ( ndims == 3 )
    {
      char ydimname[CDI_MAX_NAME+3];
      ydimname[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, ydimname);
      if ( ydimname[0] == 0 ) { ydimname[0] = 'y'; ydimname[1] = 0; }
      ydimID = checkDimName(fileID, ydimlen, ydimname);
      if ( ydimID == CDI_UNDEFID ) cdf_def_dim(fileID, ydimname, ydimlen, &ydimID);
    }

  int nvdimID = CDI_UNDEFID;
  int dimIDs[3];
  dimIDs[ndims-1] = CDI_UNDEFID;
  if ( setVdimname )
    {
      char vdimname[CDI_MAX_NAME+3]; vdimname[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_VDIMNAME, CDI_MAX_NAME, vdimname);
      if ( vdimname[0] == 0 ) strcpy(vdimname, vdimname_default);
      nvdimID = dimIDs[ndims-1] = checkDimName(fileID, nvertex, vdimname);
      if ( nvdimID == CDI_UNDEFID )
        {
          cdf_def_dim(fileID, vdimname, nvertex, dimIDs+ndims-1);
          nvdimID = dimIDs[ndims-1];
        }
    }

  if ( ndims == 3 )
    {
      dimIDs[0] = ydimID;
      dimIDs[1] = xdimID;
    }
  else /* ndims == 2 */
    {
      dimIDs[0] = xdimID;
      cdfDefGridReference(streamptr, gridID);
      cdfDefGridUUID(streamptr, gridID);
    }

  const double *xvalsPtr = gridInqXvalsPtr(gridID),
    *xboundsPtr = NULL;
  if ( xvalsPtr )
    {
      char xaxisname[CDI_MAX_NAME]; xaxisname[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, xaxisname);
      checkGridName(xaxisname, fileID);
      cdf_def_var(fileID, xaxisname, xtype, ndims-1, dimIDs, &ncxvarid);
      cdfGridCompress(fileID, ncxvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

      cdfPutGridStdAtts(fileID, ncxvarid, gridID, 'X', &gridInqsX);

      /* attribute for Panoply */
      if ( ndims == 3 )
        cdf_put_att_text(fileID, ncxvarid, "_CoordinateAxisType", 3, "Lon");

      if ( (xboundsPtr = gridInqXboundsPtr(gridID)) && nvdimID != CDI_UNDEFID )
        {
          size_t xaxisnameLen = strlen(xaxisname);
          xaxisname[xaxisnameLen] = '_';
          memcpy(xaxisname + xaxisnameLen + 1, bndsName, sizeof (bndsName));
          cdf_def_var(fileID, xaxisname, xtype, ndims, dimIDs, &ncbxvarid);
          cdfGridCompress(fileID, ncbxvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

          cdf_put_att_text(fileID, ncxvarid, "bounds", xaxisnameLen + sizeof (bndsName), xaxisname);
        }
    }

  const double *yvalsPtr = gridInqYvalsPtr(gridID),
    *yboundsPtr = NULL;
  if ( yvalsPtr )
    {
      char yaxisname[CDI_MAX_NAME];
      gridInqYname(gridID, yaxisname);
      checkGridName(yaxisname, fileID);

      cdf_def_var(fileID, yaxisname, xtype, ndims - 1, dimIDs, &ncyvarid);
      cdfGridCompress(fileID, ncyvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

      cdfPutGridStdAtts(fileID, ncyvarid, gridID, 'Y', &gridInqsY);

      /* attribute for Panoply */
      if ( ndims == 3 )
        cdf_put_att_text(fileID, ncyvarid, "_CoordinateAxisType", 3, "Lat");

      if ( (yboundsPtr = gridInqYboundsPtr(gridID)) && nvdimID != CDI_UNDEFID )
        {
          size_t yaxisnameLen = strlen(yaxisname);
          yaxisname[yaxisnameLen] = '_';
          memcpy(yaxisname + yaxisnameLen + 1, bndsName, sizeof (bndsName));
          cdf_def_var(fileID, yaxisname, xtype, ndims, dimIDs, &ncbyvarid);
          cdfGridCompress(fileID, ncbyvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

          cdf_put_att_text(fileID, ncyvarid, "bounds", yaxisnameLen + sizeof (bndsName), yaxisname);
        }
    }

  const double *areaPtr = gridInqAreaPtr(gridID);
  if ( areaPtr )
    {
      static const char yaxisname_[] = "cell_area";
      static const char units[] = "m2";
      static const char longname[] = "area of grid cell";
      static const char stdname[] = "cell_area";

      cdf_def_var(fileID, yaxisname_, xtype, ndims-1, dimIDs, &ncavarid);

      cdf_put_att_text(fileID, ncavarid, "standard_name", sizeof (stdname) - 1, stdname);
      cdf_put_att_text(fileID, ncavarid, "long_name", sizeof (longname) - 1, longname);
      cdf_put_att_text(fileID, ncavarid, "units", sizeof (units) - 1, units);
    }

  cdf_enddef(fileID);
  streamptr->ncmode = 2;

  if ( ncxvarid  != CDI_UNDEFID ) cdf_put_var_double(fileID, ncxvarid,  xvalsPtr);
  if ( ncbxvarid != CDI_UNDEFID ) cdf_put_var_double(fileID, ncbxvarid, xboundsPtr);
  if ( ncyvarid  != CDI_UNDEFID ) cdf_put_var_double(fileID, ncyvarid,  yvalsPtr);
  if ( ncbyvarid != CDI_UNDEFID ) cdf_put_var_double(fileID, ncbyvarid, yboundsPtr);
  if ( ncavarid  != CDI_UNDEFID ) cdf_put_var_double(fileID, ncavarid,  areaPtr);

  return (struct cdfDefIrregularGridCommonIDs) {
    .xdimID=xdimID, .ydimID = ydimID,
    .ncxvarid=ncxvarid, .ncyvarid=ncyvarid, .ncavarid=ncavarid
  };
}

static
void cdfDefCurvilinear(stream_t *streamptr, int gridID, int gridindex)
{
  ncgrid_t *ncgrid = streamptr->ncgrid;

  size_t dimlen = (size_t)gridInqSize(gridID);
  size_t xdimlen = (size_t)gridInqXsize(gridID);
  size_t ydimlen = (size_t)gridInqYsize(gridID);

  int xdimID = CDI_UNDEFID, ydimID = CDI_UNDEFID;
  int ncxvarid = CDI_UNDEFID, ncyvarid = CDI_UNDEFID, ncavarid = CDI_UNDEFID;
  {
    size_t ofs = 0;
    do {
      struct idSearch search
        = cdfSearchIDBySize(ofs, (size_t)gridindex, ncgrid, CDF_DIMID_X,
                            GRID_CURVILINEAR, (int)dimlen,
                            gridInqType, gridInqSize);
      size_t index = search.foundIdx;
      if ( index != SIZE_MAX )
        {
          int gridID0 = ncgrid[index].gridID;
          if (    IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0))
               && IS_EQUAL(gridInqXval(gridID0, (int)dimlen-1),
                           gridInqXval(gridID, (int)dimlen-1))
               && IS_EQUAL(gridInqYval(gridID0, 0), gridInqYval(gridID, 0))
               && IS_EQUAL(gridInqYval(gridID0, (int)dimlen-1),
                           gridInqYval(gridID, (int)dimlen-1)) )
            {
              xdimID = ncgrid[index].ncIDs[CDF_DIMID_X];
              ydimID = ncgrid[index].ncIDs[CDF_DIMID_Y];
              ncxvarid = ncgrid[index].ncIDs[CDF_VARID_X];
              ncyvarid = ncgrid[index].ncIDs[CDF_VARID_Y];
              break;
            }
          ofs = search.foundIdx;
          if ( ofs < (size_t)gridindex )
            continue;
        }
    } while (false);
  }

  if ( xdimID == CDI_UNDEFID || ydimID == CDI_UNDEFID )
    {
      struct cdfDefIrregularGridCommonIDs createdIDs
        = cdfDefIrregularGridCommon(streamptr, gridID,
                                    xdimlen, ydimlen, 3, "x", 4, "nv4",
                                    gridInqXboundsPtr(gridID)
                                    || gridInqYboundsPtr(gridID));
      xdimID = createdIDs.xdimID;
      ydimID = createdIDs.ydimID;
      ncxvarid = createdIDs.ncxvarid;
      ncyvarid = createdIDs.ncyvarid;
      ncavarid = createdIDs.ncavarid;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_X] = xdimID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_Y] = ydimID;
  ncgrid[gridindex].ncIDs[CDF_VARID_X] = ncxvarid;
  ncgrid[gridindex].ncIDs[CDF_VARID_Y] = ncyvarid;
  ncgrid[gridindex].ncIDs[CDF_VARID_A] = ncavarid;
}


static
void cdfDefUnstructured(stream_t *streamptr, int gridID, int gridindex)
{
  ncgrid_t *ncgrid = streamptr->ncgrid;

  size_t dimlen = (size_t)gridInqSize(gridID);

  int dimID = CDI_UNDEFID;
  int ncxvarid = CDI_UNDEFID, ncyvarid = CDI_UNDEFID, ncavarid = CDI_UNDEFID;
  {
    size_t ofs = 0;
    do {
      struct idSearch search
        = cdfSearchIDBySize(ofs, (size_t)gridindex, ncgrid, CDF_DIMID_X,
                            GRID_UNSTRUCTURED, (int)dimlen,
                            gridInqType, gridInqSize);
      size_t index = search.foundIdx;
      if ( index != SIZE_MAX )
        {
          int gridID0 = ncgrid[index].gridID;
          if ( gridInqNvertex(gridID0) == gridInqNvertex(gridID) &&
               IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
               IS_EQUAL(gridInqXval(gridID0, (int)dimlen-1),
                        gridInqXval(gridID, (int)dimlen-1)) &&
               IS_EQUAL(gridInqYval(gridID0, 0), gridInqYval(gridID, 0)) &&
               IS_EQUAL(gridInqYval(gridID0, (int)dimlen-1),
                        gridInqYval(gridID, (int)dimlen-1)) )
            {
              dimID = ncgrid[index].ncIDs[CDF_DIMID_X];
              ncxvarid = ncgrid[index].ncIDs[CDF_VARID_X];
              ncyvarid = ncgrid[index].ncIDs[CDF_VARID_Y];
              ncavarid = ncgrid[index].ncIDs[CDF_VARID_A];
              break;
            }
          ofs = search.foundIdx;
          if ( ofs < (size_t)gridindex )
            continue;
        }
    } while (false);
  }

  if ( dimID == CDI_UNDEFID )
    {
      size_t nvertex = (size_t)gridInqNvertex(gridID);
      struct cdfDefIrregularGridCommonIDs createdIDs
        = cdfDefIrregularGridCommon(streamptr, gridID,
                                    dimlen, 1, 2, "ncells",
                                    nvertex, "vertices", nvertex > 0);
      dimID = createdIDs.xdimID;
      ncxvarid = createdIDs.ncxvarid;
      ncyvarid = createdIDs.ncyvarid;
      ncavarid = createdIDs.ncavarid;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_X] = dimID;
  ncgrid[gridindex].ncIDs[CDF_VARID_X] = ncxvarid;
  ncgrid[gridindex].ncIDs[CDF_VARID_Y] = ncyvarid;
  ncgrid[gridindex].ncIDs[CDF_VARID_A] = ncavarid;
}

struct attTxtTab2
{
  const char *attName, *attVal;
  size_t valLen;
};

static
void cdf_def_vct_echam(stream_t *streamptr, int zaxisID)
{
  int type = zaxisInqType(zaxisID);

  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int ilev = zaxisInqVctSize(zaxisID)/2;
      if ( ilev == 0 ) return;

      int mlev = ilev - 1;

      if ( streamptr->vct.ilev > 0 )
        {
          if ( streamptr->vct.ilev != ilev )
            Error("More than one VCT for each file unsupported!");
          return;
        }

      int fileID = streamptr->fileID;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      int ncdimid = -1, ncdimid2 = -1;
      int hyaiid, hybiid, hyamid = -1, hybmid = -1;

      cdf_def_dim(fileID, "nhyi", (size_t)ilev, &ncdimid2);
      cdf_def_var(fileID, "hyai", NC_DOUBLE, 1, &ncdimid2, &hyaiid);
      cdf_def_var(fileID, "hybi", NC_DOUBLE, 1, &ncdimid2, &hybiid);
      if ( mlev > 0 )
        {
          cdf_def_dim(fileID, "nhym", (size_t)mlev, &ncdimid);
          cdf_def_var(fileID, "hyam", NC_DOUBLE, 1, &ncdimid,  &hyamid);
          cdf_def_var(fileID, "hybm", NC_DOUBLE, 1, &ncdimid,  &hybmid);
        }

      streamptr->vct.ilev   = ilev;
      streamptr->vct.mlev   = mlev;
      streamptr->vct.mlevID = ncdimid;
      streamptr->vct.ilevID = ncdimid2;

      {
        static const char lname_n[] = "long_name",
          units_n[] = "units",
          lname_v_ai[] = "hybrid A coefficient at layer interfaces",
          units_v_ai[] = "Pa",
          lname_v_bi[] = "hybrid B coefficient at layer interfaces",
          units_v_bi[] = "1";
        static const struct attTxtTab2 tab[]
          = {
          { lname_n, lname_v_ai, sizeof (lname_v_ai) - 1 },
          { units_n, units_v_ai, sizeof (units_v_ai) - 1 },
          { lname_n, lname_v_bi, sizeof (lname_v_bi) - 1 },
          { units_n, units_v_bi, sizeof (units_v_bi) - 1 },
        };
        enum { tabLen = sizeof (tab) / sizeof (tab[0]) };
        int ids[tabLen] = { hyaiid, hyaiid, hybiid, hybiid };
        for ( size_t i = 0; i < tabLen; ++i )
          cdf_put_att_text(fileID, ids[i], tab[i].attName, tab[i].valLen, tab[i].attVal);
      }

      {
        static const char lname_n[] = "long_name",
          units_n[] = "units",
          lname_v_am[] = "hybrid A coefficient at layer midpoints",
          units_v_am[] = "Pa",
          lname_v_bm[] = "hybrid B coefficient at layer midpoints",
          units_v_bm[] = "1";
        static const struct attTxtTab2 tab[]
          = {
          { lname_n, lname_v_am, sizeof (lname_v_am) - 1 },
          { units_n, units_v_am, sizeof (units_v_am) - 1 },
          { lname_n, lname_v_bm, sizeof (lname_v_bm) - 1 },
          { units_n, units_v_bm, sizeof (units_v_bm) - 1 },
        };
        enum { tabLen = sizeof (tab) / sizeof (tab[0]) };
        int ids[tabLen] = { hyamid, hyamid, hybmid, hybmid };
        for ( size_t i = 0; i < tabLen; ++i )
          cdf_put_att_text(fileID, ids[i], tab[i].attName, tab[i].valLen, tab[i].attVal);
      }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      const double *vctptr = zaxisInqVctPtr(zaxisID);

      cdf_put_var_double(fileID, hyaiid, vctptr);
      cdf_put_var_double(fileID, hybiid, vctptr+ilev);

      size_t start;
      size_t count = 1;
      double mval;
      for ( int i = 0; i < mlev; i++ )
        {
          start = (size_t)i;
          mval = (vctptr[i] + vctptr[i+1]) * 0.5;
          cdf_put_vara_double(fileID, hyamid, &start, &count, &mval);
          mval = (vctptr[ilev+i] + vctptr[ilev+i+1]) * 0.5;
          cdf_put_vara_double(fileID, hybmid, &start, &count, &mval);
        }
    }
}

static
void cdf_def_vct_cf(stream_t *streamptr, int zaxisID, int nclevID, int ncbndsID, int p0status, double p0value)
{
  int type = zaxisInqType(zaxisID);

  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int ilev = zaxisInqVctSize(zaxisID)/2;
      if ( ilev == 0 ) return;

      int mlev = ilev - 1;
      int hyaiid = 0, hybiid = 0, hyamid, hybmid;

      if ( streamptr->vct.ilev > 0 )
        {
          if ( streamptr->vct.ilev != ilev )
            Error("more than one VCT for each file unsupported!");
          return;
        }

      int fileID = streamptr->fileID;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      int dimIDs[2];
      dimIDs[0] = nclevID;
      dimIDs[1] = ncbndsID;

      streamptr->vct.mlev   = mlev;
      streamptr->vct.ilev   = ilev;
      streamptr->vct.mlevID = nclevID;
      streamptr->vct.ilevID = nclevID;

      if ( p0status == 0 )
        cdf_def_var(fileID, "a", NC_DOUBLE, 1, dimIDs,  &hyamid);
      else
        cdf_def_var(fileID, "ap", NC_DOUBLE, 1, dimIDs,  &hyamid);
      cdf_def_var(fileID, "b",  NC_DOUBLE, 1, dimIDs,  &hybmid);

      {
        static const char lname[] = "vertical coordinate formula term: ap(k)";
        cdf_put_att_text(fileID, hyamid, "long_name", sizeof (lname) - 1, lname);
      }
      {
        static const char units[] = "Pa";
        cdf_put_att_text(fileID, hyamid, "units", sizeof (units) - 1, units);
      }
      {
        static const char lname[] = "vertical coordinate formula term: b(k)";
        cdf_put_att_text(fileID, hybmid, "long_name", sizeof (lname) - 1, lname);
      }
      {
        static const char units[] = "1";
        cdf_put_att_text(fileID, hybmid, "units", sizeof (units) - 1, units);
      }

      if ( ncbndsID != -1 )
        {
          if ( p0status == 0 )
            cdf_def_var(fileID, "a_bnds", NC_DOUBLE, 2, dimIDs, &hyaiid);
          else
            cdf_def_var(fileID, "ap_bnds", NC_DOUBLE, 2, dimIDs, &hyaiid);
          cdf_def_var(fileID, "b_bnds",  NC_DOUBLE, 2, dimIDs, &hybiid);
          {
            static const char lname[] = "vertical coordinate formula term: ap(k+1/2)";
            cdf_put_att_text(fileID, hyaiid, "long_name", sizeof (lname) - 1, lname);
          }
          {
            static const char units[] = "Pa";
            cdf_put_att_text(fileID, hyaiid, "units", sizeof (units) - 1, units);
          }
          {
            static const char lname[] = "vertical coordinate formula term: b(k+1/2)";
            cdf_put_att_text(fileID, hybiid, "long_name", sizeof (lname) - 1, lname);
          }
          {
            static const char units[] = "1";
            cdf_put_att_text(fileID, hybiid, "units", sizeof (units) - 1, units);
          }
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      int vctsize = zaxisInqVctSize(zaxisID);
      double vct[vctsize];
      zaxisInqVct(zaxisID, vct);

      if ( p0status == 0 && IS_NOT_EQUAL(p0value,0) )
        for ( int i = 0; i < vctsize/2; ++i ) vct[i] /= p0value;

      double tarray[ilev*2];

      if ( ncbndsID != -1 )
        {
          for ( int i = 0; i < mlev; ++i )
            {
              tarray[2*i  ] = vct[i];
              tarray[2*i+1] = vct[i+1];
            }
          cdf_put_var_double(fileID, hyaiid, tarray);

          for ( int i = 0; i < mlev; ++i )
            {
              tarray[2*i  ] = vct[ilev+i];
              tarray[2*i+1] = vct[ilev+i+1];
            }
          cdf_put_var_double(fileID, hybiid, tarray);
        }

      for ( int i = 0; i < mlev; ++i )
        tarray[i] = (vct[i] + vct[i+1]) * 0.5;
      cdf_put_var_double(fileID, hyamid, tarray);

      for ( int i = 0; i < mlev; ++i )
        tarray[i] = (vct[ilev+i] + vct[ilev+i+1]) * 0.5;
      cdf_put_var_double(fileID, hybmid, tarray);
    }
}

struct attTxtTab { const char *txt; size_t txtLen; };

static
void cdf_def_zaxis_hybrid_echam(stream_t *streamptr, int type, int *ncvaridp, int zaxisID, int zaxisindex, int xtype, size_t dimlen, int *dimID, char *axisname)
{
  int fileID = streamptr->fileID;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  cdf_def_dim(fileID, axisname, dimlen, dimID);
  cdf_def_var(fileID, axisname, (nc_type) xtype, 1, dimID,  ncvaridp);
  int ncvarid = *ncvaridp;

  {
    static const char sname[] = "hybrid_sigma_pressure";
    cdf_put_att_text(fileID, ncvarid, "standard_name", sizeof (sname) - 1, sname);
  }
  {
    static const char *attName[] = {
      "long_name",
      "formula",
      "formula_terms"
    };
    enum { nAtt = sizeof (attName) / sizeof (attName[0]) };
    static const char lname_m[] = "hybrid level at layer midpoints",
      formula_m[] = "hyam hybm (mlev=hyam+hybm*aps)",
      fterms_m[] = "ap: hyam b: hybm ps: aps",
      lname_i[] = "hybrid level at layer interfaces",
      formula_i[] = "hyai hybi (ilev=hyai+hybi*aps)",
      fterms_i[] = "ap: hyai b: hybi ps: aps";
    static const struct attTxtTab tab[2][nAtt] = {
      {
        { lname_i, sizeof (lname_i) - 1 },
        { formula_i, sizeof (formula_i) - 1 },
        { fterms_i, sizeof (fterms_i) - 1 }
      },
      {
        { lname_m, sizeof (lname_m) - 1 },
        { formula_m, sizeof (formula_m) - 1 },
        { fterms_m, sizeof (fterms_m) - 1 }
      }
    };

    size_t tabSelect = type == ZAXIS_HYBRID;
    for (size_t i = 0; i < nAtt; ++i)
      cdf_put_att_text(fileID, ncvarid, attName[i],
                       tab[tabSelect][i].txtLen, tab[tabSelect][i].txt);
  }

  {
    static const char units[] = "level";
    cdf_put_att_text(fileID, ncvarid, "units", sizeof (units) - 1, units);
  }
  {
    static const char direction[] = "down";
    cdf_put_att_text(fileID, ncvarid, "positive", sizeof (direction) - 1, direction);
  }

  cdf_enddef(fileID);
  streamptr->ncmode = 2;

  cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));

  cdf_def_vct_echam(streamptr, zaxisID);

  if ( *dimID == CDI_UNDEFID )
    streamptr->zaxisID[zaxisindex] = type == ZAXIS_HYBRID
      ? streamptr->vct.mlevID : streamptr->vct.ilevID;
}

static
void cdf_def_zaxis_hybrid_cf(stream_t *streamptr, int type, int *ncvaridp, int zaxisID, int zaxisindex, int xtype, size_t dimlen, int *dimID, char *axisname)
{
  int fileID = streamptr->fileID;
  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  char psname[CDI_MAX_NAME]; psname[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_PSNAME, CDI_MAX_NAME, psname);
  if ( psname[0] == 0 ) strcpy(psname, "ps");

  char p0name[CDI_MAX_NAME]; p0name[0] = 0;
  double p0value = 1;
  int p0varid = CDI_UNDEFID;
  int p0status = cdiZaxisInqKeyFlt(zaxisID, CDI_KEY_P0VALUE, &p0value);
  if ( p0status == 0 )
    {
      cdiZaxisInqKeyStr(zaxisID, CDI_KEY_P0NAME, CDI_MAX_NAME, p0name);
      if ( p0name[0] == 0 ) strcpy(p0name, "p0");
      cdf_def_var(fileID, p0name, NC_DOUBLE, 0, 0,  &p0varid);
      static const char longname[] = "reference pressure";
      cdf_put_att_text(fileID, p0varid, "long_name", strlen(longname), longname);
      static const char units[] = "Pa";
      cdf_put_att_text(fileID, p0varid, "units", strlen(units), units);
    }

  char zname[CDI_MAX_NAME]; zname[0] = 0;
  char zlongname[CDI_MAX_NAME]; zlongname[0] = 0;
  char zunits[CDI_MAX_NAME]; zunits[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_NAME, CDI_MAX_NAME, zname);
  //cdiZaxisInqKeyStr(zaxisID, CDI_KEY_LONGNAME, CDI_MAX_NAME, zlongname);
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_UNITS, CDI_MAX_NAME, zunits);
  if ( zname[0] ) strcpy(axisname, zname);
  if ( zlongname[0] == 0 ) strcpy(zlongname, "hybrid sigma pressure coordinate");
  if ( zunits[0] == 0 ) strcpy(zunits, "1");

  cdf_def_dim(fileID, axisname, dimlen, dimID);
  cdf_def_var(fileID, axisname, (nc_type) xtype, 1, dimID, ncvaridp);
  int ncvarid = *ncvaridp;

  {
    static const char sname[] = "standard_name",
      sname_v[] = "atmosphere_hybrid_sigma_pressure_coordinate",
      axis[] = "axis",
      axis_v[] = "Z",
      direction[] = "positive",
      direction_v[] = "down";
    struct attTxtTab2 tab[] = {
      { sname, sname_v, sizeof (sname_v) - 1 },
      { axis, axis_v, sizeof (axis_v) - 1 },
      { direction, direction_v, sizeof (direction_v) - 1 },
    };
    enum { nAtt = sizeof (tab) / sizeof (tab[0]) };
    for ( size_t i = 0; i < nAtt; ++i )
      cdf_put_att_text(fileID, ncvarid, tab[i].attName, tab[i].valLen, tab[i].attVal);

    cdf_put_att_text(fileID, ncvarid, "long_name", strlen(zlongname), zlongname);
    cdf_put_att_text(fileID, ncvarid, "units", strlen(zunits), zunits);
  }

  size_t len = 0;
  char txt[CDI_MAX_NAME];
  if ( p0status == 0 )
    len = (size_t)(sprintf(txt, "%s%s %s%s", "a: a b: b p0: ", p0name, "ps: ", psname));
  else
    len = (size_t)(sprintf(txt, "%s%s", "ap: ap b: b ps: ", psname));
  cdf_put_att_text(fileID, ncvarid, "formula_terms", len, txt);

  int ncbvarid = CDI_UNDEFID;
  int nvdimID = CDI_UNDEFID;

  double lbounds[dimlen], ubounds[dimlen], levels[dimlen];

  if ( zaxisInqLevels(zaxisID, NULL) )
    zaxisInqLevels(zaxisID, levels);
  else
    for ( size_t i = 0; i < dimlen; ++i ) levels[i] = i+1;

  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      zaxisInqLbounds(zaxisID, lbounds);
      zaxisInqUbounds(zaxisID, ubounds);
    }
  else
    {
      for ( size_t i = 0; i < dimlen; ++i ) lbounds[i] = levels[i];
      for ( size_t i = 0; i < dimlen-1; ++i ) ubounds[i] = levels[i+1];
      ubounds[dimlen-1] = levels[dimlen-1] + 1;
    }

  //if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      size_t nvertex = 2;
      if ( dimlen > 1 && nc_inq_dimid(fileID, bndsName, &nvdimID) != NC_NOERR )
        cdf_def_dim(fileID, bndsName, nvertex, &nvdimID);

      if ( nvdimID != CDI_UNDEFID )
        {
          size_t axisnameLen = strlen(axisname);
          axisname[axisnameLen] = '_';
          memcpy(axisname + axisnameLen + 1, bndsName, sizeof (bndsName));
          axisnameLen += sizeof (bndsName);
          int dimIDs[2] = { *dimID, nvdimID };
          cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
          cdf_put_att_text(fileID, ncvarid, "bounds", axisnameLen, axisname);
          {
            static const char sname[] = "standard_name",
              sname_v[] = "atmosphere_hybrid_sigma_pressure_coordinate";
            struct attTxtTab2 tab[] = {
              { sname, sname_v, sizeof (sname_v) - 1 },
            };
            enum { nAtt = sizeof (tab) / sizeof (tab[0]) };
            for ( size_t i = 0; i < nAtt; ++i )
              cdf_put_att_text(fileID, ncbvarid, tab[i].attName, tab[i].valLen, tab[i].attVal);
            cdf_put_att_text(fileID, ncbvarid, "units", strlen(zunits), zunits);
          }

          if ( p0status == 0 )
            len = (size_t)(sprintf(txt, "%s%s %s%s", "a: a_bnds b: b_bnds p0: ", p0name, "ps: ", psname));
          else
            len = (size_t)(sprintf(txt, "%s%s", "ap: ap_bnds b: b_bnds ps: ", psname));
          cdf_put_att_text(fileID, ncbvarid, "formula_terms", len, txt);
        }
    }

  cdf_enddef(fileID);
  streamptr->ncmode = 2;

  cdf_put_var_double(fileID, ncvarid, levels);

  if ( p0varid != CDI_UNDEFID ) cdf_put_var_double(fileID, p0varid, &p0value);

  if ( ncbvarid != CDI_UNDEFID )
    {
      double zbounds[2*dimlen];
      for ( size_t i = 0; i < dimlen; ++i )
        {
          zbounds[2*i  ] = lbounds[i];
          zbounds[2*i+1] = ubounds[i];
        }
      cdf_put_var_double(fileID, ncbvarid, zbounds);
    }

  cdf_def_vct_cf(streamptr, zaxisID, *dimID, nvdimID, p0status, p0value);

  if ( *dimID == CDI_UNDEFID )
    streamptr->zaxisID[zaxisindex] = type == ZAXIS_HYBRID
      ? streamptr->vct.mlevID : streamptr->vct.ilevID;
}

static
void cdf_def_zaxis_hybrid(stream_t *streamptr, int type, int *ncvarid, int zaxisID, int zaxisindex, int xtype, size_t dimlen, int *dimID, char *axisname)
{
  void (*def_zaxis_hybrid_delegate)(stream_t *streamptr, int type, int *ncvarid, int zaxisID, int zaxisindex, int xtype, size_t dimlen, int *dimID, char *axisname)
    = ( (!CDI_cmor_mode && cdiConvention == CDI_CONVENTION_ECHAM)
        || type == ZAXIS_HYBRID_HALF )
    ? cdf_def_zaxis_hybrid_echam : cdf_def_zaxis_hybrid_cf;
  def_zaxis_hybrid_delegate(streamptr, type, ncvarid, zaxisID, zaxisindex, xtype, dimlen, dimID, axisname);
}

static
void cdfDefZaxisUUID(stream_t *streamptr, int zaxisID)
{
  unsigned char uuidOfVGrid[CDI_UUID_SIZE];
  zaxisInqUUID(zaxisID, uuidOfVGrid);

  if ( uuidOfVGrid[0] != 0 )
    {
      char uuidOfVGridStr[37];
      cdiUUID2Str(uuidOfVGrid, uuidOfVGridStr);
      if ( uuidOfVGridStr[0] != 0 && strlen(uuidOfVGridStr) == 36 )
        {
          int fileID  = streamptr->fileID;
          if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
          cdf_put_att_text(fileID, NC_GLOBAL, "uuidOfVGrid", 36, uuidOfVGridStr);
          if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
        }
    }
}

static
void cdfDefZaxisChar(stream_t *streamptr, int zaxisID, char *axisname, int *dimID, size_t dimlen, int zaxisindex)
{
  int fileID  = streamptr->fileID;
  int ncvarID = CDI_UNDEFID;
  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  /* Check StrlenID */
  char strlen[7] = "strlen\0";
  size_t clen = (size_t) zaxisInqCLen(zaxisID);
  if ( clen == 0 )
    Error("Maximal string length value is 0.\nA given character axis requires a dimension to save the maximal string length.");
  int strlenID = CDI_UNDEFID;
  strlenID = checkDimName(fileID, clen, strlen);

  if ( strlenID == CDI_UNDEFID ) cdf_def_dim(fileID, strlen, clen, &strlenID);

  /* Check 'areatype'dimID */
  char dimname[CDI_MAX_NAME+3]; dimname[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_DIMNAME, CDI_MAX_NAME, dimname);
  *dimID = checkDimName(fileID, dimlen, dimname);
  if ( !(dimlen > 0) )
    Error("No strings delivered for a character axis.");
  if ( dimname[0] == 0 ) { memcpy(dimname, "area_type", 10); dimname[10] = 0; }

  if ( *dimID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, dimlen, dimID);

  int dimIDs[2];
  dimIDs[0] = *dimID;
  dimIDs[1] = strlenID;

  /* Get Stringvalues */
  char **cvals = zaxisInqCValsPtr(zaxisID);

  if ( cvals )
    {
      /* Define variable and its attributes */
      cdf_def_var(fileID, axisname, NC_CHAR, 2, dimIDs, &ncvarID);

      cdfPutGridStdAtts(fileID, ncvarID, zaxisID, 'Z', &gridInqsZ);
      cdf_put_att_text(fileID, ncvarID, "axis", 1, "Z");
      cdfDefineAttributes(zaxisID, CDI_GLOBAL, fileID, ncvarID);

      streamptr->nczvarID[zaxisindex] = ncvarID;
      cdf_enddef(fileID);

      /* Write Stringvalues */
      size_t start[2], count[2];
      start[1] = 0;
      count[0] = 1;
      count[1] = clen;
      for ( size_t i = 0; i < dimlen; i++ )
        {
          start[0] = i;
          nc_put_vara_text(fileID, ncvarID, start, count, cvals[i]);
        }
    }

  streamptr->ncmode = 2;
}

static
void cdfDefZaxis(stream_t *streamptr, int zaxisID)
{
  /*  char zaxisname0[CDI_MAX_NAME]; */
  char axisname[CDI_MAX_NAME];
  int dimID = CDI_UNDEFID;
  int ncvarid = CDI_UNDEFID, ncbvarid = CDI_UNDEFID;
  int xtype = zaxisInqDatatype(zaxisID) == CDI_DATATYPE_FLT32 ? NC_FLOAT : NC_DOUBLE;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);

  int nzaxis = vlistNzaxis(vlistID);

  size_t dimlen = (size_t)zaxisInqSize(zaxisID);
  int type = zaxisInqType(zaxisID);

  int ndims = 1;

  if ( dimlen == 1 )
    {
      bool is_scalar = zaxisInqScalar(zaxisID) > 0;
      if ( !is_scalar && CDI_cmor_mode )
        {
          is_scalar = true;
          zaxisDefScalar(zaxisID);
        }

      if ( is_scalar ) ndims = 0;
      if ( CDI_reduce_dim ) return;

      switch (type)
        {
        case ZAXIS_SURFACE:
        case ZAXIS_CLOUD_BASE:
        case ZAXIS_CLOUD_TOP:
        case ZAXIS_ISOTHERM_ZERO:
        case ZAXIS_TOA:
        case ZAXIS_SEA_BOTTOM:
        case ZAXIS_ATMOSPHERE:
        case ZAXIS_MEANSEA:
        case ZAXIS_LAKE_BOTTOM:
        case ZAXIS_SEDIMENT_BOTTOM:
        case ZAXIS_SEDIMENT_BOTTOM_TA:
        case ZAXIS_SEDIMENT_BOTTOM_TW:
        case ZAXIS_MIX_LAYER:
          return;
        }
    }

  zaxisInqName(zaxisID, axisname);

  if ( dimID == CDI_UNDEFID )
    {
      checkZaxisName(axisname, fileID, vlistID, zaxisID, nzaxis);

      char dimname[CDI_MAX_NAME+3]; dimname[0] = 0;
      //cdiZaxisInqKeyStr(zaxisID, CDI_KEY_DIMNAME, CDI_MAX_NAME, dimname);
      if ( dimname[0] == 0 ) strcpy(dimname, axisname);

      if ( type == ZAXIS_REFERENCE ) cdfDefZaxisUUID(streamptr, zaxisID);

      if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
        {
          cdf_def_zaxis_hybrid(streamptr, type, &ncvarid, zaxisID, zaxisindex, xtype, dimlen, &dimID, axisname);

          int natts;
          cdiInqNatts(zaxisID, CDI_GLOBAL, &natts);
          if ( natts > 0 && streamptr->ncmode == 2 ) cdf_redef(fileID);
          cdfDefineAttributes(zaxisID, CDI_GLOBAL, fileID, ncvarid);
          if ( natts > 0 && streamptr->ncmode == 2 ) cdf_enddef(fileID);
        }
      else if ( type == ZAXIS_CHAR )
        cdfDefZaxisChar(streamptr, zaxisID, axisname, &dimID, dimlen, zaxisindex);
      else
        {
          dimID = checkDimName(fileID, dimlen, dimname);

          if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

          if ( ndims && dimID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, dimlen, &dimID);

          if ( zaxisInqLevels(zaxisID, NULL) )
            {
              cdf_def_var(fileID, axisname, (nc_type) xtype, ndims, &dimID, &ncvarid);

              cdfPutGridStdAtts(fileID, ncvarid, zaxisID, 'Z', &gridInqsZ);

              {
                int positive = zaxisInqPositive(zaxisID);
                static const char positive_up[] = "up",
                                  positive_down[] = "down";
                static const struct attTxtTab tab[2] = {
                  { positive_up, sizeof (positive_up) - 1 },
                  { positive_down, sizeof (positive_down) - 1 },
                };
                if ( positive == POSITIVE_UP || positive == POSITIVE_DOWN )
                  {
                    size_t select = positive == POSITIVE_DOWN;
                    cdf_put_att_text(fileID, ncvarid, "positive", tab[select].txtLen, tab[select].txt);
                  }
              }
              cdf_put_att_text(fileID, ncvarid, "axis", 1, "Z");

              if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
                {
                  int nvdimID = CDI_UNDEFID;
                  size_t nvertex = 2;
                  if ( nc_inq_dimid(fileID, bndsName, &nvdimID) != NC_NOERR )
                    cdf_def_dim(fileID, bndsName, nvertex, &nvdimID);

                  if ( nvdimID != CDI_UNDEFID )
                    {
                      size_t axisnameLen = strlen(axisname);
                      axisname[axisnameLen] = '_';
                      memcpy(axisname + axisnameLen + 1, bndsName, sizeof (bndsName));
                      int dimIDs[2];
                      dimIDs[0] = dimID;
                      dimIDs[ndims] = nvdimID;
                      cdf_def_var(fileID, axisname, (nc_type) xtype, ndims+1, dimIDs, &ncbvarid);
                      cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
                    }
                }

              cdfDefineAttributes(zaxisID, CDI_GLOBAL, fileID, ncvarid);
            }

          cdf_enddef(fileID);
          streamptr->ncmode = 2;

          if ( zaxisInqLevels(zaxisID, NULL) )
            {
              cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));

              if ( ncbvarid != CDI_UNDEFID )
                {
                  double lbounds[dimlen], ubounds[dimlen], zbounds[2*dimlen];
                  zaxisInqLbounds(zaxisID, lbounds);
                  zaxisInqUbounds(zaxisID, ubounds);
                  for ( size_t i = 0; i < dimlen; ++i )
                    {
                      zbounds[2*i  ] = lbounds[i];
                      zbounds[2*i+1] = ubounds[i];
                    }

                  cdf_put_var_double(fileID, ncbvarid, zbounds);
                }

              if ( ndims == 0 ) streamptr->nczvarID[zaxisindex] = ncvarid;
            }
        }
    }

  if ( dimID != CDI_UNDEFID )
    streamptr->zaxisID[zaxisindex] = dimID;
}

static
void cdf_def_mapping(stream_t *streamptr, int gridID)
{
  char mapping[CDI_MAX_NAME]; mapping[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
  if ( mapping[0] )
    {
      char gmapvarname[CDI_MAX_NAME]; gmapvarname[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPPING, CDI_MAX_NAME, gmapvarname);

      int fileID = streamptr->fileID;
      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      int ncvarid;
      int ncerrcode = nc_def_var(fileID, gmapvarname, (nc_type) NC_INT, 0, NULL, &ncvarid);
      if ( ncerrcode == NC_NOERR )
        cdfDefineAttributes(gridID, CDI_GLOBAL, fileID, ncvarid);

      cdf_enddef(fileID);

      if ( ncerrcode == NC_NOERR )
        {
          int dummy = 1;
          cdf_put_var_int(fileID, ncvarid, &dummy);
        }
    }
}

static
void cdfDefCharacter(stream_t *streamptr, int gridID, int gridindex, int xory, int strlen)
{
  if ( streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X] != CDI_UNDEFID ) return;

  int dimlen = ( xory == 0 ) ? gridInqXsize(gridID) : gridInqYsize(gridID);
  int dimID, strlenID;
  ncgrid_t *ncgrid = streamptr->ncgrid;

  /* Check for all grids up to gridindex whether it already is defined */

  for ( int index = 0; index < gridindex; index++ )
    {
      int gridID0 = ncgrid[index].gridID;
      int gridtype0 = gridInqType(gridID0);
      if ( gridtype0 == GRID_CHARXY )
        {
          if ( gridInqXIsc(gridID0) == strlen &&
               gridInqXsize(gridID0) == dimlen )
            return;
          else if ( gridInqYIsc(gridID0) == strlen &&
               gridInqYsize(gridID0) == dimlen )
            return;
        }
    }

  int fileID  = streamptr->fileID;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

/* Define Dims */

  char dimname[CDI_MAX_NAME+3];
  dimname[0] = 0;
  if ( xory == 0 )
    cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, dimname);
  else
    cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, dimname);
  if ( dimname[0] == 0 ) { memcpy(dimname, "region", 7); dimname[6] = 0; }
  dimID = checkDimName(fileID, dimlen, dimname);
  if ( dimID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, dimlen, &dimID);

/* Define strlength dim */

  strcpy(dimname, "strlen");
  strlenID = checkDimName(fileID, strlen, dimname);
  if ( strlenID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, strlen, &strlenID);

/* Define Variable */

  int dimIDs[2];
  dimIDs[0] = dimID;
  dimIDs[1] = strlenID;

  char axisname[CDI_MAX_NAME]; axisname[0] = 0;
  char **cvals = (char **) Malloc(dimlen * sizeof(char *));
  for ( int i = 0; i < dimlen; i++ )
    cvals[i] = Malloc(strlen * sizeof(char) );
  int ncaxisid;
  if ( xory == 0 )
    {
      cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, axisname);
      gridInqXCvals(gridID, cvals);
    }
  else
    {
      cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, axisname);
      gridInqXCvals(gridID, cvals);
    }
  int status = nc_inq_varid(fileID, axisname, &ncaxisid);
  if ( status != NC_NOERR )
    {
      cdf_def_var(fileID, axisname, NC_CHAR, 2, dimIDs, &ncaxisid);
      if ( xory == 0 )
        cdfPutGridStdAtts(fileID, ncaxisid, gridID, 'X', &gridInqsX);
      else
        cdfPutGridStdAtts(fileID, ncaxisid, gridID, 'Y', &gridInqsY);
    }
  else
    return;
  cdf_enddef(fileID);

/* Write Var */

  size_t start[2], count[2];
  start[1] = 0;
  count[0] = 1;
  count[1] = strlen;
  for (int i = 0; i < dimlen; i++)
    {
      start[0] = i;
      status = nc_put_vara_text(fileID, ncaxisid, start, count, cvals[i]);
    }

  ncgrid[gridindex].gridID = gridID;
  if ( xory == 0 )
    {
      ncgrid[gridindex].ncIDs[CDF_DIMID_X] = dimID;
      ncgrid[gridindex].ncIDs[CDF_VARID_X] = ncaxisid;
    }
  else
    {
      ncgrid[gridindex].ncIDs[CDF_DIMID_Y] = dimID;
      ncgrid[gridindex].ncIDs[CDF_VARID_Y] = ncaxisid;
    }
  streamptr->ncmode = 2;
}

static
void cdfDefRgrid(stream_t *streamptr, int gridID, int gridindex)
{
  ncgrid_t *ncgrid = streamptr->ncgrid;

  size_t dimlen = (size_t)gridInqSize(gridID);

  int iz;
  int dimID;
  {
    struct idSearch search
      = cdfSearchIDBySize(0, (size_t)gridindex, ncgrid, CDF_DIMID_X,
                          GRID_GAUSSIAN_REDUCED, (int)dimlen,
                          gridInqType, gridInqSize);
    iz = search.numNonMatching;
    dimID = search.foundID;
  }

  if ( dimID == CDI_UNDEFID )
    {
      int fileID  = streamptr->fileID;
      static bool lwarn = true;
      if ( lwarn )
        {
          Warning("Creating a NetCDF file with data on a gaussian reduced grid.");
          Warning("The further processing of the resulting file is unsupported!");
          lwarn = false;
        }

      char axisname[16] = "rgridX";
      if ( iz == 0 ) axisname[5] = '\0';
      else           sprintf(&axisname[5], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_X] = dimID;
}

static
void cdfDefGdim(stream_t *streamptr, int gridID, int gridindex)
{
  ncgrid_t *ncgrid = streamptr->ncgrid;
  int iz = 0;
  int dimID = CDI_UNDEFID;

  size_t dimlen = (size_t)gridInqSize(gridID);

  if ( gridInqYsize(gridID) == 0 )
    {
      struct idSearch search
        = cdfSearchIDBySize(0, (size_t)gridindex, ncgrid, CDF_DIMID_X,
                            GRID_GENERIC, (int)dimlen,
                            gridInqType, gridInqSize);
      iz = search.numNonMatching;
      dimID = search.foundID;
    }

  if ( gridInqXsize(gridID) == 0 )
    {
      struct idSearch search
        = cdfSearchIDBySize(0, (size_t)gridindex, ncgrid, CDF_DIMID_Y,
                            GRID_GENERIC, (int)dimlen,
                            gridInqType, gridInqSize);
      iz += search.numNonMatching;
      dimID = search.foundID;
    }

  if ( dimID == CDI_UNDEFID )
    {
      int fileID  = streamptr->fileID;
      char dimname[CDI_MAX_NAME];
      strcpy(dimname, "gsize");

      dimID = checkDimName(fileID, dimlen, dimname);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      if ( dimID == CDI_UNDEFID ) cdf_def_dim(fileID, dimname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  ncgrid[gridindex].gridID = gridID;
  ncgrid[gridindex].ncIDs[CDF_DIMID_X] = dimID;
}

static
void cdfDefGrid(stream_t *streamptr, int gridID, int gridindex)
{
  if ( streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X] != CDI_UNDEFID ) return;

  int gridtype = gridInqType(gridID);
  int size     = gridInqSize(gridID);

  if ( CDI_Debug )
    Message("gridtype = %d  size = %d", gridtype, size);

  if ( CDI_reduce_dim && size == 1 )
    {
      // no grid information
      streamptr->ncgrid[gridindex].gridID = gridID;
      return;
    }

  if ( gridtype == GRID_GAUSSIAN    ||
       gridtype == GRID_LONLAT      ||
       gridtype == GRID_PROJECTION  ||
       gridtype == GRID_GENERIC )
    {
      if ( gridtype == GRID_GENERIC )
        {
          if ( size == 1 && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
            {
              // no grid information
              streamptr->ncgrid[gridindex].gridID = gridID;
            }
          else
            {
              bool lx = false, ly = false;
              if ( gridInqXsize(gridID) > 0 /*&& gridInqXvals(gridID, NULL) > 0*/ )
                {
                  cdfDefXaxis(streamptr, gridID, gridindex, 1);
                  lx = true;
                }

              if ( gridInqYsize(gridID) > 0 /*&& gridInqYvals(gridID, NULL) > 0*/ )
                {
                  cdfDefYaxis(streamptr, gridID, gridindex, 1);
                  ly = true;
                }

              if ( !lx && !ly ) cdfDefGdim(streamptr, gridID, gridindex);
            }
        }
      else
        {
          int ndims = !(gridtype == GRID_LONLAT && size == 1 && !gridInqHasDims(gridID));

          if ( gridInqXsize(gridID) > 0 ) cdfDefXaxis(streamptr, gridID, gridindex, ndims);
          if ( gridInqYsize(gridID) > 0 ) cdfDefYaxis(streamptr, gridID, gridindex, ndims);

          cdf_def_mapping(streamptr, gridID);
        }
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      cdfDefCurvilinear(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_UNSTRUCTURED )
    {
      cdfDefUnstructured(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_GAUSSIAN_REDUCED )
    {
      cdfDefRgrid(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_SPECTRAL )
    {
      cdfDefComplex(streamptr, gridID, gridindex);
      cdfDefSP(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_FOURIER )
    {
      cdfDefComplex(streamptr, gridID, gridindex);
      cdfDefFC(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_TRAJECTORY )
    {
      cdfDefTrajLon(streamptr, gridID, gridindex);
      cdfDefTrajLat(streamptr, gridID, gridindex);
    }
  else if ( gridtype == GRID_CHARXY )
    {
      int strlen = 0;
      if ( (strlen = gridInqXIsc(gridID)) )
        cdfDefCharacter(streamptr, gridID, gridindex, 0, strlen);
      else
        if ( gridInqXsize(gridID) > 0 ) cdfDefXaxis(streamptr, gridID, gridindex, 1);
      if ( (strlen = gridInqYIsc(gridID)) )
        cdfDefCharacter(streamptr, gridID, gridindex, 1, strlen);
      else
        if ( gridInqYsize(gridID) > 0 ) cdfDefYaxis(streamptr, gridID, gridindex, 1);
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
}


void cdfDefHistory(stream_t *streamptr, int size, const char *history)
{
  int ncid = streamptr->fileID;
  cdf_put_att_text(ncid, NC_GLOBAL, "history", (size_t) size, history);
}


void cdfDefVars(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  if ( vlistID == CDI_UNDEFID )
    Error("Internal problem! vlist undefined for streamptr %p", streamptr);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);

  int ngrids = vlistNgrids(vlistID);
  if ( 2*ngrids > MAX_GRIDS_PS ) Error("Internal problem! Too many grids per stream (max=%d)\n", MAX_GRIDS_PS);
  for ( int index = 0; index < 2*ngrids; ++index )
    {
      streamptr->ncgrid[index].gridID = CDI_UNDEFID;
      for (size_t i = 0; i < CDF_SIZE_ncIDs; ++i)
        streamptr->ncgrid[index].ncIDs[i] = CDI_UNDEFID;
    }

  for ( int index = 0; index < ngrids; ++index )
    {
      int gridID = vlistGrid(vlistID, index);
      cdfDefGrid(streamptr, gridID, index);
    }
  {
    int index = ngrids-1;
    for ( int i = 0; i < ngrids; ++i )
      {
        int gridID = vlistGrid(vlistID, i);
        int projID = gridInqProj(gridID);
        if ( projID != CDI_UNDEFID ) cdfDefGrid(streamptr, projID, ++index);
      }
  }
  int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; ++index )
    {
      int zaxisID = vlistZaxis(vlistID, index);
      if ( streamptr->zaxisID[index] == CDI_UNDEFID ) cdfDefZaxis(streamptr, zaxisID);
    }

  if ( streamptr->ncmode != 2 )
    {
      cdf_enddef(streamptr->fileID);
      streamptr->ncmode = 2;
    }
}

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
