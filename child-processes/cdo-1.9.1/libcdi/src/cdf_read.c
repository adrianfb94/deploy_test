#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include <limits.h>
#include <float.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"
#include "cdf.h"
#include "cdf_int.h"
#include "vlist.h"
#include "vlist_var.h"


static
void cdfReadGridTraj(stream_t *streamptr, int gridID)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int gridindex = vlistGridIndex(vlistID, gridID);
  int lonID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
  int latID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];

  int tsID = streamptr->curTsID;
  size_t index = (size_t)tsID;

  double xlon, xlat;
  cdf_get_var1_double(fileID, lonID, &index, &xlon);
  cdf_get_var1_double(fileID, latID, &index, &xlat);

  gridDefXvals(gridID, &xlon);
  gridDefYvals(gridID, &xlat);
}

static
void cdfGetSlapDescription(stream_t *streamptr, int varID, size_t (*start)[4], size_t (*count)[4])
{
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int timetype = vlistInqVarTimetype(vlistID, varID);
  int gridindex = vlistGridIndex(vlistID, gridID);

  if ( CDI_Debug ) Message("tsID = %d", tsID);

  int xid = CDI_UNDEFID, yid = CDI_UNDEFID;
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfReadGridTraj(streamptr, gridID);
    }
  else
    {
      xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
    }
  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  int ndims = 0;
#define addDimension(startCoord, length) do \
    { \
      (*start)[ndims] = startCoord; \
      (*count)[ndims] = length; \
      ndims++; \
    } while(0)
  if ( timetype != TIME_CONSTANT ) addDimension((size_t)tsID, 1);
  if ( zid != CDI_UNDEFID ) addDimension(0, (size_t)zaxisInqSize(zaxisID));
  if ( yid != CDI_UNDEFID ) addDimension(0, (size_t)gridInqYsize(gridID));
  if ( xid != CDI_UNDEFID ) addDimension(0, (size_t)gridInqXsize(gridID));
#undef addDimension

  assert(ndims <= (int)(sizeof(*start)/sizeof(**start)));
  assert(ndims <= (int)(sizeof(*count)/sizeof(**count)));

  if ( CDI_Debug )
    for (int idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);
}

//Scans the data array for missVals, optionally applying first a scale factor and then an offset.
//Returns the number of missing + out-of-range values encountered.
static
size_t cdfDoInputDataTransformationDP(size_t valueCount, double *data, bool haveMissVal, double missVal, double scaleFactor, double offset, double validMin, double validMax)
 {
  const bool haveOffset   = IS_NOT_EQUAL(offset, 0);
  const bool haveScaleFactor = IS_NOT_EQUAL(scaleFactor, 1);
  size_t missValCount = 0;

  if ( IS_EQUAL(validMin, VALIDMISS) ) validMin = DBL_MIN;
  if ( IS_EQUAL(validMax, VALIDMISS) ) validMax = DBL_MAX;

  bool haveRangeCheck = (IS_NOT_EQUAL(validMax, DBL_MAX)) | (IS_NOT_EQUAL(validMin,DBL_MIN));
  assert(!haveRangeCheck || haveMissVal);

  switch ((int)haveMissVal | ((int)haveScaleFactor << 1)
          | ((int)haveOffset << 2) | ((int)haveRangeCheck << 3))
    {
    case 15: /* haveRangeCheck & haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] * scaleFactor + offset;
        }
      break;
    case 13: /* haveRangeCheck & haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] + offset;
        }
      break;
    case 11: /* haveRangeCheck & haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] * scaleFactor;
        }
      break;
    case 9: /* haveRangeCheck & haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal : data[i];
        }
      break;
    case 7: /* haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = data[i] * scaleFactor + offset;
      break;
    case 6: /* haveOffset & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = data[i] * scaleFactor + offset;
      break;
    case 5: /* haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] += offset;
      break;
    case 4: /* haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] += offset;
      break;
    case 3: /* haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] *= scaleFactor;
      break;
    case 2: /* haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] *= scaleFactor;
      break;
    case 1: /* haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        missValCount += (unsigned)DBL_IS_EQUAL(data[i], missVal);
      break;
    }

  return missValCount;
}

static
size_t cdfDoInputDataTransformationSP(size_t valueCount, float *data, bool haveMissVal, double missVal, double scaleFactor, double offset, double validMin, double validMax)
 {
  const bool haveOffset   = IS_NOT_EQUAL(offset, 0);
  const bool haveScaleFactor = IS_NOT_EQUAL(scaleFactor, 1);
  size_t missValCount = 0;

  if ( IS_EQUAL(validMin, VALIDMISS) ) validMin = DBL_MIN;
  if ( IS_EQUAL(validMax, VALIDMISS) ) validMax = DBL_MAX;

  bool haveRangeCheck = (IS_NOT_EQUAL(validMax, DBL_MAX)) | (IS_NOT_EQUAL(validMin,DBL_MIN));
  assert(!haveRangeCheck || haveMissVal);

  switch ((int)haveMissVal | ((int)haveScaleFactor << 1)
          | ((int)haveOffset << 2) | ((int)haveRangeCheck << 3))
    {
    case 15: /* haveRangeCheck & haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] * scaleFactor + offset);
        }
      break;
    case 13: /* haveRangeCheck & haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] + offset);
        }
      break;
    case 11: /* haveRangeCheck & haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] * scaleFactor);
        }
      break;
    case 9: /* haveRangeCheck & haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal : data[i];
        }
      break;
    case 7: /* haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] * scaleFactor + offset);
      break;
    case 6: /* haveOffset & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] * scaleFactor + offset);
      break;
    case 5: /* haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] + offset);
      break;
    case 4: /* haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] + offset);
      break;
    case 3: /* haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] * scaleFactor);
      break;
    case 2: /* haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] * scaleFactor);
      break;
    case 1: /* haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        missValCount += (unsigned)DBL_IS_EQUAL(data[i], missVal);
      break;
    }

  return missValCount;
}

static
size_t min_size(size_t a, size_t b)
{
  return a < b ? a : b;
}

static
void transpose2dArrayDP(size_t inWidth, size_t inHeight, double *data)
{
  const size_t cacheBlockSize = 256;   // Purely an optimization parameter. Current value of 32 means we are handling 8kB blocks,
                                       // which should be a decent compromise on many architectures.
#ifdef __cplusplus
  double *out[inHeight];
  double *temp[inWidth];
  temp[0] = (double *) Malloc(inHeight*inWidth*sizeof(double));
  memcpy(temp[0], data, inHeight*inWidth*sizeof(double));
  for(int i = 0; i < inHeight; i++) out[i] = data + (inWidth*i);
  for(int i = 1; i < inWidth; i++) temp[i] = temp[0] + (inHeight*i);
#else
  double (*out)[inHeight] = (double (*)[inHeight])data;
  double (*temp)[inWidth] = (double (*)[inWidth]) Malloc(inHeight*sizeof(*temp));
  memcpy(temp, data, inHeight*sizeof(*temp));
#endif

  /*
  for ( size_t y = 0; y < inHeight; ++y )
    for ( size_t x = 0; x < inWidth; ++x )
      out[x][y] = temp[y][x];
  */

  for ( size_t yBlock = 0; yBlock < inHeight; yBlock += cacheBlockSize )
    for ( size_t xBlock = 0; xBlock < inWidth; xBlock += cacheBlockSize )
      for ( size_t y = yBlock, yEnd = min_size(yBlock + cacheBlockSize, inHeight); y < yEnd; y++ )
        for ( size_t x = xBlock, xEnd = min_size(xBlock + cacheBlockSize, inWidth); x < xEnd; x++ )
          {
            out[x][y] = temp[y][x];
          }

  Free(temp[0]);
}

static
void transpose2dArraySP(size_t inWidth, size_t inHeight, float *data)
{
  const size_t cacheBlockSize = 256;   // Purely an optimization parameter. Current value of 32 means we are handling 8kB blocks,
                                       // which should be a decent compromise on many architectures.
#ifdef __cplusplus
  float *out[inHeight];
  float *temp[inWidth];
  temp[0] = (float *) Malloc(inHeight*inWidth*sizeof(float));
  memcpy(temp[0], data, inHeight*inWidth*sizeof(float));
  for(int i = 0; i < inHeight; i++) out[i] = data + (inWidth*i);
  for(int i = 1; i < inWidth; i++) temp[i] = temp[0] + (inHeight*i);
#else
  float (*out)[inHeight] = (float (*)[inHeight])data;
  float (*temp)[inWidth] = (float (*)[inWidth]) Malloc(inHeight*sizeof(*temp));
  memcpy(temp, data, inHeight*sizeof(*temp));
#endif

  /*
  for ( size_t y = 0; y < inHeight; ++y )
    for ( size_t x = 0; x < inWidth; ++x )
      out[x][y] = temp[y][x];
  */

  for ( size_t yBlock = 0; yBlock < inHeight; yBlock += cacheBlockSize )
    for ( size_t xBlock = 0; xBlock < inWidth; xBlock += cacheBlockSize )
      for ( size_t y = yBlock, yEnd = min_size(yBlock + cacheBlockSize, inHeight); y < yEnd; y++ )
        for ( size_t x = xBlock, xEnd = min_size(xBlock + cacheBlockSize, inWidth); x < xEnd; x++ )
          {
            out[x][y] = temp[y][x];
          }

  Free(temp);
}

static
void cdfInqDimIds(stream_t *streamptr, int varId, int (*outDimIds)[3])
{
  int gridId = vlistInqVarGrid(streamptr->vlistID, varId);
  int gridindex = vlistGridIndex(streamptr->vlistID, gridId);

  (*outDimIds)[0] = (*outDimIds)[1] = (*outDimIds)[2] = CDI_UNDEFID;
  switch ( gridInqType(gridId) )
    {
      case GRID_TRAJECTORY:
        cdfReadGridTraj(streamptr, gridId);
        break;

      case GRID_UNSTRUCTURED:
        (*outDimIds)[0] = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
        break;

      default:
        (*outDimIds)[0] = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
        (*outDimIds)[1] = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
        break;
    }

  int zaxisID = vlistInqVarZaxis(streamptr->vlistID, varId);
  int zaxisindex = vlistZaxisIndex(streamptr->vlistID, zaxisID);
  (*outDimIds)[2] = streamptr->zaxisID[zaxisindex];
}

static
int cdfGetSkipDim(int fileId, int ncvarid, int (*dimIds)[3])
{
  if((*dimIds)[0] != CDI_UNDEFID) return 0;
  if((*dimIds)[1] != CDI_UNDEFID) return 0;
  int nvdims;
  cdf_inq_varndims(fileId, ncvarid, &nvdims);
  if(nvdims != 3) return 0;

  int varDimIds[3];
  cdf_inq_vardimid(fileId, ncvarid, varDimIds);
  size_t size = 0;
  if ( (*dimIds)[2] == varDimIds[2] )
    {
      cdf_inq_dimlen(fileId, varDimIds[1], &size);
      if ( size == 1 ) return 1;
    }
  else if ( (*dimIds)[2] == varDimIds[1] )
    {
      cdf_inq_dimlen(fileId, varDimIds[2], &size);
      if ( size == 1 ) return 2;
    }
  return 0;
}

static
void cdfGetSliceSlapDescription(stream_t *streamptr, int varId, int levelId, bool *outSwapXY, size_t (*start)[4], size_t (*count)[4])
{
  int tsID = streamptr->curTsID;
  if ( CDI_Debug ) Message("tsID = %d", tsID);

  int fileId = streamptr->fileID;
  int vlistId = streamptr->vlistID;
  int ncvarid = streamptr->vars[varId].ncvarid;

  int gridId = vlistInqVarGrid(vlistId, varId);
  int timetype = vlistInqVarTimetype(vlistId, varId);
  int gridsize = gridInqSize(gridId);

  streamptr->numvals += gridsize;

  int dimIds[3];    //this array joins the old variables xid, yid, and zid
  cdfInqDimIds(streamptr, varId, &dimIds);

  int skipdim = cdfGetSkipDim(fileId, ncvarid, &dimIds);

  int dimorder[3];
  vlistInqVarDimorder(vlistId, varId, &dimorder);

  *outSwapXY = (dimorder[2] == 2 || dimorder[0] == 1) && dimIds[0] != CDI_UNDEFID && dimIds[1] != CDI_UNDEFID ;

  int ndims = 0;

#define addDimension(startIndex, extent) do {   \
      (*start)[ndims] = startIndex; \
      (*count)[ndims] = extent; \
      ndims++; \
  } while(0)

  if ( timetype != TIME_CONSTANT ) addDimension((size_t)tsID, 1);
  if ( skipdim == 1 ) addDimension(0, 1);

  for ( int id = 0; id < 3; ++id )
    {
      size_t size;
      int curDimId = dimIds[dimorder[id]-1];
      if ( curDimId == CDI_UNDEFID ) continue;
      switch ( dimorder[id] )
        {
          case 1:
          case 2:
            cdf_inq_dimlen(fileId, curDimId, &size);
            addDimension(0, size);
            break;

          case 3:
            addDimension((size_t)levelId, 1);
            break;

          default:
            Error("Internal errror: Malformed dimension order encountered. Please report this bug.\n");
        }
    }

  if ( skipdim == 2 ) addDimension(0, 1);

  assert(ndims <= (int)(sizeof(*start)/sizeof(**start)));
  assert(ndims <= (int)(sizeof(*count)/sizeof(**count)));

#undef addDimension

  if ( CDI_Debug )
    for (int idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, (*start)[idim], (*count)[idim]);

  int nvdims;
  cdf_inq_varndims(fileId, ncvarid, &nvdims);

  if ( nvdims != ndims )
    Error("Internal error, variable %s has an unsupported array structure!", vlistInqVarNamePtr(vlistId, varId));
}

static
void cdfReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int ncvarid = streamptr->vars[varID].ncvarid;

  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  size_t start[4];
  size_t count[4];
  cdfGetSlapDescription(streamptr, varID, &start, &count);

  cdf_get_vara_double(fileID, ncvarid, start, count, data);

  size_t size = (size_t)gridInqSize(gridID)*(size_t)zaxisInqSize(zaxisID);
  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationDP(size, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}

static
void cdfReadVarSP(stream_t *streamptr, int varID, float *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int ncvarid = streamptr->vars[varID].ncvarid;

  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  size_t start[4];
  size_t count[4];
  cdfGetSlapDescription(streamptr, varID, &start, &count);

  cdf_get_vara_float(fileID, ncvarid, start, count, data);

  size_t size = (size_t)gridInqSize(gridID) * (size_t)zaxisInqSize(zaxisID);
  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationSP(size, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}


void cdf_read_var(stream_t *streamptr, int varID, int memtype, void *data, int *nmiss)
{
  if ( memtype == MEMTYPE_DOUBLE )
    cdfReadVarDP(streamptr, varID, (double*) data, nmiss);
  else
    cdfReadVarSP(streamptr, varID, (float*) data, nmiss);
}

static
void cdfReadVarSliceDP(stream_t *streamptr, int varID, int levelID, double *data, int *nmiss)
{
  if ( CDI_Debug )
    Message("streamID = %d  varID = %d  levelID = %d", streamptr->self, varID, levelID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  size_t start[4], count[4];
  bool swapxy;
  cdfGetSliceSlapDescription(streamptr, varID, levelID, &swapxy, &start, &count);

  int ncvarid = streamptr->vars[varID].ncvarid;
  int gridId = vlistInqVarGrid(vlistID, varID);
  size_t gridsize = (size_t)gridInqSize(gridId);
  size_t xsize = (size_t)gridInqXsize(gridId);
  size_t ysize = (size_t)gridInqYsize(gridId);

  if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_FLT32 )
    {
      float *data_fp = (float *) Malloc(gridsize*sizeof(*data_fp));
      cdf_get_vara_float(fileID, ncvarid, start, count, data_fp);
      for ( size_t i = 0; i < gridsize; i++ )
        data[i] = (double) data_fp[i];
      Free(data_fp);
    }
  else
    {
      cdf_get_vara_double(fileID, ncvarid, start, count, data);
      
      if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_UINT8 )
        {
          nc_type xtype;
          cdf_inq_vartype(fileID, ncvarid, &xtype);
          if ( xtype == NC_BYTE )
            {
              for ( size_t i = 0; i < gridsize; i++ )
                if ( data[i] < 0 ) data[i] += 256;
            }
        }
    }

  if ( swapxy ) transpose2dArrayDP(ysize, xsize, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationDP(gridsize, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}

static
void cdfReadVarSliceSP(stream_t *streamptr, int varID, int levelID, float *data, int *nmiss)
{
  if ( CDI_Debug )
    Message("streamID = %d  varID = %d  levelID = %d", streamptr->self, varID, levelID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  size_t start[4], count[4];
  bool swapxy;
  cdfGetSliceSlapDescription(streamptr, varID, levelID, &swapxy, &start, &count);

  int ncvarid = streamptr->vars[varID].ncvarid;
  int gridId = vlistInqVarGrid(vlistID, varID);
  size_t gridsize = (size_t)gridInqSize(gridId);
  size_t xsize = (size_t)gridInqXsize(gridId);
  size_t ysize = (size_t)gridInqYsize(gridId);

  if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_FLT64 )
    {
      double *data_dp = (double *) Malloc(gridsize*sizeof(*data_dp));
      cdf_get_vara_double(fileID, ncvarid, start, count, data_dp);
      for ( size_t i = 0; i < gridsize; i++ )
        data[i] = (float) data_dp[i];
      Free(data_dp);
    }
  else
    {
      cdf_get_vara_float(fileID, ncvarid, start, count, data);

      if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_UINT8 )
        {
          nc_type xtype;
          cdf_inq_vartype(fileID, ncvarid, &xtype);
          if ( xtype == NC_BYTE )
            {
              for ( size_t i = 0; i < gridsize; i++ )
                if ( data[i] < 0 ) data[i] += 256;
            }
        }
    }

  if ( swapxy ) transpose2dArraySP(ysize, xsize, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationSP(gridsize, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}


void cdf_read_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, void *data, int *nmiss)
{
  if ( memtype == MEMTYPE_DOUBLE )
    cdfReadVarSliceDP(streamptr, varID, levelID, (double*) data, nmiss);
  else
    cdfReadVarSliceSP(streamptr, varID, levelID, (float*) data, nmiss);
}


void cdf_read_record(stream_t *streamptr, int memtype, void *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d", streamptr->self);

  int tsID    = streamptr->curTsID;
  int vrecID  = streamptr->tsteps[tsID].curRecID;
  int recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID   = streamptr->tsteps[tsID].records[recID].varID;
  int levelID = streamptr->tsteps[tsID].records[recID].levelID;

  if ( memtype == MEMTYPE_DOUBLE )
    cdfReadVarSliceDP(streamptr, varID, levelID, (double*) data, nmiss);
  else
    cdfReadVarSliceSP(streamptr, varID, levelID, (float*) data, nmiss);
}

#endif
