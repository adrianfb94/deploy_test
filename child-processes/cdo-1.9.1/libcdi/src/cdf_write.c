#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"
#include "cdf.h"
#include "cdf_int.h"
#include "vlist.h"
#include "vlist_var.h"


void cdfDefVarDeflate(int ncid, int ncvarid, int deflate_level)
{
#if  defined(HAVE_NETCDF4)
  int retval;
  /* Set chunking, shuffle, and deflate. */
  int shuffle = 1;
  int deflate = 1;

  if ( deflate_level < 1 || deflate_level > 9 ) deflate_level = 1;

  if ( (retval = nc_def_var_deflate(ncid, ncvarid, shuffle, deflate, deflate_level)) )
    {
      Error("nc_def_var_deflate failed, status = %d", retval);
    }
#else
  static bool lwarn = true;
  if ( lwarn )
    {
      lwarn = false;
      Warning("Deflate compression failed, NetCDF4 not available!");
    }
#endif
}


int cdfDefDatatype(int datatype, int filetype)
{
  int xtype = NC_FLOAT;

  if ( datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64 )
    Error("CDI/NetCDF library does not support complex numbers!");

  if ( filetype == CDI_FILETYPE_NC4 )
    {
      if      ( datatype == CDI_DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == CDI_DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == CDI_DATATYPE_INT32  ) xtype = NC_INT;
#if  defined  (HAVE_NETCDF4)
      else if ( datatype == CDI_DATATYPE_UINT8  ) xtype = NC_UBYTE;
      else if ( datatype == CDI_DATATYPE_UINT16 ) xtype = NC_USHORT;
      else if ( datatype == CDI_DATATYPE_UINT32 ) xtype = NC_UINT;
#else
      else if ( datatype == CDI_DATATYPE_UINT8  ) xtype = NC_SHORT;
      else if ( datatype == CDI_DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == CDI_DATATYPE_UINT32 ) xtype = NC_INT;
#endif
      else if ( datatype == CDI_DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                        xtype = NC_FLOAT;
    }
  else
    {
      if      ( datatype == CDI_DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == CDI_DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == CDI_DATATYPE_INT32  ) xtype = NC_INT;
      else if ( datatype == CDI_DATATYPE_UINT8  ) xtype = NC_SHORT;
      else if ( datatype == CDI_DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == CDI_DATATYPE_UINT32 ) xtype = NC_INT;
      else if ( datatype == CDI_DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                        xtype = NC_FLOAT;
    }

  return xtype;
}

static
void cdfDefVarMissval(stream_t *streamptr, int varID, int dtype, int lcheck)
{
  if ( streamptr->vars[varID].defmiss == false )
    {
      int vlistID = streamptr->vlistID;
      int fileID  = streamptr->fileID;
      int ncvarid = streamptr->vars[varID].ncvarid;
      double missval = vlistInqVarMissval(vlistID, varID);

      if ( lcheck && streamptr->ncmode == 2 ) cdf_redef(fileID);

      int xtype = cdfDefDatatype(dtype, streamptr->filetype);

      if ( xtype == NC_BYTE && missval > 127 && missval < 256 ) xtype = NC_INT;

      if ( lcheck == 0 ||
           streamptr->ncmode != 2 ||
           streamptr->filetype == CDI_FILETYPE_NC ||
           streamptr->filetype == CDI_FILETYPE_NC2||
           streamptr->filetype == CDI_FILETYPE_NC5 )
        cdf_put_att_double(fileID, ncvarid, "_FillValue", (nc_type) xtype, 1, &missval);

      cdf_put_att_double(fileID, ncvarid, "missing_value", (nc_type) xtype, 1, &missval);

      if ( lcheck && streamptr->ncmode == 2 ) cdf_enddef(fileID);

      streamptr->vars[varID].defmiss = true;
    }
}

static
void cdfDefInstitut(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int instID  = vlistInqInstitut(vlistID);

  if ( instID != CDI_UNDEFID )
    {
      const char *longname = institutInqLongnamePtr(instID);
      if ( longname )
	{
	  size_t len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "institution", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
}

static
void cdfDefSource(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int modelID = vlistInqModel(vlistID);

  if ( modelID != CDI_UNDEFID )
    {
      const char *longname = modelInqNamePtr(modelID);
      if ( longname )
	{
          size_t len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "source", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
}

static inline
void *resizeBuf(void **buf, size_t *bufSize, size_t reqSize)
{
  if ( reqSize > *bufSize )
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}


void cdfDefineAttributes(int cdiID, int varID, int fileID, int ncvarID)
{
  int atttype, attlen;
  size_t len;
  char attname[CDI_MAX_NAME+1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for ( int iatt = 0; iatt < natts; ++iatt )
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      if ( attlen == 0 ) continue;

      if ( atttype == CDI_DATATYPE_TXT )
        {
          size_t attSize = (size_t)attlen*sizeof(char);
          char *atttxt = (char *)resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          len = (size_t)attlen;
          cdf_put_att_text(fileID, ncvarID, attname, len, atttxt);
        }
      else if ( atttype == CDI_DATATYPE_INT8  || atttype == CDI_DATATYPE_UINT8  ||
                atttype == CDI_DATATYPE_INT16 || atttype == CDI_DATATYPE_UINT16 ||
                atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32 )
        {
          size_t attSize = (size_t)attlen*sizeof(int);
          int *attint = (int *)resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          len = (size_t)attlen;
          nc_type xtype = (atttype == CDI_DATATYPE_INT8)  ? NC_BYTE :
                          (atttype == CDI_DATATYPE_INT16) ? NC_SHORT :
#if  defined  (HAVE_NETCDF4)
                          (atttype == CDI_DATATYPE_UINT8)  ? NC_UBYTE :
                          (atttype == CDI_DATATYPE_UINT16) ? NC_USHORT :
                          (atttype == CDI_DATATYPE_UINT32) ? NC_UINT :
#endif
                          NC_INT;
          cdf_put_att_int(fileID, ncvarID, attname, xtype, len, attint);
        }
      else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
        {
          size_t attSize = (size_t)attlen * sizeof(double);
          double *attflt = (double *)resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          len = (size_t)attlen;
          if ( atttype == CDI_DATATYPE_FLT32 )
            {
              float attflt_sp[len];
              for ( size_t i = 0; i < len; ++i ) attflt_sp[i] = (float)attflt[i];
              cdf_put_att_float(fileID, ncvarID, attname, NC_FLOAT, len, attflt_sp);
            }
          else
            cdf_put_att_double(fileID, ncvarID, attname, NC_DOUBLE, len, attflt);
        }
    }

  Free(attBuf);
}

static
void cdfDefGlobalAtts(stream_t *streamptr)
{
  if ( streamptr->globalatts ) return;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  cdfDefSource(streamptr);
  cdfDefInstitut(streamptr);

  int natts;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_redef(fileID);

  cdfDefineAttributes(vlistID, CDI_GLOBAL, fileID, NC_GLOBAL);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_enddef(fileID);

  streamptr->globalatts = 1;
}

static
void cdfDefLocalAtts(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  if ( streamptr->localatts ) return;
  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID ) return;

  streamptr->localatts = 1;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  for ( int varID = 0; varID < streamptr->nvars; varID++ )
    {
      int instID = vlistInqVarInstitut(vlistID, varID);
      if ( instID != CDI_UNDEFID )
	{
          int ncvarid = streamptr->vars[varID].ncvarid;
  	  const char *name = institutInqNamePtr(instID);
	  if ( name )
	    {
              size_t len = strlen(name);
	      cdf_put_att_text(fileID, ncvarid, "institution", len, name);
	    }
	}
      }

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
}

static
void cdf_get_gmapvarname(int gridID, char *gmapvarname)
{
  int pgridID = gridID;
  char mapping[CDI_MAX_NAME]; mapping[0] = 0;
  cdiGridInqKeyStr(pgridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);

  if ( !mapping[0] )
    {
      int projID = gridInqProj(gridID);
      if ( projID != CDI_UNDEFID )
        {
          pgridID = projID;
          cdiGridInqKeyStr(pgridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
        }
    }

  if ( mapping[0] )
    cdiGridInqKeyStr(pgridID, CDI_KEY_MAPPING, CDI_MAX_NAME, gmapvarname);
}

static
int nc_grid_index(stream_t *streamptr, int gridID)
{
  int index = 0;
  int vlistID = streamptr->vlistID;
  int ngrids = vlistNgrids(vlistID);
  for ( index = 0; index < ngrids; ++index )
    if ( streamptr->ncgrid[index].gridID == gridID ) break;

  assert(index < ngrids);

  return index;
}

static
int cdfDefVar(stream_t *streamptr, int varID)
{
  int ncvarid = -1;
  int xid = CDI_UNDEFID, yid = CDI_UNDEFID;
  size_t xsize = 0, ysize = 0;
  int dims[4];
  size_t chunks[4] = {0,0,0,0};
  int ndims = 0;
  int tablenum;
  int dimorder[3];
  size_t iax = 0;
  char axis[5];
  int ensID, ensCount, forecast_type;

  int fileID  = streamptr->fileID;

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d, varID = %d", streamptr->self, fileID, varID);

  if ( streamptr->vars[varID].ncvarid != CDI_UNDEFID )
    return streamptr->vars[varID].ncvarid;

  int vlistID  = streamptr->vlistID;
  int gridID   = vlistInqVarGrid(vlistID, varID);
  int zaxisID  = vlistInqVarZaxis(vlistID, varID);
  int timetype = vlistInqVarTimetype(vlistID, varID);
  int code     = vlistInqVarCode(vlistID, varID);
  int param    = vlistInqVarParam(vlistID, varID);
  int pnum, pcat, pdis;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  int chunktype = vlistInqVarChunkType(vlistID, varID);

  vlistInqVarDimorder(vlistID, varID, &dimorder);

  size_t gridsize  = (size_t)(gridInqSize(gridID));
  bool lchunk = (gridsize >= 16);
  int gridtype  = gridInqType(gridID);
  int gridindex = nc_grid_index(streamptr, gridID);
  if ( gridtype != GRID_TRAJECTORY )
    {
      xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
      if ( xid != CDI_UNDEFID ) cdf_inq_dimlen(fileID, xid, &xsize);
      if ( yid != CDI_UNDEFID ) cdf_inq_dimlen(fileID, yid, &ysize);
    }

  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];
  bool zaxis_is_scalar = false;
  if ( zid == CDI_UNDEFID ) zaxis_is_scalar = zaxisInqScalar(zaxisID) > 0;

  if ( dimorder[0] != 3 ) lchunk = false; /* ZYX and ZXY */

  if ( ((dimorder[0]>0)+(dimorder[1]>0)+(dimorder[2]>0)) < ((xid!=CDI_UNDEFID)+(yid!=CDI_UNDEFID)+(zid!=CDI_UNDEFID)) )
    {
      printf("zid=%d  yid=%d  xid=%d\n", zid, yid, xid);
      Error("Internal problem, dimension order missing!");
    }

  int tid = streamptr->basetime.ncdimid;

  if ( vlistHasTime(vlistID) && timetype != TIME_CONSTANT )
    {
      if ( tid == CDI_UNDEFID ) Error("Internal problem, time undefined!");
      chunks[ndims] = 1;
      dims[ndims++] = tid;
      axis[iax++] = 'T';
    }
  /*
  if ( zid != CDI_UNDEFID ) axis[iax++] = 'Z';
  if ( zid != CDI_UNDEFID ) chunks[ndims] = 1;
  if ( zid != CDI_UNDEFID ) dims[ndims++] = zid;

  if ( yid != CDI_UNDEFID ) chunks[ndims] = ysize;
  if ( yid != CDI_UNDEFID ) dims[ndims++] = yid;

  if ( xid != CDI_UNDEFID ) chunks[ndims] = xsize;
  if ( xid != CDI_UNDEFID ) dims[ndims++] = xid;
  */
  size_t chunk_size_max = 65536;
  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != CDI_UNDEFID )
        {
          axis[iax++] = 'Z';
          chunks[ndims] = 1;
          dims[ndims] = zid;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != CDI_UNDEFID )
        {
          if ( chunktype == CDI_CHUNK_AUTO )
            chunks[ndims] = (chunk_size_max > gridsize) ? ysize : chunk_size_max/xsize;
          else
            chunks[ndims] = (chunktype == CDI_CHUNK_LINES) ? 1 : ysize;
          dims[ndims] = yid;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != CDI_UNDEFID )
        {
          if ( chunktype == CDI_CHUNK_AUTO && yid == CDI_UNDEFID )
            chunks[ndims] = (chunk_size_max > xsize) ? xsize : chunk_size_max;
          else
            chunks[ndims] = xsize;
          dims[ndims] = xid;
          ndims++;
        }
    }

  if ( CDI_Debug )
    fprintf(stderr, "chunktype %d  chunks %d %d %d %d\n", chunktype, (int)chunks[0], (int)chunks[1], (int)chunks[2], (int)chunks[3]);

  char varname[CDI_MAX_NAME];
  char name[CDI_MAX_NAME]; name[0] = 0;
  char longname[CDI_MAX_NAME]; longname[0] = 0;
  char stdname[CDI_MAX_NAME]; stdname[0] = 0;
  char units[CDI_MAX_NAME]; units[0] = 0;
  if ( vlistInqVarNamePtr(vlistID, varID) ) vlistInqVarName(vlistID, varID, name);
  vlistInqVarLongname(vlistID, varID, longname);
  vlistInqVarStdname(vlistID, varID, stdname);
  vlistInqVarUnits(vlistID, varID, units);

  int tableID  = vlistInqVarTable(vlistID, varID);
  if ( !name[0] ) tableInqEntry(tableID, code, -1, name, longname, units);
  if ( name[0] )
    {
      sprintf(varname, "%s", name);

      bool checkname = true;
      int iz = 0;

      while ( checkname )
        {
          if ( iz ) sprintf(varname, "%s_%d", name, iz+1);

          int status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR ) checkname = false;

          if ( checkname ) iz++;

          if ( iz >= CDI_MAX_NAME ) Error("Double entry of variable name '%s'!", name);
        }

      if ( strcmp(name, varname) != 0 )
        {
          if ( iz == 1 )
            Warning("Changed double entry of variable name '%s' to '%s'!", name, varname);
          else
            Warning("Changed multiple entry of variable name '%s' to '%s'!", name, varname);
        }

      strcpy(name, varname);
    }
  else
    {
      if ( code < 0 ) code = -code;
      if ( pnum < 0 ) pnum = -pnum;

      if ( pdis == 255 )
	sprintf(varname, "var%d", code);
      else
	sprintf(varname, "param%d.%d.%d", pnum, pcat, pdis);

      char *varname2 = varname+strlen(varname);

      bool checkname = true;
      int iz = 0;

      while ( checkname )
        {
          if ( iz ) sprintf(varname2, "_%d", iz+1);

          int status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR ) checkname = false;

          if ( checkname ) iz++;

          if ( iz >= CDI_MAX_NAME ) break;
        }

      strcpy(name, varname);
      code = 0;
      pdis = 255;
    }

  /* if ( streamptr->ncmode == 2 ) cdf_redef(fileID); */

  int dtype = vlistInqVarDatatype(vlistID, varID);
  int xtype = cdfDefDatatype(dtype, streamptr->filetype);

  cdf_def_var(fileID, name, (nc_type) xtype, ndims, dims, &ncvarid);

#if  defined  (HAVE_NETCDF4)
  if ( lchunk && (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C) )
    cdf_def_var_chunking(fileID, ncvarid, NC_CHUNKED, chunks);
#endif

  if ( streamptr->comptype == CDI_COMPRESS_ZIP )
    {
      if ( lchunk && (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C) )
        {
          cdfDefVarDeflate(fileID, ncvarid, streamptr->complevel);
        }
      else
        {
          if ( lchunk )
            {
              static bool lwarn = true;
              if ( lwarn )
                {
                  lwarn = false;
                  Warning("Deflate compression is only available for NetCDF4!");
                }
            }
        }
    }

  if ( *stdname )  cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(stdname), stdname);
  if ( *longname ) cdf_put_att_text(fileID, ncvarid, "long_name", strlen(longname), longname);
  if ( *units )    cdf_put_att_text(fileID, ncvarid, "units", strlen(units), units);

  if ( code > 0 && pdis == 255 )
    cdf_put_att_int(fileID, ncvarid, "code", NC_INT, 1, &code);

  if ( pdis != 255 )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdf_put_att_text(fileID, ncvarid, "param", strlen(paramstr), paramstr);
    }

  if ( tableID != CDI_UNDEFID )
    {
      tablenum = tableInqNum(tableID);
      if ( tablenum > 0 )
        cdf_put_att_int(fileID, ncvarid, "table", NC_INT, 1, &tablenum);
    }

  char coordinates[CDI_MAX_NAME]; coordinates[0] = 0;

  if ( zaxis_is_scalar || zaxisInqType(zaxisID) == ZAXIS_CHAR )
    {
      int nczvarID = streamptr->nczvarID[zaxisindex];
      if ( nczvarID != CDI_UNDEFID )
        {
          size_t len = strlen(coordinates);
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, nczvarID, coordinates+len);
        }
    }

  if ( gridtype != GRID_GENERIC && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN &&
       gridtype != GRID_PROJECTION && gridtype != GRID_CURVILINEAR && gridtype != GRID_CHARXY )
    {
      size_t len = strlen(gridNamePtr(gridtype));
      if ( len > 0 )
        cdf_put_att_text(fileID, ncvarid, "CDI_grid_type", len, gridNamePtr(gridtype));
    }

  char gmapvarname[CDI_MAX_NAME]; gmapvarname[0] = 0;

  cdf_get_gmapvarname(gridID, gmapvarname);

  if ( gmapvarname[0] ) cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(gmapvarname), gmapvarname);

  if ( gridtype == GRID_TRAJECTORY )
    {
      cdf_put_att_text(fileID, ncvarid, "coordinates", 9, "lon lat" );
    }
  else if ( gridtype == GRID_LONLAT && xid == CDI_UNDEFID && yid == CDI_UNDEFID && gridsize == 1 )
    {
      int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      if ( ncyvarID != CDI_UNDEFID )
        {
          size_t len = strlen(coordinates);
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncyvarID, coordinates+len);
        }
      if ( ncxvarID != CDI_UNDEFID )
        {
          size_t len = strlen(coordinates);
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncxvarID, coordinates+len);
        }
    }
  else if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    {
      char cellarea[CDI_MAX_NAME] = "area: ";
      int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      int ncavarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_A];
      // CMOR order: coordinates = "lat lon"
      if ( cdiCoordinatesLonLat )
        {
          if ( ncxvarID != CDI_UNDEFID )
            {
              size_t len = strlen(coordinates);
              if ( len ) coordinates[len++] = ' ';
              cdf_inq_varname(fileID, ncxvarID, coordinates+len);
            }
          if ( ncyvarID != CDI_UNDEFID )
            {
              size_t len = strlen(coordinates);
              if ( len ) coordinates[len++] = ' ';
              cdf_inq_varname(fileID, ncyvarID, coordinates+len);
            }
        }
      else
        {
          if ( ncyvarID != CDI_UNDEFID )
            {
              size_t len = strlen(coordinates);
              if ( len ) coordinates[len++] = ' ';
              cdf_inq_varname(fileID, ncyvarID, coordinates+len);
            }
          if ( ncxvarID != CDI_UNDEFID )
            {
              size_t len = strlen(coordinates);
              if ( len ) coordinates[len++] = ' ';
              cdf_inq_varname(fileID, ncxvarID, coordinates+len);
            }
        }

      if ( ncavarID != CDI_UNDEFID )
        {
          size_t len = strlen(cellarea);
          cdf_inq_varname(fileID, ncavarID, cellarea+len);
          len = strlen(cellarea);
          cdf_put_att_text(fileID, ncvarid, "cell_measures", len, cellarea);
        }

      if ( gridtype == GRID_UNSTRUCTURED )
        {
          int position = gridInqPosition(gridID);
          if ( position > 0 )
            cdf_put_att_int(fileID, ncvarid, "number_of_grid_in_reference", NC_INT, 1, &position);
        }
    }
  else if ( gridtype == GRID_SPECTRAL || gridtype == GRID_FOURIER )
    {
      int gridTruncation = gridInqTrunc(gridID);
      axis[iax++] = '-';
      axis[iax++] = '-';
      cdf_put_att_text(fileID, ncvarid, "axis", iax, axis);
      cdf_put_att_int(fileID, ncvarid, "truncation", NC_INT, 1, &gridTruncation);
    }
  else if ( gridtype == GRID_CHARXY )
    {
      if ( gridInqXIsc(gridID) )
        {
          int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
          size_t len = strlen(coordinates);
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncxvarID, coordinates+len);
        }
      else if ( gridInqYIsc(gridID) )
        {
          int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
          size_t len = strlen(coordinates);
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncyvarID, coordinates+len);
        }
    }

  size_t len = strlen(coordinates);
  if ( len ) cdf_put_att_text(fileID, ncvarid, "coordinates", len, coordinates);

  /*  if ( xtype == NC_BYTE || xtype == NC_SHORT || xtype == NC_INT ) */
    {
      int astype = NC_DOUBLE;

      double addoffset   = vlistInqVarAddoffset(vlistID, varID);
      double scalefactor = vlistInqVarScalefactor(vlistID, varID);
      bool laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      bool lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      if ( laddoffset || lscalefactor )
        {
          if ( IS_EQUAL(addoffset,   (double) ((float) addoffset)) &&
               IS_EQUAL(scalefactor, (double) ((float) scalefactor)) )
            {
              astype = NC_FLOAT;
            }

          if ( xtype == (int) NC_FLOAT ) astype = NC_FLOAT;

          cdf_put_att_double(fileID, ncvarid, "add_offset",   (nc_type) astype, 1, &addoffset);
          cdf_put_att_double(fileID, ncvarid, "scale_factor", (nc_type) astype, 1, &scalefactor);
        }
    }

  if ( dtype == CDI_DATATYPE_UINT8 && xtype == NC_BYTE )
    {
      int validrange[2] = {0, 255};
      cdf_put_att_int(fileID, ncvarid, "valid_range", NC_SHORT, 2, validrange);
      cdf_put_att_text(fileID, ncvarid, "_Unsigned", 4, "true");
    }

  streamptr->vars[varID].ncvarid = ncvarid;

  if ( vlistInqVarMissvalUsed(vlistID, varID) )
    cdfDefVarMissval(streamptr, varID, vlistInqVarDatatype(vlistID, varID), 0);

  if ( zid == -1 )
    {
      if ( zaxisInqType(zaxisID) == ZAXIS_CLOUD_BASE          ||
           zaxisInqType(zaxisID) == ZAXIS_CLOUD_TOP           ||
           zaxisInqType(zaxisID) == ZAXIS_ISOTHERM_ZERO       ||
           zaxisInqType(zaxisID) == ZAXIS_TOA                 ||
           zaxisInqType(zaxisID) == ZAXIS_SEA_BOTTOM          ||
           zaxisInqType(zaxisID) == ZAXIS_LAKE_BOTTOM         ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM     ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM_TA  ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM_TW  ||
           zaxisInqType(zaxisID) == ZAXIS_MIX_LAYER           ||
           zaxisInqType(zaxisID) == ZAXIS_ATMOSPHERE )
        {
          zaxisInqName(zaxisID, varname);
          cdf_put_att_text(fileID, ncvarid, "level_type", strlen(varname), varname);
        }
    }

  if ( vlistInqVarEnsemble( vlistID,  varID, &ensID, &ensCount, &forecast_type ) )
    {
      /* void cdf_put_att_int(  int ncid, int varid, const char *name, nc_type xtype,
	                        size_t len, const int *ip )
       */
	cdf_put_att_int(fileID, ncvarid, "realization", NC_INT, 1, &ensID);
	cdf_put_att_int(fileID, ncvarid, "ensemble_members", NC_INT, 1, &ensCount);
	cdf_put_att_int(fileID, ncvarid, "forecast_init_type", NC_INT, 1, &forecast_type);
    }

  /* Attributes */
  cdfDefineAttributes(vlistID, varID, fileID, ncvarid);

  /* if ( streamptr->ncmode == 2 ) cdf_enddef(fileID); */

  return ncvarid;
}


void cdfEndDef(stream_t *streamptr)
{
  cdfDefGlobalAtts(streamptr);
  cdfDefLocalAtts(streamptr);

  if ( streamptr->accessmode == 0 )
    {
      int fileID  = streamptr->fileID;
      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      int nvars =  streamptr->nvars;
      for ( int varID = 0; varID < nvars; varID++ )
	cdfDefVar(streamptr, varID);

      if ( streamptr->ncmode == 2 )
        {
          if ( CDI_netcdf_hdr_pad == 0UL )
            cdf_enddef(fileID);
          else
            cdf__enddef(fileID, CDI_netcdf_hdr_pad);
        }

      streamptr->accessmode = 1;
    }
}

static
void cdfWriteGridTraj(stream_t *streamptr, int gridID)
{
  int gridindex = nc_grid_index(streamptr, gridID);
  int lonID = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
  int latID = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];

  double xlon = gridInqXval(gridID, 0);
  double xlat = gridInqYval(gridID, 0);
  int tsID = streamptr->curTsID;
  size_t index = (size_t)tsID;

  int fileID = streamptr->fileID;
  cdf_put_var1_double(fileID, lonID, &index, &xlon);
  cdf_put_var1_double(fileID, latID, &index, &xlat);
}

static
void cdf_write_var_data(int fileID, int vlistID, int varID, int ncvarid, int dtype, size_t nvals, size_t xsize, size_t ysize,
                        bool swapxy, size_t *start, size_t *count, int memtype, const void *data, int nmiss)
{
  const double *pdata_dp = (const double *) data;
  double *mdata_dp = NULL;
  double *sdata_dp = NULL;
  const float *pdata_sp = (const float *) data;
  float *mdata_sp = NULL;
  float *sdata_sp = NULL;

  /*  if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 ) */
    {
      double missval      = vlistInqVarMissval(vlistID, varID);
      double addoffset    = vlistInqVarAddoffset(vlistID, varID);
      double scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      bool laddoffset     = IS_NOT_EQUAL(addoffset, 0);
      bool lscalefactor   = IS_NOT_EQUAL(scalefactor, 1);

      if ( laddoffset || lscalefactor )
        {
          if ( memtype == MEMTYPE_FLOAT )
            {
              mdata_sp = (float *) Malloc(nvals*sizeof(float));
              memcpy(mdata_sp, pdata_sp, nvals*sizeof(float));
              pdata_sp = mdata_sp;

              if ( nmiss > 0 )
                {
                  for ( size_t i = 0; i < nvals; i++ )
                    {
                      double temp = mdata_sp[i];
                      if ( !DBL_IS_EQUAL(temp, missval) )
                        {
                          if ( laddoffset )   temp -= addoffset;
                          if ( lscalefactor ) temp /= scalefactor;
                          mdata_sp[i] = (float)temp;
                        }
                    }
                }
              else
                {
                  for ( size_t i = 0; i < nvals; i++ )
                    {
                      double temp = mdata_sp[i];
                      if ( laddoffset )   temp -= addoffset;
                      if ( lscalefactor ) temp /= scalefactor;
                      mdata_sp[i] = (float)temp;
                    }
                }
            }
          else
            {
              mdata_dp = (double *) Malloc(nvals*sizeof(double));
              memcpy(mdata_dp, pdata_dp, nvals*sizeof(double));
              pdata_dp = mdata_dp;

              if ( nmiss > 0 )
                {
                  for ( size_t i = 0; i < nvals; i++ )
                    {
                      if ( !DBL_IS_EQUAL(mdata_dp[i], missval) )
                        {
                          if ( laddoffset )   mdata_dp[i] -= addoffset;
                          if ( lscalefactor ) mdata_dp[i] /= scalefactor;
                        }
                    }
                }
              else
                {
                  for ( size_t i = 0; i < nvals; i++ )
                    {
                      if ( laddoffset )   mdata_dp[i] -= addoffset;
                      if ( lscalefactor ) mdata_dp[i] /= scalefactor;
                    }
                }
            }
        }

      if ( dtype == CDI_DATATYPE_UINT8 || dtype == CDI_DATATYPE_INT8 ||
           dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
        {
          if ( memtype == MEMTYPE_FLOAT )
            {
              if ( mdata_sp == NULL )
                {
                  mdata_sp = (float *) Malloc(nvals*sizeof(float));
                  memcpy(mdata_sp, pdata_sp, nvals*sizeof(float));
                  pdata_sp = mdata_sp;
                }

              for ( size_t i = 0; i < nvals; i++ ) mdata_sp[i] = roundf(mdata_sp[i]);

              if ( dtype == CDI_DATATYPE_UINT8 )
                {
                  nc_type xtype;
                  cdf_inq_vartype(fileID, ncvarid, &xtype);
                  if ( xtype == NC_BYTE )
                    {
                      for ( size_t i = 0; i < nvals; ++i )
                        if ( mdata_sp[i] > 127 ) mdata_sp[i] -= 256;
                    }
                }
            }
          else
            {
              if ( mdata_dp == NULL )
                {
                  mdata_dp = (double *) Malloc(nvals*sizeof(double));
                  memcpy(mdata_dp, pdata_dp, nvals*sizeof(double));
                  pdata_dp = mdata_dp;
                }

              for ( size_t i = 0; i < nvals; i++ ) mdata_dp[i] = round(mdata_dp[i]);

              if ( dtype == CDI_DATATYPE_UINT8 )
                {
                  nc_type xtype;
                  cdf_inq_vartype(fileID, ncvarid, &xtype);
                  if ( xtype == NC_BYTE )
                    {
                      for ( size_t i = 0; i < nvals; ++i )
                        if ( mdata_dp[i] > 127 ) mdata_dp[i] -= 256;
                    }
                }
            }
        }

      if ( CDF_Debug && memtype != MEMTYPE_FLOAT )
        {
          double fmin =  1.0e200;
          double fmax = -1.0e200;
          for ( size_t i = 0; i < nvals; ++i )
            {
              if ( !DBL_IS_EQUAL(pdata_dp[i], missval) )
                {
                  if ( pdata_dp[i] < fmin ) fmin = pdata_dp[i];
                  if ( pdata_dp[i] > fmax ) fmax = pdata_dp[i];
                }
            }
          Message("nvals = %zu, nmiss = %d, missval = %g, minval = %g, maxval = %g",
                  nvals, nmiss, missval, fmin, fmax);
        }
    }

  if ( swapxy ) // implemented only for cdf_write_var_slice()
    {
      size_t gridsize = xsize*ysize;
      if ( memtype == MEMTYPE_FLOAT )
        {
          sdata_sp = (float *) Malloc(gridsize*sizeof(float));
          for ( size_t j = 0; j < ysize; ++j )
            for ( size_t i = 0; i < xsize; ++i )
              sdata_sp[i*ysize+j] = pdata_sp[j*xsize+i];
          pdata_sp = sdata_sp;
        }
      else
        {
          sdata_dp = (double *) Malloc(gridsize*sizeof (double));
          for ( size_t j = 0; j < ysize; ++j )
            for ( size_t i = 0; i < xsize; ++i )
              sdata_dp[i*ysize+j] = pdata_dp[j*xsize+i];
          pdata_dp = sdata_dp;
        }
    }

  if ( memtype == MEMTYPE_FLOAT )
    cdf_put_vara_float(fileID, ncvarid, start, count, pdata_sp);
  else
    cdf_put_vara_double(fileID, ncvarid, start, count, pdata_dp);

  if ( mdata_dp ) Free(mdata_dp);
  if ( sdata_dp ) Free(sdata_dp);
  if ( mdata_sp ) Free(mdata_sp);
  if ( sdata_sp ) Free(sdata_sp);
}


void cdf_write_var(stream_t *streamptr, int varID, int memtype, const void *data, int nmiss)
{
  if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);

  size_t xsize = 0, ysize = 0;
  size_t size;
  size_t start[5];
  size_t count[5];
  bool swapxy = false;
  size_t ndims = 0;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %ld", ntsteps);

  int ncvarid = cdfDefVar(streamptr, varID);

  int gridID   = vlistInqVarGrid(vlistID, varID);
  int zaxisID  = vlistInqVarZaxis(vlistID, varID);
  int timetype = vlistInqVarTimetype(vlistID, varID);

  int xid = CDI_UNDEFID, yid = CDI_UNDEFID;
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      int gridindex = nc_grid_index(streamptr, gridID);
      xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
    }

  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  if ( timetype != TIME_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  if ( zid != CDI_UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = (size_t)zaxisInqSize(zaxisID);
      ndims++;
    }

  if ( yid != CDI_UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, yid, &size);
      /*      count[ndims] = gridInqYsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if ( xid != CDI_UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, xid, &size);
      /*      count[ndims] = gridInqXsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if ( CDI_Debug )
    for (size_t idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = (size_t)(gridInqSize(gridID)) * (size_t)(zaxisInqSize(zaxisID));

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}


void cdf_write_var_chunk(stream_t *streamptr, int varID, int memtype,
                         const int rect[][2], const void *data, int nmiss)
{
  if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);

  int xid = CDI_UNDEFID, yid = CDI_UNDEFID;
  size_t xsize = 0, ysize = 0;
  size_t start[5];
  size_t count[5];
  bool swapxy = false;
  size_t ndims = 0;
  int streamID = streamptr->self;

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  int vlistID = streamInqVlist(streamID);
  int fileID  = streamInqFileID(streamID);

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %ld", ntsteps);

  int ncvarid = cdfDefVar(streamptr, varID);

  int gridID   = vlistInqVarGrid(vlistID, varID);
  int zaxisID  = vlistInqVarZaxis(vlistID, varID);
  int timetype = vlistInqVarTimetype(vlistID, varID);

  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      int gridindex = nc_grid_index(streamptr, gridID);
      xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
    }

  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  if ( timetype != TIME_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }
  if ( zid != CDI_UNDEFID )
    {
      int size = zaxisInqSize(zaxisID);
      xassert(rect[2][0] >= 0 && rect[2][0] <= rect[2][1]
              && rect[2][1] <= size);
      start[ndims] = (size_t)rect[2][0];
      count[ndims] = (size_t)rect[2][1] - (size_t)rect[2][0] + 1;
      ndims++;
    }
  if ( yid != CDI_UNDEFID )
    {
      size_t size;
      cdf_inq_dimlen(fileID, yid, &size);
      xassert(rect[1][0] >= 0 && rect[1][0] <= rect[1][1]
              && (size_t)rect[1][1] <= size);
      start[ndims] = (size_t)rect[1][0];
      count[ndims] = (size_t)rect[1][1] - (size_t)rect[1][0] + 1;
      ndims++;
    }
  if ( xid != CDI_UNDEFID )
    {
      size_t size;
      cdf_inq_dimlen(fileID, xid, &size);
      xassert(rect[0][0] >= 0 && rect[0][0] <= rect[0][1]
              && (size_t)rect[0][1] <= size);
      start[ndims] = (size_t)rect[0][0];
      count[ndims] = (size_t)rect[0][1] - (size_t)rect[0][0] + 1;
      ndims++;
    }

  if ( CDI_Debug )
    for (size_t idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = (size_t)(gridInqSize(gridID)) * (size_t)(zaxisInqSize(zaxisID));

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals,
                     xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}


void cdf_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, int nmiss)
{
  if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);

  size_t xsize = 0, ysize = 0;
  size_t start[5];
  size_t count[5];
  int dimorder[3];
  int xid = CDI_UNDEFID, yid = CDI_UNDEFID;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %ld", ntsteps);

  int ncvarid = cdfDefVar(streamptr, varID);

  int gridID   = vlistInqVarGrid(vlistID, varID);
  int zaxisID  = vlistInqVarZaxis(vlistID, varID);
  int timetype = vlistInqVarTimetype(vlistID, varID);
  vlistInqVarDimorder(vlistID, varID, &dimorder);

  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      int gridindex = nc_grid_index(streamptr, gridID);
      xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
    }

  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  bool swapxy = (dimorder[2] == 2 || dimorder[0] == 1) && xid != CDI_UNDEFID && yid != CDI_UNDEFID;

  size_t ndims = 0;
  if ( timetype != TIME_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != CDI_UNDEFID )
        {
          start[ndims] = (size_t)levelID;
          count[ndims] = 1;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != CDI_UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, yid, &ysize);
          count[ndims] = ysize;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != CDI_UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, xid, &xsize);
          count[ndims] = xsize;
          ndims++;
        }
    }

  if ( CDI_Debug )
    for (size_t idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = (size_t)(gridInqSize(gridID));

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}


void cdf_write_record(stream_t *streamptr, int memtype, const void *data, int nmiss)
{
  int varID   = streamptr->record->varID;
  int levelID = streamptr->record->levelID;

  cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
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
