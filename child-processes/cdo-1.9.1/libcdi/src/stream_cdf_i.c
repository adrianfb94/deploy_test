#if defined (HAVE_CONFIG_H)
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

//#define TEST_GROUPS 1

#include <ctype.h>
#include <limits.h>

#include "dmemory.h"
#include "gaussgrid.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "stream_cdf.h"
#include "cdf_int.h"
#include "varscan.h"
#include "vlist.h"
#include "cdf_util.h"
#include "cdf_lazy_grid.h"


#define  X_AXIS  1
#define  Y_AXIS  2
#define  Z_AXIS  3
#define  T_AXIS  4

#define  POSITIVE_UP    1
#define  POSITIVE_DOWN  2

typedef struct {
  int     ncvarid;
  int     dimtype;
  size_t  len;
  char    name[CDI_MAX_NAME];
}
ncdim_t;
#define  MAX_COORDVARS  4
#define  MAX_AUXVARS    4

typedef struct {
  int      ncid;
  int      isvar;
  bool     ignore;
  bool     isx;
  bool     isy;
  bool     isc;
  bool     islon;
  bool     islat;
  bool     islev;
  bool     istime;
  bool     warn;
  bool     calendar;
  bool     climatology;
  bool     lformulaterms;
  int      param;
  int      code;
  int      tabnum;
  int      bounds;
  int      gridID;
  int      zaxisID;
  int      timetype;
  int      gridtype;
  int      zaxistype;
  int      xdim;
  int      ydim;
  int      zdim;
  int      xvarid;
  int      yvarid;
  int      zvarid;
  int      cvarids[MAX_COORDVARS];
  int      tvarid;
  int      psvarid;
  int      p0varid;
  int      ncoordvars;
  int      coordvarids[MAX_COORDVARS];
  int      nauxvars;
  int      auxvarids[MAX_AUXVARS];
  int      cellarea;
  int      tableID;
  int      truncation;
  int      position;
  bool     defmissval;
  bool     deffillval;
  int      xtype;
  int      gmapid;
  int      positive;
  int      ndims;
  int      dimids[8];
  int      dimtype[8];
  int      chunks[8];
  int      chunked;
  int      chunktype;
  int      natts;
  int      deflate;
  bool     lunsigned;
  bool     lvalidrange;
  int     *atts;
  size_t   vctsize;
  double  *vct;
  double   missval;
  double   fillval;
  double   addoffset;
  double   scalefactor;
  double   validrange[2];
  char     name[CDI_MAX_NAME];
  char     longname[CDI_MAX_NAME];
  char     stdname[CDI_MAX_NAME];
  char     units[CDI_MAX_NAME];
  char     extra[CDI_MAX_NAME];
  ensinfo_t   *ensdata;    /* Ensemble information */
}
ncvar_t;


static
void scanTimeString(const char *ptu, int *rdate, int *rtime)
{
  int year = 1, month = 1, day = 1;
  int hour = 0, minute = 0, second = 0;
  int v1 = 1, v2 = 1, v3 = 1;

  *rdate = 0;
  *rtime = 0;

  if ( *ptu )
    {
      v1 = atoi(ptu);
      if ( v1 < 0 ) ptu++;
      while ( isdigit((int) *ptu) ) ptu++;
      if ( *ptu )
        {
          v2 = atoi(++ptu);
          while ( isdigit((int) *ptu) ) ptu++;
          if ( *ptu )
            {
              v3 = atoi(++ptu);
              while ( isdigit((int) *ptu) ) ptu++;
            }
        }
    }

  if ( v3 > 999 && v1 < 32 )
    { year = v3; month = v2; day = v1; }
  else
    { year = v1; month = v2; day = v3; }

  while ( isspace((int) *ptu) ) ptu++;

  if ( *ptu )
    {
      while ( ! isdigit((int) *ptu) ) ptu++;

      hour = atoi(ptu);
      while ( isdigit((int) *ptu) ) ptu++;
      if ( *ptu == ':' )
        {
          ptu++;
          minute = atoi(ptu);
          while ( isdigit((int) *ptu) ) ptu++;
          if ( *ptu == ':' )
            {
              ptu++;
              second = atoi(ptu);
            }
        }
    }

  *rdate = cdiEncodeDate(year, month, day);
  *rtime = cdiEncodeTime(hour, minute, second);
}

static
int scanTimeUnit(const char *unitstr)
{
  size_t len = strlen(unitstr);
  int timeunit = get_timeunit(len, unitstr);
  if ( timeunit == -1 )
    Message("Unsupported TIMEUNIT: %s!", unitstr);

  return timeunit;
}

static
void setForecastTime(const char *timestr, taxis_t *taxis)
{
  size_t len = strlen(timestr);
  if ( len != 0 )
    scanTimeString(timestr, &taxis->fdate, &taxis->ftime);
  else
    taxis->fdate = taxis->ftime = 0;
}

static
int setBaseTime(const char *timeunits, taxis_t *taxis)
{
  int taxistype = TAXIS_ABSOLUTE;
  int rdate = -1, rtime = -1;

  size_t len = strlen(timeunits);
  while ( isspace(*timeunits) && len ) { timeunits++; len--; }

  char *restrict tu = (char *)Malloc((len+1) * sizeof(char));

  for ( size_t i = 0; i < len; i++ ) tu[i] = (char)tolower((int)timeunits[i]);
  tu[len] = 0;

  int timeunit = get_timeunit(len, tu);
  if ( timeunit == -1 )
    {
      Message("Unsupported TIMEUNIT: %s!", timeunits);
      return 1;
    }

  size_t pos = 0;
  while ( pos < len && !isspace(tu[pos]) ) ++pos;
  if ( tu[pos] )
    {
      while ( isspace(tu[pos]) ) ++pos;

      if ( str_is_equal(tu+pos, "since") )
        taxistype = TAXIS_RELATIVE;

      while ( pos < len && !isspace(tu[pos]) ) ++pos;
      if ( tu[pos] )
        {
          while ( isspace(tu[pos]) ) ++pos;

          if ( taxistype == TAXIS_ABSOLUTE )
            {
              if ( timeunit == TUNIT_DAY )
                {
                  if ( !str_is_equal(tu+pos, "%y%m%d.%f") )
                    {
                      Message("Unsupported format %s for TIMEUNIT day!", tu+pos);
                      timeunit = -1;
                    }
                }
              else if ( timeunit == TUNIT_MONTH )
                {
                  if ( !str_is_equal(tu+pos, "%y%m.%f") )
                    {
                      Message("Unsupported format %s for TIMEUNIT month!", tu+pos);
                      timeunit = -1;
                    }
                }
            }
          else if ( taxistype == TAXIS_RELATIVE )
            {
              scanTimeString(tu+pos, &rdate, &rtime);

              taxis->rdate = rdate;
              taxis->rtime = rtime;

              if ( CDI_Debug )
                Message("rdate = %d  rtime = %d", rdate, rtime);
            }
        }
    }

  taxis->type = taxistype;
  taxis->unit = timeunit;

  Free(tu);

  if ( CDI_Debug )
    Message("taxistype = %d  unit = %d", taxistype, timeunit);

  return 0;
}

static
bool xtypeIsText(int xtype)
{
  bool isText = ( xtype == NC_CHAR )
#if  defined  (HAVE_NETCDF4)
    || ( xtype == NC_STRING )
#endif
    ;
  return isText;
}

static
bool xtypeIsFloat(nc_type xtype)
{
  bool isFloat = xtype == NC_FLOAT || xtype == NC_DOUBLE;

  return isFloat;
}

static
bool xtypeIsInt(nc_type xtype)
{
  bool isInt = xtype == NC_SHORT || xtype == NC_INT
            || xtype == NC_BYTE
#if  defined  (HAVE_NETCDF4)
            || xtype == NC_USHORT || xtype == NC_UINT
            || xtype == NC_UBYTE
#endif
             ;

  return isInt;
}

static
int cdfInqDatatype(int xtype, bool lunsigned)
{
  int datatype = -1;

#if  defined  (HAVE_NETCDF4)
  if ( xtype == NC_BYTE && lunsigned ) xtype = NC_UBYTE;
#endif

  if      ( xtype == NC_BYTE   )  datatype = CDI_DATATYPE_INT8;
  else if ( xtype == NC_CHAR   )  datatype = CDI_DATATYPE_UINT8;
  else if ( xtype == NC_SHORT  )  datatype = CDI_DATATYPE_INT16;
  else if ( xtype == NC_INT    )  datatype = CDI_DATATYPE_INT32;
  else if ( xtype == NC_FLOAT  )  datatype = CDI_DATATYPE_FLT32;
  else if ( xtype == NC_DOUBLE )  datatype = CDI_DATATYPE_FLT64;
#if  defined  (HAVE_NETCDF4)
  else if ( xtype == NC_UBYTE  )  datatype = CDI_DATATYPE_UINT8;
  else if ( xtype == NC_LONG   )  datatype = CDI_DATATYPE_INT32;
  else if ( xtype == NC_USHORT )  datatype = CDI_DATATYPE_UINT16;
  else if ( xtype == NC_UINT   )  datatype = CDI_DATATYPE_UINT32;
  else if ( xtype == NC_INT64  )  datatype = CDI_DATATYPE_FLT64;
  else if ( xtype == NC_UINT64 )  datatype = CDI_DATATYPE_FLT64;
#endif

  return datatype;
}

static
void cdfGetAttInt(int fileID, int ncvarid, const char *attname, size_t attlen, int *attint)
{
  *attint = 0;

  nc_type atttype;
  size_t nc_attlen;
  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( xtypeIsFloat(atttype) || xtypeIsInt(atttype) )
    {
      bool lalloc = nc_attlen > attlen;
      int *pintatt = lalloc ? (int *)(Malloc(nc_attlen*sizeof(int))) : attint;
      cdf_get_att_int(fileID, ncvarid, attname, pintatt);
      if ( lalloc )
        {
          memcpy(attint, pintatt, attlen*sizeof(int));
          Free(pintatt);
        }
    }
}

static
void cdfGetAttDouble(int fileID, int ncvarid, char *attname, size_t attlen, double *attdouble)
{
  *attdouble = 0;

  nc_type atttype;
  size_t nc_attlen;
  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( xtypeIsFloat(atttype) || xtypeIsInt(atttype) )
    {
      bool lalloc = nc_attlen > attlen;
      double *pdoubleatt = lalloc ? (double*)Malloc(nc_attlen*sizeof(double)) : attdouble;
      cdf_get_att_double(fileID, ncvarid, attname, pdoubleatt);
      if ( lalloc )
        {
          memcpy(attdouble, pdoubleatt, attlen*sizeof(double));
          Free(pdoubleatt);
        }
    }
}

static
bool cdfCheckAttText(int fileID, int ncvarid, const char *attname)
{
  bool status = false;
  nc_type atttype;

  int status_nc = nc_inq_atttype(fileID, ncvarid, attname, &atttype);

  if ( status_nc == NC_NOERR
       && (atttype == NC_CHAR
#if  defined  (HAVE_NETCDF4)
           || atttype == NC_STRING
#endif
           ) )
    {
      status = true;
    }

  return status;
}

static
void cdfGetAttText(int fileID, int ncvarid, const char *attname, size_t attlen, char *atttext)
{
  nc_type atttype;
  size_t nc_attlen;

  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( atttype == NC_CHAR )
    {
      char attbuf[65636];
      if ( nc_attlen < sizeof(attbuf) )
        {
          cdf_get_att_text(fileID, ncvarid, attname, attbuf);

          if ( nc_attlen > (attlen-1) ) nc_attlen = (attlen-1);

          attbuf[nc_attlen++] = 0;
          memcpy(atttext, attbuf, nc_attlen);
        }
      else
        {
          atttext[0] = 0;
        }
    }
#if  defined  (HAVE_NETCDF4)
  else if ( atttype == NC_STRING )
    {
      if ( nc_attlen == 1 )
        {
          char *attbuf = NULL;
          cdf_get_att_string(fileID, ncvarid, attname, &attbuf);

          size_t ssize = strlen(attbuf) + 1;

          if ( ssize > attlen ) ssize = attlen;
          memcpy(atttext, attbuf, ssize);
          atttext[ssize - 1] = 0;
          Free(attbuf);
        }
      else
        {
          atttext[0] = 0;
        }
    }
#endif
}


void cdf_scale_add(size_t size, double *data, double addoffset, double scalefactor)
{
  bool laddoffset   = IS_NOT_EQUAL(addoffset, 0);
  bool lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

  if ( laddoffset && lscalefactor )
    {
      for (size_t i = 0; i < size; ++i )
        data[i] = data[i] * scalefactor + addoffset;
    }
  else if (lscalefactor)
    {
      for (size_t i = 0; i < size; ++i )
        data[i] *= scalefactor;
    }
  else if (laddoffset)
    {
      for (size_t i = 0; i < size; ++i )
        data[i] += addoffset;
    }
}

static
void cdfCreateRecords(stream_t *streamptr, int tsID)
{
  if ( tsID < 0 || (tsID >= streamptr->ntsteps && tsID > 0) ) return;

  if ( streamptr->tsteps[tsID].nallrecs > 0 ) return;

  int vlistID  = streamptr->vlistID;

  tsteps_t* sourceTstep = streamptr->tsteps;
  tsteps_t* destTstep = sourceTstep + tsID;

  int nvars = vlistNvars(vlistID);
  int nrecs = vlistNrecs(vlistID);

  if ( nrecs <= 0 ) return;

  if ( tsID == 0 )
    {
      int nvrecs = nrecs; /* use all records at first timestep */

      streamptr->nrecs += nrecs;

      destTstep->records    = (record_t *) Malloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = CDI_UNDEFID;
      destTstep->recIDs     = (int *) Malloc((size_t)nvrecs*sizeof (int));;
      for ( int recID = 0; recID < nvrecs; recID++ ) destTstep->recIDs[recID] = recID;

      record_t *records = destTstep->records;

      for ( int varID = 0, recID = 0; varID < nvars; varID++ )
        {
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          int nlev    = zaxisInqSize(zaxisID);
          for ( int levelID = 0; levelID < nlev; levelID++ )
            {
              recordInitEntry(&records[recID]);
              records[recID].varID   = (short)varID;
              records[recID].levelID = (short)levelID;
              recID++;
            }
        }
    }
  else if ( tsID == 1 )
    {
      int nvrecs = 0;
      for ( int varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT )
            {
              int zaxisID = vlistInqVarZaxis(vlistID, varID);
              nvrecs += zaxisInqSize(zaxisID);
            }
        }

      streamptr->nrecs += nvrecs;

      destTstep->records    = (record_t *) Malloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nvrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = CDI_UNDEFID;

      memcpy(destTstep->records, sourceTstep->records, (size_t)nrecs*sizeof(record_t));

      if ( nvrecs )
        {
          destTstep->recIDs = (int *) Malloc((size_t)nvrecs * sizeof (int));
          for ( int recID = 0, vrecID = 0; recID < nrecs; recID++ )
            {
              int varID = destTstep->records[recID].varID;
              if ( vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT )
                {
                  destTstep->recIDs[vrecID++] = recID;
                }
            }
        }
    }
  else
    {
      if ( streamptr->tsteps[1].records == 0 ) cdfCreateRecords(streamptr, 1);

      int nvrecs = streamptr->tsteps[1].nrecs;

      streamptr->nrecs += nvrecs;

      destTstep->records    = (record_t *) Malloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nvrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = CDI_UNDEFID;

      memcpy(destTstep->records, sourceTstep->records, (size_t)nrecs*sizeof(record_t));

      destTstep->recIDs     = (int *) Malloc((size_t)nvrecs * sizeof(int));

      memcpy(destTstep->recIDs, streamptr->tsteps[1].recIDs, (size_t)nvrecs*sizeof(int));
    }
}

static
int cdf_time_dimid(int fileID, int ndims, int nvars)
{
  char dimname[80];
  for ( int dimid = 0; dimid < ndims; ++dimid )
    {
      dimname[0] = 0;
      cdf_inq_dimname(fileID, dimid, dimname);
      if ( str_is_equal(dimname, "time") || str_is_equal(dimname, "Time") ) return dimid;
    }

  for ( int varid = 0; varid < nvars; ++varid )
    {
      nc_type xtype;
      int nvdims, nvatts, dimids[9];
      cdf_inq_var(fileID, varid, NULL, &xtype, &nvdims, dimids, &nvatts);
      if ( nvdims == 1 )
        {
          char sbuf[CDI_MAX_NAME];
          for ( int iatt = 0; iatt < nvatts; ++iatt )
            {
              sbuf[0] = 0;
              cdf_inq_attname(fileID, varid, iatt, sbuf);
              if ( strncmp(sbuf, "units", 5) == 0 )
                {
                  cdfGetAttText(fileID, varid, "units", sizeof(sbuf), sbuf);
                  str_tolower(sbuf);

                  if ( is_time_units(sbuf) ) return dimids[0];
                }
            }
        }
    }

  return CDI_UNDEFID;
}

static
void init_ncdims(long ndims, ncdim_t *ncdims)
{
  for ( long ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncdims[ncdimid].ncvarid      = CDI_UNDEFID;
      ncdims[ncdimid].dimtype      = CDI_UNDEFID;
      ncdims[ncdimid].len          = 0;
      ncdims[ncdimid].name[0]      = 0;
    }
}

static
void init_ncvars(long nvars, ncvar_t *ncvars)
{
  for ( long ncvarid = 0; ncvarid < nvars; ++ncvarid )
    {
      ncvars[ncvarid].ncid            = CDI_UNDEFID;
      ncvars[ncvarid].isvar           = CDI_UNDEFID;
      ncvars[ncvarid].ignore          = false;
      ncvars[ncvarid].isx             = false;
      ncvars[ncvarid].isy             = false;
      ncvars[ncvarid].isc             = false;
      ncvars[ncvarid].islon           = false;
      ncvars[ncvarid].islat           = false;
      ncvars[ncvarid].islev           = false;
      ncvars[ncvarid].istime          = false;
      ncvars[ncvarid].warn            = false;
      ncvars[ncvarid].calendar        = false;
      ncvars[ncvarid].climatology     = false;
      ncvars[ncvarid].lformulaterms   = false;
      ncvars[ncvarid].timetype        = TIME_CONSTANT;
      ncvars[ncvarid].param           = CDI_UNDEFID;
      ncvars[ncvarid].code            = CDI_UNDEFID;
      ncvars[ncvarid].tabnum          = 0;
      ncvars[ncvarid].bounds          = CDI_UNDEFID;
      ncvars[ncvarid].gridID          = CDI_UNDEFID;
      ncvars[ncvarid].zaxisID         = CDI_UNDEFID;
      ncvars[ncvarid].gridtype        = CDI_UNDEFID;
      ncvars[ncvarid].zaxistype       = CDI_UNDEFID;
      ncvars[ncvarid].xdim            = CDI_UNDEFID;
      ncvars[ncvarid].ydim            = CDI_UNDEFID;
      ncvars[ncvarid].zdim            = CDI_UNDEFID;
      ncvars[ncvarid].xvarid          = CDI_UNDEFID;
      ncvars[ncvarid].yvarid          = CDI_UNDEFID;
      ncvars[ncvarid].zvarid          = CDI_UNDEFID;
      ncvars[ncvarid].tvarid          = CDI_UNDEFID;
      ncvars[ncvarid].psvarid         = CDI_UNDEFID;
      ncvars[ncvarid].p0varid         = CDI_UNDEFID;
      ncvars[ncvarid].ncoordvars      = 0;
      for ( int i = 0; i < MAX_COORDVARS; ++i )
        {
          ncvars[ncvarid].coordvarids[i]  = CDI_UNDEFID;
          ncvars[ncvarid].cvarids[i]      = CDI_UNDEFID;
        }
      ncvars[ncvarid].nauxvars      = 0;
      for ( int i = 0; i < MAX_AUXVARS; ++i )
        ncvars[ncvarid].auxvarids[i]  = CDI_UNDEFID;
      ncvars[ncvarid].cellarea        = CDI_UNDEFID;
      ncvars[ncvarid].tableID         = CDI_UNDEFID;
      ncvars[ncvarid].xtype           = 0;
      ncvars[ncvarid].ndims           = 0;
      ncvars[ncvarid].gmapid          = CDI_UNDEFID;
      ncvars[ncvarid].vctsize         = 0;
      ncvars[ncvarid].vct             = NULL;
      ncvars[ncvarid].truncation      = 0;
      ncvars[ncvarid].position        = 0;
      ncvars[ncvarid].positive        = 0;
      ncvars[ncvarid].chunked         = 0;
      ncvars[ncvarid].chunktype       = CDI_UNDEFID;
      ncvars[ncvarid].defmissval      = false;
      ncvars[ncvarid].deffillval      = false;
      ncvars[ncvarid].missval         = 0;
      ncvars[ncvarid].fillval         = 0;
      ncvars[ncvarid].addoffset       = 0;
      ncvars[ncvarid].scalefactor     = 1;
      ncvars[ncvarid].natts           = 0;
      ncvars[ncvarid].atts            = NULL;
      ncvars[ncvarid].deflate         = 0;
      ncvars[ncvarid].lunsigned       = false;
      ncvars[ncvarid].lvalidrange     = false;
      ncvars[ncvarid].validrange[0]   = VALIDMISS;
      ncvars[ncvarid].validrange[1]   = VALIDMISS;
      ncvars[ncvarid].ensdata         = NULL;
      memset(ncvars[ncvarid].name, 0, CDI_MAX_NAME);
      memset(ncvars[ncvarid].longname, 0, CDI_MAX_NAME);
      memset(ncvars[ncvarid].stdname, 0, CDI_MAX_NAME);
      memset(ncvars[ncvarid].units, 0, CDI_MAX_NAME);
      memset(ncvars[ncvarid].extra, 0, CDI_MAX_NAME);
    }
}

static
void cdf_set_var(ncvar_t *ncvars, int ncvarid, short isvar)
{
  if ( ncvars[ncvarid].isvar != CDI_UNDEFID &&
       ncvars[ncvarid].isvar != isvar   &&
       ncvars[ncvarid].warn  == false )
    {
      if ( ! ncvars[ncvarid].ignore )
        Warning("Inconsistent variable definition for %s!", ncvars[ncvarid].name);

      ncvars[ncvarid].warn = true;
      isvar = FALSE;
    }

  ncvars[ncvarid].isvar = isvar;
}

static
void cdf_set_dim(ncvar_t *ncvars, int ncvarid, int dimid, int dimtype)
{
  if ( ncvars[ncvarid].dimtype[dimid] != CDI_UNDEFID &&
       ncvars[ncvarid].dimtype[dimid] != dimtype )
    {
      Warning("Inconsistent dimension definition for %s! dimid = %d;  type = %d;  newtype = %d",
              ncvars[ncvarid].name, dimid, ncvars[ncvarid].dimtype[dimid], dimtype);
    }

  ncvars[ncvarid].dimtype[dimid] = dimtype;
}

static
void scan_hybrid_formulaterms(int ncid, int ncfvarid, int *avarid, int *bvarid, int *psvarid, int *p0varid)
{
  *avarid  = -1;
  *bvarid  = -1;
  *psvarid = -1;
  *p0varid = -1;

  char attstring[1024];
  cdfGetAttText(ncid, ncfvarid, "formula_terms", sizeof(attstring), attstring);
  char *pstring = attstring;

  bool lstop = false;
  for ( int i = 0; i < 4; i++ )
    {
      while ( isspace((int) *pstring) ) pstring++;
      if ( *pstring == 0 ) break;
      char *tagname = pstring;
      while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
      if ( *pstring == 0 ) lstop = true;
      *pstring++ = 0;

      while ( isspace((int) *pstring) ) pstring++;
      if ( *pstring == 0 ) break;
      char *varname = pstring;
      while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
      if ( *pstring == 0 ) lstop = true;
      *pstring++ = 0;

      int dimvarid;
      int status_nc = nc_inq_varid(ncid, varname, &dimvarid);
      if ( status_nc == NC_NOERR )
        {
          if      ( strcmp(tagname, "ap:") == 0 ) *avarid  = dimvarid;
          else if ( strcmp(tagname, "a:")  == 0 ) *avarid  = dimvarid;
          else if ( strcmp(tagname, "b:")  == 0 ) *bvarid  = dimvarid;
          else if ( strcmp(tagname, "ps:") == 0 ) *psvarid = dimvarid;
          else if ( strcmp(tagname, "p0:") == 0 ) *p0varid = dimvarid;
        }
      else if ( strcmp(tagname, "ps:") != 0 )
        {
          Warning("%s - %s", nc_strerror(status_nc), varname);
        }

      if ( lstop ) break;
    }
}

static
bool isHybridSigmaPressureCoordinate(int ncid, int ncvarid, ncvar_t *ncvars, const ncdim_t *ncdims)
{
  bool status = false;
  ncvar_t *ncvar = &ncvars[ncvarid];

  if ( strcmp(ncvar->stdname, "atmosphere_hybrid_sigma_pressure_coordinate") == 0 )
    {
      cdiConvention = CDI_CONVENTION_CF;

      status = true;
      ncvar->zaxistype = ZAXIS_HYBRID;
      //int ndims = ncvar->ndims;
      int dimid = ncvar->dimids[0];
      size_t dimlen = ncdims[dimid].len;
      int avarid1 = -1, bvarid1 = -1, psvarid1 = -1, p0varid1 = -1;
      int ncfvarid = ncvarid;
      if ( ncvars[ncfvarid].lformulaterms )
        scan_hybrid_formulaterms(ncid, ncfvarid, &avarid1, &bvarid1, &psvarid1, &p0varid1);
      // printf("avarid1, bvarid1, psvarid1, p0varid1 %d %d %d %d\n", avarid1, bvarid1, psvarid1, p0varid1);
      if ( avarid1  != -1 ) ncvars[avarid1].isvar = FALSE;
      if ( bvarid1  != -1 ) ncvars[bvarid1].isvar = FALSE;
      if ( psvarid1 != -1 ) ncvar->psvarid = psvarid1;
      if ( p0varid1 != -1 ) ncvar->p0varid = p0varid1;

      if ( ncvar->bounds != CDI_UNDEFID && ncvars[ncvar->bounds].lformulaterms )
        {
          ncfvarid = ncvar->bounds;
          int avarid2 = -1, bvarid2 = -1, psvarid2 = -1, p0varid2 = -1;
          if ( ncvars[ncfvarid].lformulaterms )
            scan_hybrid_formulaterms(ncid, ncfvarid, &avarid2, &bvarid2, &psvarid2, &p0varid2);
          // printf("avarid2, bvarid2, psvarid2, p0varid2 %d %d %d %d\n", avarid2, bvarid2, psvarid2, p0varid2);
          if ( avarid2 != -1 && bvarid2 != -1 )
            {
              ncvars[avarid2].isvar = FALSE;
              ncvars[bvarid2].isvar = FALSE;

              int ndims2 = ncvars[avarid2].ndims;
              int dimid2 = ncvars[avarid2].dimids[0];
              size_t dimlen2 = ncdims[dimid2].len;

              if ( (ndims2 == 2 && dimid == ncvars[avarid2].dimids[0] ) ||
                   (ndims2 == 1 && dimlen == dimlen2-1 ) )
                {
                  double px = 1;
                  if ( p0varid1 != -1 && p0varid1 == p0varid2 )
                    cdf_get_var_double(ncid, p0varid2, &px);

                  double abuf[dimlen*2], bbuf[dimlen*2];
                  cdf_get_var_double(ncid, avarid2, abuf);
                  cdf_get_var_double(ncid, bvarid2, bbuf);

                  size_t vctsize = (dimlen+1)*2;
                  double *vct = (double *) Malloc(vctsize*sizeof(double));
                  if ( ndims2 == 2 )
                    {
                      for ( size_t i = 0; i < dimlen; ++i )
                        {
                          vct[i] = abuf[i*2];
                          vct[i+dimlen+1] = bbuf[i*2];
                        }
                      vct[dimlen]     = abuf[dimlen*2-1];
                      vct[dimlen*2+1] = bbuf[dimlen*2-1];
                    }
                  else
                    {
                       for ( size_t i = 0; i < dimlen2; ++i )
                        {
                          vct[i] = abuf[i];
                          vct[i+dimlen+1] = bbuf[i];
                        }
                    }

                  if ( p0varid1 != -1 && IS_NOT_EQUAL(px, 1) )
                    for ( size_t i = 0; i < dimlen+1; ++i ) vct[i] *= px;

                  ncvar->vct = vct;
                  ncvar->vctsize = vctsize;
                }
            }
        }
    }

  return status;
}

static
void cdf_set_cdi_attr(int ncid, int ncvarid, int attnum, int cdiID, int varID)
{
  nc_type atttype;
  size_t attlen;
  char attname[CDI_MAX_NAME];

  cdf_inq_attname(ncid, ncvarid, attnum, attname);
  cdf_inq_attlen(ncid, ncvarid, attname, &attlen);
  cdf_inq_atttype(ncid, ncvarid, attname, &atttype);
  if ( xtypeIsInt(atttype) )
    {
      int attint[attlen];
      cdfGetAttInt(ncid, ncvarid, attname, attlen, attint);
      int datatype = (atttype == NC_SHORT)  ? CDI_DATATYPE_INT16 :
                     (atttype == NC_BYTE)   ? CDI_DATATYPE_INT8 :
#if  defined  (HAVE_NETCDF4)
                     (atttype == NC_UBYTE)  ? CDI_DATATYPE_UINT8 :
                     (atttype == NC_USHORT) ? CDI_DATATYPE_UINT16 :
                     (atttype == NC_UINT)   ? CDI_DATATYPE_UINT32 :
#endif
                     CDI_DATATYPE_INT32;
      cdiDefAttInt(cdiID, varID, attname, datatype, (int)attlen, attint);
    }
  else if ( xtypeIsFloat(atttype) )
    {
      double attflt[attlen];
      cdfGetAttDouble(ncid, ncvarid, attname, attlen, attflt);
      int datatype = (atttype == NC_FLOAT) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
      cdiDefAttFlt(cdiID, varID, attname, datatype, (int)attlen, attflt);
    }
  else if ( xtypeIsText(atttype) )
    {
      char attstring[8192];
      cdfGetAttText(ncid, ncvarid, attname, sizeof(attstring), attstring);
      cdiDefAttTxt(cdiID, varID, attname, (int)attlen, attstring);
    }
}

static
void cdf_print_vars(const ncvar_t *ncvars, int nvars, const char *oname)
{
  char axis[7];
  static const char iaxis[] = {'t', 'z', 'y', 'x'};

  fprintf(stderr, "%s:\n", oname);

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      int ndim = 0;
      if ( ncvars[ncvarid].isvar )
        {
          axis[ndim++] = 'v';
          axis[ndim++] = ':';
          for ( int i = 0; i < ncvars[ncvarid].ndims; i++ )
            {/*
              if      ( ncvars[ncvarid].tvarid != -1 ) axis[ndim++] = iaxis[0];
              else if ( ncvars[ncvarid].zvarid != -1 ) axis[ndim++] = iaxis[1];
              else if ( ncvars[ncvarid].yvarid != -1 ) axis[ndim++] = iaxis[2];
              else if ( ncvars[ncvarid].xvarid != -1 ) axis[ndim++] = iaxis[3];
              else
             */
              if      ( ncvars[ncvarid].dimtype[i] == T_AXIS ) axis[ndim++] = iaxis[0];
              else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) axis[ndim++] = iaxis[1];
              else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) axis[ndim++] = iaxis[2];
              else if ( ncvars[ncvarid].dimtype[i] == X_AXIS ) axis[ndim++] = iaxis[3];
              else                                             axis[ndim++] = '?';
            }
        }
      else
        {
          axis[ndim++] = 'c';
          axis[ndim++] = ':';
          if      ( ncvars[ncvarid].istime ) axis[ndim++] = iaxis[0];
          else if ( ncvars[ncvarid].islev  ) axis[ndim++] = iaxis[1];
          else if ( ncvars[ncvarid].islat  ) axis[ndim++] = iaxis[2];
          else if ( ncvars[ncvarid].isy    ) axis[ndim++] = iaxis[2];
          else if ( ncvars[ncvarid].islon  ) axis[ndim++] = iaxis[3];
          else if ( ncvars[ncvarid].isx    ) axis[ndim++] = iaxis[3];
          else                               axis[ndim++] = '?';
        }

      axis[ndim++] = 0;

      fprintf(stderr, "%3d %3d  %-6s %s\n", ncvarid, ndim-3, axis, ncvars[ncvarid].name);
    }
}

static
void cdf_scan_attr_axis(ncvar_t *ncvars, ncdim_t *ncdims, int ncvarid, const char *attstring, size_t attlen,
                        int nvdims, int *dimidsp, const char *name)
{
  int i;
  for ( i = 0; i < (int)attlen; ++i )
    {
      if ( attstring[i] != '-' && attstring[i] != 't' && attstring[i] != 'z' &&
           attstring[i] != 'y' && attstring[i] != 'x' )
        {
          Warning("Unexpected character in axis attribute for %s, ignored!", name);
          break;
        }
    }

  if ( i == (int) attlen && (int) attlen == nvdims )
    {
      while ( attlen-- )
        {
          if ( (int) attstring[attlen] == 't' )
            {
              if ( attlen != 0 ) Warning("axis attribute 't' not on first position");
              cdf_set_dim(ncvars, ncvarid, (int)attlen, T_AXIS);
            }
          else if ( (int) attstring[attlen] == 'z' )
            {
              ncvars[ncvarid].zdim = dimidsp[attlen];
              cdf_set_dim(ncvars, ncvarid, (int)attlen, Z_AXIS);

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdf_set_var(ncvars, ncvarid, FALSE);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
                }
            }
          else if ( (int) attstring[attlen] == 'y' )
            {
              ncvars[ncvarid].ydim = dimidsp[attlen];
              cdf_set_dim(ncvars, ncvarid, (int)attlen, Y_AXIS);

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdf_set_var(ncvars, ncvarid, FALSE);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Y_AXIS;
                }
            }
          else if ( (int) attstring[attlen] == 'x' )
            {
              ncvars[ncvarid].xdim = dimidsp[attlen];
              cdf_set_dim(ncvars, ncvarid, (int)attlen, X_AXIS);

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdf_set_var(ncvars, ncvarid, FALSE);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = X_AXIS;
                }
            }
        }
    }
}

static
int cdf_get_cell_varid(char *attstring, int ncid)
{
  int nc_cell_id = CDI_UNDEFID;

  char *pstring = attstring;
  while ( isspace((int) *pstring) ) pstring++;
  char *cell_measures = pstring;
  while ( isalnum((int) *pstring) ) pstring++;
  *pstring++ = 0;
  while ( isspace((int) *pstring) ) pstring++;
  char *cell_var = pstring;
  while ( ! isspace((int) *pstring) && *pstring != 0 ) pstring++;
  *pstring++ = 0;
  /*
    printf("cell_measures >%s<\n", cell_measures);
    printf("cell_var >%s<\n", cell_var);
  */
  if ( str_is_equal(cell_measures, "area") )
    {
      int nc_var_id;
      int status = nc_inq_varid(ncid, cell_var, &nc_var_id);
      if ( status == NC_NOERR )
        nc_cell_id = nc_var_id;
      /*
      else
        Warning("%s - %s", nc_strerror(status), cell_var);
      */
    }

  return nc_cell_id;
}

static
void cdf_scan_var_attr(int nvars, ncvar_t *ncvars, ncdim_t *ncdims, int timedimid, int modelID, int format)
{
  int ncdimid;
  int nvdims, nvatts;
  int iatt;
  nc_type xtype, atttype;
  size_t attlen;
  char name[CDI_MAX_NAME];
  char attname[CDI_MAX_NAME];
  char attstring[8192];

  int nchecked_vars = 0;
  enum { max_check_vars = 9 };
  char *checked_vars[max_check_vars];
  for ( int i = 0; i < max_check_vars; ++i ) checked_vars[i] = NULL;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      int ncid    = ncvars[ncvarid].ncid;
      int *dimidsp = ncvars[ncvarid].dimids;

      cdf_inq_var(ncid, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts);
      strcpy(ncvars[ncvarid].name, name);

      for ( ncdimid = 0; ncdimid < nvdims; ncdimid++ )
        ncvars[ncvarid].dimtype[ncdimid] = -1;

      ncvars[ncvarid].xtype = xtype;
      ncvars[ncvarid].ndims = nvdims;

#if  defined  (HAVE_NETCDF4)
      if ( format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4 )
        {
          int shuffle, deflate, deflate_level;
          size_t chunks[nvdims];
          int storage_in;
          nc_inq_var_deflate(ncid, ncvarid, &shuffle, &deflate, &deflate_level);
          if ( deflate > 0 ) ncvars[ncvarid].deflate = 1;
          /*
          size_t cache_size, nelems;
          float preemption;
          nc_get_chunk_cache(&cache_size, &nelems, &preemption);
          printf("cache_size %lu nelems %lu preemption %g\n", cache_size, nelems, preemption);
          nc_get_var_chunk_cache(ncid, ncvarid, &cache_size, &nelems, &preemption);
          printf("varid %d cache_size %lu nelems %lu preemption %g\n", ncvarid, cache_size, nelems, preemption);
          */
          if ( nc_inq_var_chunking(ncid, ncvarid, &storage_in, chunks) == NC_NOERR )
            {
              if ( storage_in == NC_CHUNKED )
                {
                  ncvars[ncvarid].chunked = 1;
                  for ( int i = 0; i < nvdims; ++i ) ncvars[ncvarid].chunks[i] = (int)chunks[i];
                  if ( CDI_Debug )
                    {
                      fprintf(stderr, "%s: chunking %d %d %d  chunks ", name, storage_in, NC_CONTIGUOUS, NC_CHUNKED);
                      for ( int i = 0; i < nvdims; ++i ) fprintf(stderr, "%ld ", chunks[i]);
                      fprintf(stderr, "\n");
                    }
                  {
                    char *buf = ncvars[ncvarid].extra;
                    size_t pos = strlen(buf);
                    static const char prefix[] = "chunks=";
                    memcpy(buf + pos, prefix, sizeof (prefix));
                    pos += sizeof (prefix) - 1;
                    for ( int i = nvdims-1; i >= 0; --i )
                      {
                        pos += (size_t)(sprintf(buf + pos, "%zu%s", chunks[i],
                                                i > 0 ? "x" : ""));
                      }
                    buf[pos] = ' '; buf[pos + 1] = 0;
                  }
                }
            }
        }
#endif

      if ( nvdims > 0 )
        {
          if ( timedimid == dimidsp[0] )
            {
              ncvars[ncvarid].timetype = TIME_VARYING;
              cdf_set_dim(ncvars, ncvarid, 0, T_AXIS);
            }
          else
            {
              for ( ncdimid = 1; ncdimid < nvdims; ncdimid++ )
                {
                  if ( timedimid == dimidsp[ncdimid] )
                    {
                      Warning("Time must be the first dimension! Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
                      ncvars[ncvarid].isvar = FALSE;
                    }
                }
            }
        }

      for ( iatt = 0; iatt < nvatts; iatt++ )
        {
          int nc_cell_id = CDI_UNDEFID;

          cdf_inq_attname(ncid, ncvarid, iatt, attname);
          cdf_inq_atttype(ncid, ncvarid, attname, &atttype);
          cdf_inq_attlen(ncid, ncvarid, attname, &attlen);

          size_t attstringsize = sizeof(attstring);
          bool isText = xtypeIsText(atttype);
          bool isNumber = xtypeIsFloat(atttype) || xtypeIsInt(atttype);
          if ( isText )
            {
              cdfGetAttText(ncid, ncvarid, attname, sizeof(attstring), attstring);
              attstringsize = strlen(attstring) + 1;
              if ( attstringsize > CDI_MAX_NAME ) attstringsize = CDI_MAX_NAME;
            }

          if ( isText && strcmp(attname, "long_name") == 0 )
            {
              memcpy(ncvars[ncvarid].longname, attstring, attstringsize);
            }
          else if ( isText && strcmp(attname, "standard_name") == 0 )
            {
              memcpy(ncvars[ncvarid].stdname, attstring, attstringsize);
            }
          else if ( isText && strcmp(attname, "units") == 0 )
            {
              memcpy(ncvars[ncvarid].units, attstring, attstringsize);
            }
          else if ( strcmp(attname, "calendar") == 0 )
            {
              ncvars[ncvarid].calendar = true;
            }
          else if ( isText && strcmp(attname, "param") == 0 )
            {
	      int pnum = 0, pcat = 255, pdis = 255;
	      sscanf(attstring, "%d.%d.%d", &pnum, &pcat, &pdis);
	      ncvars[ncvarid].param = cdiEncodeParam(pnum, pcat, pdis);
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isNumber && strcmp(attname, "code") == 0 )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].code);
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isNumber && strcmp(attname, "table") == 0 )
            {
              int tablenum;
              cdfGetAttInt(ncid, ncvarid, attname, 1, &tablenum);
              if ( tablenum > 0 )
                {
                  ncvars[ncvarid].tabnum = tablenum;
                  ncvars[ncvarid].tableID = tableInq(modelID, tablenum, NULL);
                  if ( ncvars[ncvarid].tableID == CDI_UNDEFID )
                    ncvars[ncvarid].tableID = tableDef(modelID, tablenum, NULL);
                }
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isText && strcmp(attname, "trunc_type") == 0 )
            {
              if ( str_is_equal(attstring, "Triangular") )
                ncvars[ncvarid].gridtype = GRID_SPECTRAL;
            }
          else if ( isText && (strcmp(attname, "grid_type") == 0 || strcmp(attname, "CDI_grid_type") == 0) )
            {
              str_tolower(attstring);
              set_gridtype(attstring, &ncvars[ncvarid].gridtype);
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isText && strcmp(attname, "level_type") == 0 )
            {
              str_tolower(attstring);
              set_zaxistype(attstring, &ncvars[ncvarid].zaxistype);
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isNumber && strcmp(attname, "trunc_count") == 0 )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( isNumber && strcmp(attname, "truncation") == 0 )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( isNumber && strcmp(attname, "number_of_grid_in_reference") == 0 )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].position);
            }
          else if ( isNumber && strcmp(attname, "add_offset") == 0 )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].addoffset);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].addoffset != 0 )
		Warning("attribute add_offset not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdf_set_var(ncvars, ncvarid, TRUE); */
            }
          else if ( isNumber && strcmp(attname, "scale_factor") == 0 )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].scalefactor);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].scalefactor != 1 )
		Warning("attribute scale_factor not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdf_set_var(ncvars, ncvarid, TRUE); */
            }
          else if ( isText && strcmp(attname, "climatology") == 0 )
            {
              int ncboundsid;
              int status = nc_inq_varid(ncid, attstring, &ncboundsid);
              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].climatology = true;
                  ncvars[ncvarid].bounds = ncboundsid;
                  cdf_set_var(ncvars, ncvars[ncvarid].bounds, FALSE);
                  cdf_set_var(ncvars, ncvarid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);
            }
          else if ( isText && strcmp(attname, "bounds") == 0 )
            {
              int ncboundsid;
              int status = nc_inq_varid(ncid, attstring, &ncboundsid);
              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].bounds = ncboundsid;
                  cdf_set_var(ncvars, ncvars[ncvarid].bounds, FALSE);
                  cdf_set_var(ncvars, ncvarid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);
            }
          else if ( isText &&  strcmp(attname, "formula_terms") == 0 )
            {
              ncvars[ncvarid].lformulaterms = true;
            }
          else if ( isText && strcmp(attname, "cell_measures") == 0 && (nc_cell_id=cdf_get_cell_varid(attstring, ncid)) != CDI_UNDEFID )
            {
              ncvars[ncvarid].cellarea = nc_cell_id;
              ncvars[nc_cell_id].isvar = FALSE;
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          /*
          else if ( strcmp(attname, "coordinates") == 0 )
            {
              char *pstring, *xvarname = NULL, *yvarname = NULL;
              pstring = attstring;

              while ( isspace((int) *pstring) ) pstring++;
              xvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;
              while ( isspace((int) *pstring) ) pstring++;
              yvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;

              cdf_inq_varid(ncid, xvarname, &ncvars[ncvarid].xvarid);
              cdf_inq_varid(ncid, yvarname, &ncvars[ncvarid].yvarid);

              cdf_set_var(ncvars, ncvars[ncvarid].xvarid, FALSE);
              cdf_set_var(ncvars, ncvars[ncvarid].yvarid, FALSE);
              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          */
          else if ( isText && (strcmp(attname, "associate")  == 0 || strcmp(attname, "coordinates") == 0) )
            {
              bool lstop = false;
              char *pstring = attstring;

              for ( int i = 0; i < MAX_COORDVARS; i++ )
                {
                  while ( isspace((int) *pstring) ) pstring++;
                  if ( *pstring == 0 ) break;
                  char *varname = pstring;
                  while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
                  if ( *pstring == 0 ) lstop = true;
                  *pstring++ = 0;

                  int dimvarid;
                  int status = nc_inq_varid(ncid, varname, &dimvarid);
                  if ( status == NC_NOERR )
                    {
                      cdf_set_var(ncvars, dimvarid, FALSE);
                      if ( !cdiIgnoreAttCoordinates )
                        {
                          ncvars[ncvarid].coordvarids[i] = dimvarid;
                          ncvars[ncvarid].ncoordvars++;
                        }
                    }
                  else
                    {
                      int k;
                      for ( k = 0; k < nchecked_vars; ++k )
                        if ( strcmp(checked_vars[k], varname) == 0 ) break;

                      if ( k == nchecked_vars )
                        {
                          if ( nchecked_vars < max_check_vars ) checked_vars[nchecked_vars++] = strdup(varname);
                          Warning("%s - %s", nc_strerror(status), varname);
                        }
                    }

                  if ( lstop ) break;
                }

              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isText && strcmp(attname, "auxiliary_variable") == 0 )
            {
              bool lstop = false;
              char *pstring = attstring;

              for ( int i = 0; i < MAX_AUXVARS; i++ )
                {
                  while ( isspace((int) *pstring) ) pstring++;
                  if ( *pstring == 0 ) break;
                  char *varname = pstring;
                  while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
                  if ( *pstring == 0 ) lstop = true;
                  *pstring++ = 0;

                  int dimvarid;
                  int status = nc_inq_varid(ncid, varname, &dimvarid);
                  if ( status == NC_NOERR )
                    {
                      cdf_set_var(ncvars, dimvarid, FALSE);
                      //  if ( !cdiIgnoreAttCoordinates )
                        {
                          ncvars[ncvarid].auxvarids[i] = dimvarid;
                          ncvars[ncvarid].nauxvars++;
                        }
                    }
                  else
                    Warning("%s - %s", nc_strerror(status), varname);

                  if ( lstop ) break;
                }

              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isText && strcmp(attname, "grid_mapping") == 0 )
            {
              int nc_gmap_id;
              int status = nc_inq_varid(ncid, attstring, &nc_gmap_id);
              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].gmapid = nc_gmap_id;
                  cdf_set_var(ncvars, ncvars[ncvarid].gmapid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);

              cdf_set_var(ncvars, ncvarid, TRUE);
            }
          else if ( isText && strcmp(attname, "positive") == 0 )
            {
              str_tolower(attstring);

              if      ( str_is_equal(attstring, "down") ) ncvars[ncvarid].positive = POSITIVE_DOWN;
              else if ( str_is_equal(attstring, "up")   ) ncvars[ncvarid].positive = POSITIVE_UP;

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdf_set_var(ncvars, ncvarid, FALSE);
                  cdf_set_dim(ncvars, ncvarid, 0, Z_AXIS);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
                }
            }
          else if ( isNumber && strcmp(attname, "_FillValue") == 0 )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].fillval);
	      ncvars[ncvarid].deffillval = true;
	      /* cdf_set_var(ncvars, ncvarid, TRUE); */
            }
          else if ( isNumber && strcmp(attname, "missing_value") == 0 )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].missval);
	      ncvars[ncvarid].defmissval = true;
	      /* cdf_set_var(ncvars, ncvarid, TRUE); */
            }
          else if ( isNumber && strcmp(attname, "valid_range") == 0 && attlen == 2 )
            {
              if ( ncvars[ncvarid].lvalidrange == false )
                {
                  bool lignore = xtypeIsFloat(atttype) != xtypeIsFloat(xtype);
                  if ( !cdiIgnoreValidRange && lignore == false )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 2, ncvars[ncvarid].validrange);
                      ncvars[ncvarid].lvalidrange = true;
                      if ( ((int)ncvars[ncvarid].validrange[0]) == 0 && ((int)ncvars[ncvarid].validrange[1]) == 255 )
                        ncvars[ncvarid].lunsigned = true;
                      /* cdf_set_var(ncvars, ncvarid, TRUE); */
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_range, ignored!", name);
                    }
                }
            }
          else if ( isNumber && strcmp(attname, "valid_min") == 0 && attlen == 1 )
            {
              if ( ncvars[ncvarid].lvalidrange == false )
                {
                  bool lignore = xtypeIsFloat(atttype) != xtypeIsFloat(xtype);
                  if ( !cdiIgnoreValidRange && lignore == false )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[0]);
                      ncvars[ncvarid].lvalidrange = true;
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_min, ignored!", name);
                    }
                }
            }
          else if ( isNumber && strcmp(attname, "valid_max") == 0 && attlen == 1 )
            {
              if ( ncvars[ncvarid].lvalidrange == false )
                {
                  bool lignore = xtypeIsFloat(atttype) != xtypeIsFloat(xtype);
                  if ( !cdiIgnoreValidRange && lignore == false )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[1]);
                      ncvars[ncvarid].lvalidrange = true;
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_max, ignored!", name);
                    }
                }
            }
          else if ( isText && strcmp(attname, "_Unsigned") == 0 )
            {
              str_tolower(attstring);

              if ( str_is_equal(attstring, "true") )
                {
                  ncvars[ncvarid].lunsigned = true;
                  /*
                  ncvars[ncvarid].lvalidrange = true;
                  ncvars[ncvarid].validrange[0] = 0;
                  ncvars[ncvarid].validrange[1] = 255;
                  */
                }
	      /* cdf_set_var(ncvars, ncvarid, TRUE); */
            }
          else if ( isText && strcmp(attname, "cdi") == 0 )
            {
	      str_tolower(attstring);

	      if ( str_is_equal(attstring, "ignore") )
		{
		  ncvars[ncvarid].ignore = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		}
            }
          else if ( isText && strcmp(attname, "axis") == 0 )
            {
	      attlen = strlen(attstring);

	      if ( (int) attlen > nvdims && nvdims > 0 && attlen > 1 )
		{
		    Warning("Unexpected axis attribute length for %s, ignored!", name);
		}
              else if ( nvdims == 0 && attlen == 1 )
                {
                  if ( attstring[0] == 'z' || attstring[0] == 'Z' )
                    {
                      cdf_set_var(ncvars, ncvarid, FALSE);
                      ncvars[ncvarid].islev = true;
                    }
                }
	      else
		{
		  str_tolower(attstring);
                  cdf_scan_attr_axis(ncvars, ncdims, ncvarid, attstring, attlen, nvdims, dimidsp, name);
		}
	    }
	  else if ( isNumber &&
                    (strcmp(attname, "realization") == 0       ||
                     strcmp(attname, "ensemble_members") == 0  ||
                     strcmp(attname, "forecast_init_type") == 0) )
	    {
	      int temp;

	      if( ncvars[ncvarid].ensdata == NULL )
		ncvars[ncvarid].ensdata = (ensinfo_t *) Malloc( sizeof( ensinfo_t ) );

	      cdfGetAttInt(ncid, ncvarid, attname, 1, &temp);

	      if( strcmp(attname, "realization") == 0 )
		ncvars[ncvarid].ensdata->ens_index = temp;
	      else if( strcmp(attname, "ensemble_members") == 0 )
		ncvars[ncvarid].ensdata->ens_count = temp;
	      else if( strcmp(attname, "forecast_init_type") == 0 )
		ncvars[ncvarid].ensdata->forecast_init_type = temp;

	      cdf_set_var(ncvars, ncvarid, TRUE);
	    }
	  else
	    {
	      if ( ncvars[ncvarid].natts == 0 )
		ncvars[ncvarid].atts = (int*) Malloc((size_t)nvatts*sizeof(int));

	      ncvars[ncvarid].atts[ncvars[ncvarid].natts++] = iatt;
	      /*
	      int attrint;
	      double attrflt;
	      nc_type atttype;
	      cdf_inq_attlen(ncid, ncvarid, attname, &attlen);
	      cdf_inq_atttype(ncid, ncvarid, attname, &atttype);
	      if ( attlen == 1 && (atttype == NC_INT || atttype == NC_SHORT) )
		{
		  cdfGetAttInt(ncid, ncvarid, attname, 1, &attrint);
		  printf("int: %s.%s = %d\n", ncvars[ncvarid].name, attname, attrint);
		}
	      else if ( attlen == 1 && (atttype == NC_FLOAT || atttype == NC_DOUBLE) )
		{
		  cdfGetAttDouble(ncid, ncvarid, attname, 1, &attrflt);
		  printf("flt: %s.%s = %g\n", ncvars[ncvarid].name, attname, attrflt);
		}
	      else if ( atttype == NC_CHAR )
		{
		  attstring[attlen] = 0;
		  printf("txt: %s.%s = %s\n", ncvars[ncvarid].name, attname, attstring);
		}
	      else
		printf("att: %s.%s = unknown\n", ncvars[ncvarid].name, attname);
              */
	    }
	}
    }

  for ( int i = 0; i < max_check_vars; ++i ) if ( checked_vars[i] ) Free(checked_vars[i]);
}

static
void cdf_set_dimtype(int nvars, ncvar_t *ncvars, ncdim_t *ncdims)
{
  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE )
	{
	  int ndims = ncvars[ncvarid].ndims;
	  for ( int i = 0; i < ndims; i++ )
	    {
	      int ncdimid = ncvars[ncvarid].dimids[i];
              int dimtype = ncdims[ncdimid].dimtype;
	      if ( dimtype >= X_AXIS && dimtype <= T_AXIS )
                cdf_set_dim(ncvars, ncvarid, i, dimtype);
	    }

	  if ( CDI_Debug )
	    {
	      Message("var %d %s", ncvarid, ncvars[ncvarid].name);
	      for ( int i = 0; i < ndims; i++ )
		printf("  dim%d type=%d  ", i, ncvars[ncvarid].dimtype[i]);
	      printf("\n");
	    }
          }
      }

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE )
	{
	  bool lxdim = false, lydim = false, lzdim = false/* , ltdim = false */;
          int lcdim = 0;
	  int ndims = ncvars[ncvarid].ndims;
	  for ( int i = 0; i < ndims; i++ )
	    {
              int dimtype = ncvars[ncvarid].dimtype[i];
              lxdim = lxdim | (dimtype == X_AXIS);
	      lydim = lydim | (dimtype == Y_AXIS);
	      lzdim = lzdim | (dimtype == Z_AXIS);
              if ( ncvars[ncvarid].cvarids[i] != CDI_UNDEFID ) lcdim++;
	      /* else if ( ncvars[ncvarid].dimtype[i] == T_AXIS ) ltdim = true; */
	    }

          int allcdims = lcdim;

          if ( !lxdim && ncvars[ncvarid].xvarid != CDI_UNDEFID )
            {
              if ( ncvars[ncvars[ncvarid].xvarid].ndims == 0 ) lxdim = true;
            }

          if ( !lydim && ncvars[ncvarid].yvarid != CDI_UNDEFID )
            {
              if ( ncvars[ncvars[ncvarid].yvarid].ndims == 0 ) lydim = true;
            }

          if ( lxdim && (lydim || ncvars[ncvarid].gridtype == GRID_UNSTRUCTURED) )
            for ( int i = ndims-1; i >= 0; i-- )
              {
                if ( ncvars[ncvarid].dimtype[i] == -1 )
                  {
                    if ( !lzdim )
                      {
                        if ( lcdim )
                          {
                            int cdimvar = ncvars[ncvarid].cvarids[allcdims-lcdim];
                            ncvars[ncvarid].zvarid = cdimvar;
                            lcdim--;
		            ncvars[cdimvar].zaxistype = ZAXIS_CHAR;
                          }
                        cdf_set_dim(ncvars, ncvarid, i, Z_AXIS);
                        lzdim = true;
                        int ncdimid = ncvars[ncvarid].dimids[i];
                        if ( ncdims[ncdimid].dimtype == CDI_UNDEFID )
                          ncdims[ncdimid].dimtype = Z_AXIS;
                      }
                  }
              }
	}
    }

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      int ndims = ncvars[ncvarid].ndims;
      for ( int i = 0; i < ndims; i++ )
        {
          if ( ncvars[ncvarid].dimtype[i] == CDI_UNDEFID )
            {
              int ncdimid = ncvars[ncvarid].dimids[i];
              if ( ncdims[ncdimid].dimtype == Z_AXIS )
                {
                  ncvars[ncvarid].islev = true;
                  cdf_set_dim(ncvars, ncvarid, i, Z_AXIS);
                }
            }
        }
    }

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE )
	{
	  bool lxdim = false, lydim = false, lzdim = false/* , ltdim = false */;
          int lcdim = 0;
	  int ndims = ncvars[ncvarid].ndims;
	  for ( int i = 0; i < ndims; i++ )
	    {
	      if      ( ncvars[ncvarid].dimtype[i] == X_AXIS ) lxdim = true;
	      else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) lydim = true;
	      else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) lzdim = true;
              else if ( ncvars[ncvarid].cvarids[i] != CDI_UNDEFID ) lcdim++;
	      /* else if ( ncvars[ncvarid].dimtype[i] == T_AXIS ) ltdim = true; */
	    }

          int allcdims = lcdim;

          if ( !lxdim && ncvars[ncvarid].xvarid != CDI_UNDEFID )
            {
              if ( ncvars[ncvars[ncvarid].xvarid].ndims == 0 ) lxdim = true;
            }

          if ( !lydim && ncvars[ncvarid].yvarid != CDI_UNDEFID )
            {
              if ( ncvars[ncvars[ncvarid].yvarid].ndims == 0 ) lydim = true;
            }

          //   if ( ndims > 1 )
            for ( int i = ndims-1; i >= 0; i-- )
              {
                if ( ncvars[ncvarid].dimtype[i] == -1 )
                  {
                    if ( !lxdim )
                      {
                        if ( lcdim && ncvars[ncvarid].xvarid == CDI_UNDEFID )
                          {
                            int cdimvar = ncvars[ncvarid].cvarids[allcdims-lcdim];
                            ncvars[ncvarid].xvarid = cdimvar;
                            lcdim--;
                          }
                        cdf_set_dim(ncvars, ncvarid, i, X_AXIS);
                        lxdim = true;
                      }
                    else if ( !lydim && ncvars[ncvarid].gridtype != GRID_UNSTRUCTURED )
                      {
                        if ( lcdim && ncvars[ncvarid].yvarid == CDI_UNDEFID )
                          {
                            int cdimvar = ncvars[ncvarid].cvarids[allcdims-lcdim];
                            ncvars[ncvarid].yvarid = cdimvar;
                            lcdim--;
                          }
                        cdf_set_dim(ncvars, ncvarid, i, Y_AXIS);
                        lydim = true;
                      }
                    else if ( !lzdim )
                      {
                        if ( lcdim > 0 )
                          {
                            int cdimvar = ncvars[ncvarid].cvarids[allcdims-lcdim];
                            ncvars[ncvarid].zvarid = cdimvar;
                            lcdim--;
		            ncvars[cdimvar].zaxistype = ZAXIS_CHAR;
                          }
                        cdf_set_dim(ncvars, ncvarid, i, Z_AXIS);
                        lzdim = true;
                      }
                  }
              }
	}
    }
}

/* verify coordinate vars - first scan (dimname == varname) */
static
void verify_coordinate_vars_1(int ncid, int ndims, ncdim_t *ncdims, ncvar_t *ncvars, int timedimid, bool *lhybrid_cf)
{
  for ( int ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      int ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  if ( ncvars[ncvarid].dimids[0] == timedimid )
	    {
              ncvars[ncvarid].istime = true;
	      ncdims[ncdimid].dimtype = T_AXIS;
	      continue;
	    }

          if ( isHybridSigmaPressureCoordinate(ncid, ncvarid, ncvars, ncdims) )
            {
              *lhybrid_cf = true;
              continue;
            }

	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( is_lon_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		}
	      else if ( is_lat_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		}
	      else if ( is_x_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].isx = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		}
	      else if ( is_y_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].isy = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		}
	      else if ( is_pressure_units(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0 )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		}
	      else if ( is_DBL_axis(ncvars[ncvarid].longname) )
                {
                  ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		}
	      else if ( is_height_units(ncvars[ncvarid].units) )
		{
		  if ( is_depth_axis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( is_height_axis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		}
	    }
          else
            {
              if ( (strcmp(ncvars[ncvarid].longname, "generalized_height") == 0 ||
                    strcmp(ncvars[ncvarid].longname, "generalized height") == 0) &&
                   strcmp(ncvars[ncvarid].stdname, "height") == 0 )
                  ncvars[ncvarid].zaxistype = ZAXIS_REFERENCE;
            }

	  if ( !ncvars[ncvarid].islon && ncvars[ncvarid].longname[0] != 0 &&
               !ncvars[ncvarid].islat && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( str_is_equal(ncvars[ncvarid].longname+1, "ongitude") )
		{
		  ncvars[ncvarid].islon = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		  continue;
		}
	      else if ( str_is_equal(ncvars[ncvarid].longname+1, "atitude") )
		{
		  ncvars[ncvarid].islat = true;
		  cdf_set_var(ncvars, ncvarid, FALSE);
		  cdf_set_dim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		  continue;
		}
	    }

	  if ( ncvars[ncvarid].zaxistype != CDI_UNDEFID )
	    {
              ncvars[ncvarid].islev = true;
	      cdf_set_var(ncvars, ncvarid, FALSE);
	      cdf_set_dim(ncvars, ncvarid, 0, Z_AXIS);
	      ncdims[ncdimid].dimtype = Z_AXIS;
	    }
	}
    }
}

/* verify coordinate vars - second scan (all other variables) */
static
void verify_coordinate_vars_2(int nvars, ncvar_t *ncvars)
{
  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == 0 )
	{
	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( is_lon_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = true;
		  continue;
		}
	      else if ( is_lat_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = true;
		  continue;
		}
	      else if ( is_x_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].isx = true;
		  continue;
		}
	      else if ( is_y_axis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].isy = true;
		  continue;
		}
	      else if ( ncvars[ncvarid].zaxistype == CDI_UNDEFID &&
                        (strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0) )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		  continue;
		}
	      else if ( ncvars[ncvarid].zaxistype == CDI_UNDEFID && is_pressure_units(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		  continue;
		}
	      else if ( is_DBL_axis(ncvars[ncvarid].longname) )
		{
                  ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		  continue;
		}
	      else if ( is_height_units(ncvars[ncvarid].units) )
		{
		  if ( is_depth_axis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( is_height_axis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		  continue;
		}
            }
          else if ( strcmp(ncvars[ncvarid].stdname, "region") == 0  ||
                    strcmp(ncvars[ncvarid].stdname, "area_type") == 0 ||
                    cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned) == CDI_DATATYPE_UINT8 )
            {
              ncvars[ncvarid].isc = true;
            }

	  /* not needed anymore for rotated grids */
	  if ( !ncvars[ncvarid].islon && ncvars[ncvarid].longname[0] != 0 &&
               !ncvars[ncvarid].islat && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( str_is_equal(ncvars[ncvarid].longname+1, "ongitude") )
		{
		  ncvars[ncvarid].islon = true;
		  continue;
		}
	      else if ( str_is_equal(ncvars[ncvarid].longname+1, "atitude") )
		{
		  ncvars[ncvarid].islat = true;
		  continue;
		}
	    }
	}
    }
}

static
void grid_set_chunktype(grid_t *grid, ncvar_t *ncvar)
{
  if ( ncvar->chunked )
    {
      int ndims = ncvar->ndims;

      if ( grid->type == GRID_UNSTRUCTURED )
        {
          ncvar->chunktype = ncvar->chunks[ndims-1] == grid->size
            ? CDI_CHUNK_GRID : CDI_CHUNK_AUTO;
        }
      else
        {
          if ( grid->x.size > 1 && grid->y.size > 1 && ndims > 1 &&
               grid->x.size == ncvar->chunks[ndims-1] &&
               grid->y.size == ncvar->chunks[ndims-2] )
            ncvar->chunktype = CDI_CHUNK_GRID;
          else if ( grid->x.size > 1 && grid->x.size == ncvar->chunks[ndims-1] )
            ncvar->chunktype = CDI_CHUNK_LINES;
          else
            ncvar->chunktype = CDI_CHUNK_AUTO;
        }
    }
}

/* define all input grids */
static
void cdf_load_vals(size_t size, int ndims, int varid, ncvar_t *ncvar, double **gridvals, struct xyValGet *valsGet,
                   int ntdims, size_t *start, size_t *count)
{
  if ( CDI_netcdf_lazy_grid_load )
    {
      *valsGet = (struct xyValGet){
        .scalefactor = ncvar->scalefactor,
        .addoffset = ncvar->addoffset,
        .start = { start[0], start[1], start[2] },
        .count = { count[0], count[1], count[2] },
        .size = size,
        .datasetNCId = ncvar->ncid,
        .varNCId = varid,
        .ndims = (short)ndims,
      };
      *gridvals = cdfPendingLoad;
    }
  else
    {
      *gridvals = (double*) Malloc(size*sizeof(double));
      if ( ntdims == 1 )
        cdf_get_vara_double(ncvar->ncid, varid, start, count, *gridvals);
      else
        cdf_get_var_double(ncvar->ncid, varid, *gridvals);
      cdf_scale_add(size, *gridvals, ncvar->addoffset, ncvar->scalefactor);
    }
}

static
void cdf_load_cvals(size_t size, int varid, ncvar_t *ncvar, char ***gridvals, size_t dimlength)
{
  size_t startc[] = {0, 0};
  size_t countc[] = {1, size/dimlength};
  *gridvals = (char **) Malloc(dimlength * sizeof(char *));
  for ( size_t i = 0; i < dimlength; i++ )
    {
      (*gridvals)[i] = (char*) Malloc((size/dimlength) * sizeof(char));
      cdf_get_vara_text(ncvar->ncid, varid, startc, countc, (*gridvals)[i]);
      startc[0] = i+1;
    }
}

static
void cdf_load_bounds(size_t size, ncvar_t *ncvar, double **gridbounds, struct cdfLazyGridIds *cellBoundsGet)
{
  if ( CDI_netcdf_lazy_grid_load )
    {
      cellBoundsGet->datasetNCId = ncvar->ncid;
      cellBoundsGet->varNCId  = ncvar->bounds;
      *gridbounds = cdfPendingLoad;
    }
  else
    {
      *gridbounds = (double*) Malloc(size*sizeof(double));
      cdf_get_var_double(ncvar->ncid, ncvar->bounds, *gridbounds);
    }
}

static
void cdf_load_cellarea(size_t size, ncvar_t *ncvar, double **gridarea, struct cdfLazyGridIds *cellAreaGet)
{
  if ( CDI_netcdf_lazy_grid_load )
    {
      cellAreaGet->datasetNCId = ncvar->ncid;
      cellAreaGet->varNCId = ncvar->cellarea;
      *gridarea = cdfPendingLoad;
    }
  else
    {
      *gridarea = (double*) Malloc(size*sizeof(double));
      cdf_get_var_double(ncvar->ncid, ncvar->cellarea, *gridarea);
    }
}

static
void cdf_copy_axis_attr(ncvar_t *ncvar, struct gridaxis_t *gridaxis)
{
  strcpy(gridaxis->name, ncvar->name);
  strcpy(gridaxis->longname, ncvar->longname);
  strcpy(gridaxis->units, ncvar->units);
  if ( gridaxis->cvals )
    gridaxis->stdname = ncvar->stdname;
}

static
int cdf_get_xydimid(int ndims, int *dimids, int *dimtype, int *xdimid, int *ydimid)
{
  int nxdims = 0, nydims = 0;
  int xdimids[2] = {-1,-1}, ydimids[2] = {-1,-1};

  for ( int i = 0; i < ndims; i++ )
    {
      if ( dimtype[i] == X_AXIS && nxdims < 2 )
        {
          xdimids[nxdims] = dimids[i];
          nxdims++;
        }
      else if ( dimtype[i] == Y_AXIS && nydims < 2 )
        {
          ydimids[nydims] = dimids[i];
          nydims++;
        }
    }

  if ( nxdims == 2 )
    {
      *xdimid = xdimids[1];
      *ydimid = xdimids[0];
    }
  else if ( nydims == 2 )
    {
      *xdimid = ydimids[1];
      *ydimid = ydimids[0];
    }
  else
    {
      *xdimid = xdimids[0];
      *ydimid = ydimids[0];
    }

  return nydims;
}

static
void cdf_check_gridtype(int *gridtype, bool islon, bool islat, size_t xsize, size_t ysize, grid_t *grid)
{
  if ( islat && (islon || xsize == 0) )
    {
      double yinc = 0;
      if ( islon && (int) ysize > 1 )
        {
          yinc = fabs(grid->y.vals[0] - grid->y.vals[1]);
          for ( size_t i = 2; i < ysize; i++ )
            if ( (fabs(grid->y.vals[i-1] - grid->y.vals[i]) - yinc) > (yinc/1000) )
              {
                yinc = 0;
                break;
              }
        }
      if ( ysize < 10000 && isGaussGrid(ysize, yinc, grid->y.vals) )
        {
          *gridtype = GRID_GAUSSIAN;
          grid->np = (int)(ysize/2);
        }
      else
        *gridtype = GRID_LONLAT;
    }
  else if ( islon && !islat && ysize == 0 )
    {
      *gridtype = GRID_LONLAT;
    }
  else
    *gridtype = GRID_GENERIC;
}

static
bool cdf_read_xcoord(struct cdfLazyGrid *restrict lazyGrid, ncdim_t *ncdims, ncvar_t *ncvar, int xvarid, ncvar_t *axisvar,
                     size_t *xsize, size_t ysize, int ntdims, size_t *start, size_t *count, bool *islon)
{
  grid_t *grid = &lazyGrid->base;
  bool skipvar = true;
  *islon = axisvar->islon;
  int ndims = axisvar->ndims;
  size_t size = 0;
  int datatype = cdfInqDatatype(axisvar->xtype, axisvar->lunsigned);

  if ( (ndims - ntdims) == 2 )
    {
      /* Check size of 2 dimensional coordinate variables */
      int dimid = axisvar->dimids[ndims-2];
      size_t dimsize1 = ncdims[dimid].len;
      dimid = axisvar->dimids[ndims-1];
      size_t dimsize2 = ncdims[dimid].len;

      if ( datatype == CDI_DATATYPE_UINT8 )
        {
          ncvar->gridtype = GRID_CHARXY;
          size = dimsize1*dimsize2;
          skipvar = dimsize1 != *xsize;
        }
      else
        {
          ncvar->gridtype = GRID_CURVILINEAR;
          size = (*xsize)*ysize;
          skipvar = dimsize1*dimsize2 != size;
        }
    }
  else if ( (ndims - ntdims) == 1 )
    {
      size = *xsize;
      /* Check size of 1 dimensional coordinate variables */
      int dimid = axisvar->dimids[ndims-1];
      size_t dimsize = ncdims[dimid].len;
      skipvar = dimsize != size;
    }
  else if ( ndims == 0 && *xsize == 0 )
    {
      size = *xsize = 1;
      skipvar = false;
    }

  if ( skipvar )
    {
      Warning("Unsupported array structure, skipped variable %s!", ncvar->name);
      ncvar->isvar = -1;
      return true;
    }

  if ( datatype != -1 )  grid->datatype = datatype;

  if ( datatype == CDI_DATATYPE_UINT8 && !CDI_netcdf_lazy_grid_load )
    {
      cdf_load_cvals(size, xvarid, axisvar, &grid->x.cvals, *xsize);
      grid->x.clength = size / (*xsize) ;
    }
  else
    cdf_load_vals(size, ndims, xvarid, axisvar, &grid->x.vals, &lazyGrid->xValsGet, ntdims, start, count);

  cdf_copy_axis_attr(axisvar, &grid->x);

  return false;
}

static
bool cdf_read_ycoord(struct cdfLazyGrid *restrict lazyGrid, ncdim_t *ncdims, ncvar_t *ncvar, int yvarid, ncvar_t *axisvar,
                     size_t xsize, size_t *ysize, int ntdims, size_t *start, size_t *count, bool *islat)
{
  grid_t *grid = &lazyGrid->base;
  bool skipvar = true;
  *islat = axisvar->islat;
  int ndims = axisvar->ndims;
  size_t size = 0;
  int datatype = cdfInqDatatype(axisvar->xtype, axisvar->lunsigned);

  if ( (ndims - ntdims) == 2 )
    {
      /* Check size of 2 dimensional coordinate variables */
      int dimid = axisvar->dimids[ndims-2];
      size_t dimsize1 = ncdims[dimid].len;
      dimid = axisvar->dimids[ndims-1];
      size_t dimsize2 = ncdims[dimid].len;

      if ( datatype == CDI_DATATYPE_UINT8 )
        {
          ncvar->gridtype = GRID_CHARXY;
          size = dimsize1*dimsize2;
          skipvar = dimsize1 != *ysize;
        }
      else
        {
          ncvar->gridtype = GRID_CURVILINEAR;
          size = xsize*(*ysize);
          skipvar = dimsize1*dimsize2 != size;
        }
    }
  else if ( (ndims - ntdims) == 1 )
    {
      if ( (int) *ysize == 0 ) size = xsize;
      else                    size = *ysize;

      int dimid = axisvar->dimids[ndims-1];
      size_t dimsize = ncdims[dimid].len;
      skipvar = dimsize != size;
    }
  else if ( ndims == 0 && *ysize == 0 )
    {
      size = *ysize = 1;
      skipvar = false;
    }

  if ( skipvar )
    {
      Warning("Unsupported array structure, skipped variable %s!", ncvar->name);
      ncvar->isvar = -1;
      return true;
    }

  if ( datatype != -1 )  grid->datatype = datatype;

  if ( datatype == CDI_DATATYPE_UINT8 && !CDI_netcdf_lazy_grid_load )
    {
      cdf_load_cvals(size, yvarid, axisvar, &grid->y.cvals, *ysize);
      grid->y.clength = size / (*ysize) ;
    }
  else
    cdf_load_vals(size, ndims, yvarid, axisvar, &grid->y.vals, &lazyGrid->yValsGet, ntdims, start, count);

  cdf_copy_axis_attr(axisvar, &grid->y);

  return false;
}

static
bool cdf_read_coordinates(struct cdfLazyGrid *restrict lazyGrid, ncvar_t *ncvar, ncvar_t *ncvars, ncdim_t *ncdims,
                          int timedimid, int xvarid, int yvarid, size_t xsize, size_t ysize, int *vdimid)
{
  grid_t *grid = &lazyGrid->base;
  size_t size = 0;

  grid->datatype = CDI_DATATYPE_FLT64;

  if ( ncvar->gridtype == GRID_TRAJECTORY )
    {
      if ( ncvar->xvarid == CDI_UNDEFID ) Error("Longitude coordinate undefined for %s!", ncvar->name);
      if ( ncvar->yvarid == CDI_UNDEFID ) Error("Latitude coordinate undefined for %s!", ncvar->name);
    }
  else
    {
      static bool ltwarn = true;
      size_t start[3], count[3];
      int ntdims = 0;

      if ( xvarid != CDI_UNDEFID && yvarid != CDI_UNDEFID )
        {
          int ndims = ncvars[xvarid].ndims;
          if ( ndims != ncvars[yvarid].ndims && !ncvars[xvarid].isc && !ncvars[yvarid].isc )
            {
              Warning("Inconsistent grid structure for variable %s!", ncvar->name);
              ncvar->xvarid = xvarid = CDI_UNDEFID;
              ncvar->yvarid = yvarid = CDI_UNDEFID;
            }
          if ( ndims > 1 )
            {
              if ( ndims <= 3 )
                {
                  if ( ncvars[xvarid].dimids[0] == timedimid && ncvars[yvarid].dimids[0] == timedimid )
                    {
                      size_t ntsteps = 0;
                      cdf_inq_dimlen(ncvar->ncid, timedimid, &ntsteps);
                      if ( ltwarn && ntsteps > 1 ) Warning("Time varying grids unsupported, using grid at time step 1!");
                      ltwarn = false;
                      ntdims = 1;
                      start[0] = start[1] = start[2] = 0;
                      count[0] = 1; count[1] = ysize; count[ndims-1] = xsize;
                    }
                }
              else
                {
                  Warning("Unsupported grid structure for variable %s (grid dims > 2)!", ncvar->name);
                  ncvar->xvarid = xvarid = CDI_UNDEFID;
                  ncvar->yvarid = yvarid = CDI_UNDEFID;
                }
            }
        }

      if ( xvarid != CDI_UNDEFID )
        {
          if ( (ncvars[xvarid].ndims - ntdims) > 2 )
            {
              Warning("Coordinate variable %s has to many dimensions (%d), skipped!", ncvars[xvarid].name, ncvars[xvarid].ndims);
              //ncvar->xvarid = CDI_UNDEFID;
              xvarid = CDI_UNDEFID;
            }
        }

      if ( yvarid != CDI_UNDEFID )
        {
          if ( (ncvars[yvarid].ndims - ntdims) > 2 )
            {
              Warning("Coordinate variable %s has to many dimensions (%d), skipped!", ncvars[yvarid].name, ncvars[yvarid].ndims);
              //ncvar->yvarid = CDI_UNDEFID;
              yvarid = CDI_UNDEFID;
            }
        }

      bool islon = false, islat = false;

      if ( xvarid != CDI_UNDEFID )
        if ( cdf_read_xcoord(lazyGrid, ncdims, ncvar, xvarid, &ncvars[xvarid],
                             &xsize, ysize, ntdims, start, count, &islon) )
          return true;

      if ( yvarid != CDI_UNDEFID )
        if ( cdf_read_ycoord(lazyGrid, ncdims, ncvar, yvarid, &ncvars[yvarid],
                             xsize, &ysize, ntdims, start, count, &islat) )
          return true;

      if      ( (int) ysize == 0 ) size = xsize;
      else if ( (int) xsize == 0 ) size = ysize;
      else if ( ncvar->gridtype == GRID_UNSTRUCTURED ) size = xsize;
      else                         size = xsize*ysize;

      if ( ncvar->gridtype == CDI_UNDEFID || ncvar->gridtype == GRID_GENERIC )
        cdf_check_gridtype(&ncvar->gridtype, islon, islat, xsize, ysize, grid);
    }

  int gridtype = grid->type;
  if ( gridtype != GRID_PROJECTION ) gridtype = ncvar->gridtype;
  else if ( gridtype == GRID_PROJECTION && ncvar->gridtype == GRID_LONLAT )
    {
      int gmapvarid = ncvar->gmapid;
      if ( gmapvarid != CDI_UNDEFID && cdfCheckAttText(ncvar->ncid, gmapvarid, "grid_mapping_name") )
        {
          char attstring[CDI_MAX_NAME];
          cdfGetAttText(ncvar->ncid, gmapvarid, "grid_mapping_name", CDI_MAX_NAME, attstring);
          if ( strcmp(attstring, "latitude_longitude") == 0 ) gridtype = ncvar->gridtype;
        }
    }

  switch (gridtype)
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_UNSTRUCTURED:
    case GRID_CURVILINEAR:
    case GRID_PROJECTION:
      {
        grid->size  = (int)size;
        grid->x.size = (int)xsize;
        grid->y.size = (int)ysize;
        if ( xvarid != CDI_UNDEFID )
          {
            grid->x.flag = 1;
            int bvarid = ncvars[xvarid].bounds;
            if ( bvarid != CDI_UNDEFID )
              {
                int nbdims = ncvars[bvarid].ndims;
                if ( nbdims == 2 || nbdims == 3 )
                  {
                    *vdimid = ncvars[bvarid].dimids[nbdims-1];
                    grid->nvertex = (int)ncdims[*vdimid].len;
                    cdf_load_bounds(size*(size_t)grid->nvertex, &ncvars[xvarid], &grid->x.bounds, &lazyGrid->xBoundsGet);
                  }
              }
          }
        if ( yvarid != CDI_UNDEFID )
          {
            grid->y.flag = 1;
            int bvarid = ncvars[yvarid].bounds;
            if ( bvarid != CDI_UNDEFID )
              {
                int nbdims = ncvars[bvarid].ndims;
                if ( nbdims == 2 || nbdims == 3 )
                  {
                    if ( *vdimid == CDI_UNDEFID )
                      {
                        *vdimid = ncvars[bvarid].dimids[nbdims-1];
                        grid->nvertex = (int)ncdims[*vdimid].len;
                      }
                    cdf_load_bounds(size*(size_t)grid->nvertex, &ncvars[yvarid], &grid->y.bounds, &lazyGrid->yBoundsGet);
                  }
              }
          }

        if ( ncvar->cellarea != CDI_UNDEFID )
          cdf_load_cellarea(size, ncvar, &grid->area, &lazyGrid->cellAreaGet);

        break;
      }
    case GRID_SPECTRAL:
      {
        grid->size = (int)size;
        grid->lcomplex = 1;
        grid->trunc = ncvar->truncation;
        break;
      }
    case GRID_FOURIER:
      {
        grid->size = (int)size;
        grid->trunc = ncvar->truncation;
        break;
      }
    case GRID_TRAJECTORY:
      {
        grid->size = 1;
        break;
      }
    case GRID_CHARXY:
      {
        grid->size = (int)size;
        grid->x.size = (int)xsize;
        grid->y.size = (int)ysize;
        break;
      }
    }

  // if ( grid->type != GRID_PROJECTION && grid->type != ncvar->gridtype )
  if ( grid->type != gridtype )
    {
      // int gridtype = ncvar->gridtype;
      grid->type = gridtype;
      cdiGridTypeInit(grid, gridtype, grid->size);
    }

  if ( grid->size == 0 )
    {
      int ndims = ncvar->ndims;
      int *dimtype = ncvar->dimtype;
      if ( ndims == 0 ||
           (ndims == 1 && dimtype[0] == T_AXIS) ||
           (ndims == 1 && dimtype[0] == Z_AXIS) ||
           (ndims == 2 && dimtype[0] == T_AXIS && dimtype[1] == Z_AXIS) )
        {
          grid->type  = GRID_GENERIC;
          grid->size  = 1;
          grid->x.size = 0;
          grid->y.size = 0;
        }
      else
        {
          Warning("Unsupported grid, skipped variable %s!", ncvar->name);
          ncvar->isvar = -1;
          return true;
        }
    }

  return false;
}

static
bool cdf_set_unstructured_par(ncvar_t *ncvar, grid_t *grid, int *xdimid, int *ydimid, int number_of_grid_used, unsigned char *uuidOfHGrid)
{
  int ndims = ncvar->ndims;
  int *dimtype = ncvar->dimtype;

  int zdimid = CDI_UNDEFID;
  int xdimidx = CDI_UNDEFID, ydimidx = CDI_UNDEFID;

  for ( int i = 0; i < ndims; i++ )
    {
      if      ( dimtype[i] == X_AXIS ) xdimidx = i;
      else if ( dimtype[i] == Y_AXIS ) ydimidx = i;
      else if ( dimtype[i] == Z_AXIS ) zdimid = ncvar->dimids[i];
    }

  if ( *xdimid != CDI_UNDEFID && *ydimid != CDI_UNDEFID && zdimid == CDI_UNDEFID )
    {
      if ( grid->x.size > grid->y.size && grid->y.size < 1000 )
        {
          dimtype[ydimidx] = Z_AXIS;
          *ydimid = CDI_UNDEFID;
          grid->size  = grid->x.size;
          grid->y.size = 0;
        }
      else if ( grid->y.size > grid->x.size && grid->x.size < 1000 )
        {
          dimtype[xdimidx] = Z_AXIS;
          *xdimid = *ydimid;
          *ydimid = CDI_UNDEFID;
          grid->size  = grid->y.size;
          grid->x.size = grid->y.size;
          grid->y.size = 0;
        }
    }

  if ( grid->size != grid->x.size )
    {
      Warning("Unsupported array structure, skipped variable %s!", ncvar->name);
      ncvar->isvar = -1;
      return true;
    }

  if ( number_of_grid_used != CDI_UNDEFID ) grid->number = number_of_grid_used;
  if ( ncvar->position > 0 ) grid->position = ncvar->position;
  if ( uuidOfHGrid[0] != 0 ) memcpy(grid->uuid, uuidOfHGrid, 16);

  return false;
}

static
void cdf_read_mapping_atts(int ncid, int gmapvarid, int projID, const char *varname)
{
  if ( cdfCheckAttText(ncid, gmapvarid, "grid_mapping_name") )
    {
      char attstring[CDI_MAX_NAME];
      cdfGetAttText(ncid, gmapvarid, "grid_mapping_name", CDI_MAX_NAME, attstring);
      cdiGridDefKeyStr(projID, CDI_KEY_MAPNAME, (int)(strlen(attstring)+1), attstring);
    }
  else
    {
      Warning("Text attribute %s:grid_mapping_name missing!", varname);
    }

  int nvatts;
  cdf_inq_varnatts(ncid, gmapvarid, &nvatts);
  for ( int attnum = 0; attnum < nvatts; ++attnum )
    cdf_set_cdi_attr(ncid, gmapvarid, attnum, projID, CDI_GLOBAL);
}

static
void cdf_set_grid_to_similar_vars(ncvar_t *ncvar1, ncvar_t *ncvar2, int gridtype, int xdimid, int ydimid)
{
  if ( ncvar2->isvar == TRUE && ncvar2->gridID == CDI_UNDEFID )
    {
      int xdimid2 = CDI_UNDEFID, ydimid2 = CDI_UNDEFID, zdimid2 = CDI_UNDEFID;
      int xdimidx = CDI_UNDEFID, ydimidx = CDI_UNDEFID;
      int ndims2 = ncvar2->ndims;

      int *dimtype2 = ncvar2->dimtype;
      int *dimids2 = ncvar2->dimids;
      for ( int i = 0; i < ndims2; i++ )
        {
          if      ( dimtype2[i] == X_AXIS ) { xdimid2 = dimids2[i]; xdimidx = i; }
          else if ( dimtype2[i] == Y_AXIS ) { ydimid2 = dimids2[i]; ydimidx = i; }
          else if ( dimtype2[i] == Z_AXIS ) { zdimid2 = dimids2[i]; }
        }

      if ( ncvar2->gridtype == CDI_UNDEFID && gridtype == GRID_UNSTRUCTURED )
        {
          if ( xdimid == xdimid2 && ydimid2 != CDI_UNDEFID && zdimid2 == CDI_UNDEFID )
            {
              ncvar2->dimtype[ydimidx] = Z_AXIS;
              ydimid2 = CDI_UNDEFID;
            }

          if ( xdimid == ydimid2 && xdimid2 != CDI_UNDEFID && zdimid2 == CDI_UNDEFID )
            {
              ncvar2->dimtype[xdimidx] = Z_AXIS;
              xdimid2 = ydimid2;
              ydimid2 = CDI_UNDEFID;
            }
        }

      if ( xdimid == xdimid2 && (ydimid == ydimid2 || (xdimid == ydimid && ydimid2 == CDI_UNDEFID)) )
        {
          bool same_grid = ncvar1->xvarid == ncvar2->xvarid
                        && ncvar1->yvarid == ncvar2->yvarid
                        && ncvar1->position == ncvar2->position;
          /*
            if ( xvarid != -1 && ncvar2->xvarid != CDI_UNDEFID &&
            xvarid != ncvar2->xvarid ) same_grid = false;

            if ( yvarid != -1 && ncvar2->yvarid != CDI_UNDEFID &&
            yvarid != ncvar2->yvarid ) same_grid = false;
          */

          if ( same_grid )
            {
              if ( CDI_Debug ) Message("Same gridID %d %s", ncvar1->gridID, ncvar2->name);
              ncvar2->gridID = ncvar1->gridID;
              ncvar2->chunktype = ncvar1->chunktype;
            }
        }
    }
}

static
int cdf_define_all_grids(ncgrid_t *ncgrid, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars,
                         int timedimid, unsigned char *uuidOfHGrid, char *gridfile, int number_of_grid_used)
{
  for ( int ncvarid = 0; ncvarid < nvars; ++ncvarid )
    {
      ncvar_t *ncvar = &ncvars[ncvarid];
      if ( ncvar->isvar && ncvar->gridID == CDI_UNDEFID )
	{
          int ndims = ncvar->ndims;
          int *dimtype = ncvar->dimtype;
          int vdimid = CDI_UNDEFID;
          struct addIfNewRes projAdded = { .Id = CDI_UNDEFID, .isNew = 0 },
                             gridAdded = { .Id = CDI_UNDEFID, .isNew = 0 };
	  int xdimid = CDI_UNDEFID, ydimid = CDI_UNDEFID;
          int nydims = cdf_get_xydimid(ndims, ncvar->dimids, dimtype, &xdimid, &ydimid);

          int xaxisid = (xdimid != CDI_UNDEFID) ? ncdims[xdimid].ncvarid : CDI_UNDEFID;
          int yaxisid = (ydimid != CDI_UNDEFID) ? ncdims[ydimid].ncvarid : CDI_UNDEFID;
          int xvarid = (ncvar->xvarid != CDI_UNDEFID) ? ncvar->xvarid : xaxisid;
          int yvarid = (ncvar->yvarid != CDI_UNDEFID) ? ncvar->yvarid : yaxisid;

	  size_t xsize = (xdimid != CDI_UNDEFID) ? ncdims[xdimid].len : 0;
	  size_t ysize = (ydimid != CDI_UNDEFID) ? ncdims[ydimid].len : 0;

	  if ( ydimid == CDI_UNDEFID && yvarid != CDI_UNDEFID )
	    {
	      if ( ncvars[yvarid].ndims == 1 )
		{
		  ydimid = ncvars[yvarid].dimids[0];
		  ysize  = ncdims[ydimid].len;
		}
	    }

          if ( xsize > INT_MAX )
            {
              Warning("Size limit exceeded for x-axis dimension (limit=%d)!", INT_MAX);
              return CDI_EDIMSIZE;
            }
          if ( ysize > INT_MAX )
            {
              Warning("Size limit exceeded for y-axis dimension (limit=%d)!", INT_MAX);
              return CDI_EDIMSIZE;
            }

          int gmapvarid = ncvar->gmapid;
          bool lproj = gmapvarid != CDI_UNDEFID;

          if ( !lproj && xaxisid != CDI_UNDEFID && xaxisid != xvarid && yaxisid != CDI_UNDEFID && yaxisid != yvarid )
            {
              lproj = true;
            }

          bool lgrid = !(lproj && ncvar->xvarid == CDI_UNDEFID);

          bool lunstructured = xdimid != CDI_UNDEFID && xdimid == ydimid && nydims == 0;
	  if ( (ncvar->gridtype == CDI_UNDEFID || ncvar->gridtype == GRID_GENERIC) && lunstructured )
            ncvar->gridtype = GRID_UNSTRUCTURED;

          struct cdfLazyGrid *restrict lazyGrid = NULL, *restrict lazyProj = NULL;

          {
            int gridtype = !lgrid ? GRID_PROJECTION : ncvar->gridtype;
            if ( CDI_netcdf_lazy_grid_load )
              {
                cdfLazyGridRenew(&lazyGrid, gridtype);
                if ( lgrid && lproj ) cdfLazyGridRenew(&lazyProj, GRID_PROJECTION);
              }
            else
              {
                cdfBaseGridRenew(&lazyGrid, gridtype);
                if ( lgrid && lproj ) cdfBaseGridRenew(&lazyProj, GRID_PROJECTION);
              }
          }
          grid_t *grid = &lazyGrid->base;
          grid_t *proj = ( lgrid && lproj ) ? &lazyProj->base : NULL;

          xaxisid = (xdimid != CDI_UNDEFID) ? ncdims[xdimid].ncvarid : CDI_UNDEFID;
          yaxisid = (ydimid != CDI_UNDEFID) ? ncdims[ydimid].ncvarid : CDI_UNDEFID;


          if ( cdf_read_coordinates(lazyGrid, ncvar, ncvars, ncdims,
                                    timedimid, xvarid, yvarid, xsize, ysize, &vdimid) )
            continue;

	  if ( number_of_grid_used != CDI_UNDEFID &&
               (grid->type == CDI_UNDEFID || grid->type == GRID_GENERIC) &&
               xdimid != CDI_UNDEFID && xsize > 9999 )
            grid->type = GRID_UNSTRUCTURED;

          if ( grid->type == GRID_UNSTRUCTURED )
            if ( cdf_set_unstructured_par(ncvar, grid, &xdimid, &ydimid, number_of_grid_used, uuidOfHGrid) )
              continue;

          if ( lproj && lgrid )
            {
              int dumid;
              cdf_read_coordinates(lazyProj, ncvar, ncvars, ncdims, timedimid,
                                   xaxisid, yaxisid, xsize, ysize, &dumid);
	    }

	  if ( CDI_Debug )
	    {
	      Message("grid: type = %d, size = %d, nx = %d, ny %d",
		      grid->type, grid->size, grid->x.size, grid->y.size);
              if ( proj )
                Message("proj: type = %d, size = %d, nx = %d, ny %d",
                        proj->type, proj->size, proj->x.size, proj->y.size);
	    }


          if ( lgrid && lproj )
            {
              projAdded = cdiVlistAddGridIfNew(vlistID, proj, 2);
              grid->proj = projAdded.Id;
            }

          gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 1);
          ncvar->gridID = gridAdded.Id;

          int gridID = ncvar->gridID;

          if ( lproj && gmapvarid != CDI_UNDEFID )
            {
              int projID = lgrid ? grid->proj : gridID;
              int ncid = ncvars[gmapvarid].ncid;
              const char *gmapname = ncvars[gmapvarid].name;
              cdf_read_mapping_atts(ncid, gmapvarid, projID, gmapname);
              cdiGridDefKeyStr(projID, CDI_KEY_MAPPING, (int)(strlen(gmapname)+1), gmapname);
              gridVerifyProj(projID);
            }

          if ( grid->type == GRID_UNSTRUCTURED && gridfile[0] != 0 )
            gridDefReference(gridID, gridfile);

          if ( ncvar->chunked ) grid_set_chunktype(grid, ncvar);

	  int gridindex = vlistGridIndex(vlistID, gridID);
          ncgrid[gridindex].gridID = gridID;
          ncgrid[gridindex].ncIDs[CDF_DIMID_X] = xdimid;
          ncgrid[gridindex].ncIDs[CDF_DIMID_Y] = ydimid;
          if ( grid->type == GRID_TRAJECTORY )
            {
              ncgrid[gridindex].ncIDs[CDF_VARID_X] = xvarid;
              ncgrid[gridindex].ncIDs[CDF_VARID_Y] = yvarid;
            }

          if ( xdimid == CDI_UNDEFID && ydimid == CDI_UNDEFID && grid->size == 1 )
            gridDefHasDims(gridID, FALSE);

          if ( xdimid != CDI_UNDEFID ) cdiGridDefKeyStr(gridID, CDI_KEY_XDIMNAME, (int)(strlen(ncdims[xdimid].name)+1), ncdims[xdimid].name);
          if ( ydimid != CDI_UNDEFID ) cdiGridDefKeyStr(gridID, CDI_KEY_YDIMNAME, (int)(strlen(ncdims[ydimid].name)+1), ncdims[ydimid].name);
          if ( vdimid != CDI_UNDEFID ) cdiGridDefKeyStr(gridID, CDI_KEY_VDIMNAME, (int)(strlen(ncdims[vdimid].name)+1), ncdims[vdimid].name);

	  if ( CDI_Debug ) Message("gridID %d %d %s", gridID, ncvarid, ncvar->name);

	  for ( int ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
            cdf_set_grid_to_similar_vars(ncvar, &ncvars[ncvarid2], grid->type, xdimid, ydimid);

          if ( gridAdded.isNew ) lazyGrid = NULL;
          if ( projAdded.isNew ) lazyProj = NULL;

          if ( lazyGrid )
            {
              if ( CDI_netcdf_lazy_grid_load ) cdfLazyGridDestroy(lazyGrid);
              if ( grid ) { grid_free(grid); Free(grid); }
            }

          if ( lazyProj )
            {
              if ( CDI_netcdf_lazy_grid_load ) cdfLazyGridDestroy(lazyProj);
              if ( proj ) { grid_free(proj); Free(proj); }
            }
	}
    }

  return 0;
}

/* define all input zaxes */
static
int cdf_define_all_zaxes(stream_t *streamptr, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars,
                         size_t vctsize_echam, double *vct_echam, unsigned char *uuidOfVGrid)
{
  char *pname, *plongname, *punits;
  size_t vctsize = vctsize_echam;
  double *vct = vct_echam;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ncvar_t *ncvar = &ncvars[ncvarid];
      if ( ncvar->isvar == TRUE && ncvar->zaxisID == CDI_UNDEFID )
	{
          bool is_scalar = false;
	  bool with_bounds = false;
	  int zdimid = CDI_UNDEFID;
	  int zvarid = CDI_UNDEFID;
	  size_t zsize = 1;
          int psvarid = -1;
          int p0varid = -1;

          int positive = 0;
	  int ndims = ncvar->ndims;

          if ( ncvar->zvarid != -1 && ncvars[ncvar->zvarid].ndims == 0 )
            {
              zvarid = ncvar->zvarid;
              is_scalar = true;
            }
          else
            {
              for ( int i = 0; i < ndims; i++ )
                {
                  if ( ncvar->dimtype[i] == Z_AXIS )
                    zdimid = ncvar->dimids[i];
                }

              if ( zdimid != CDI_UNDEFID )
                {
                  // zvarid = ncdims[zdimid].ncvarid;
                  zvarid = (ncvar->zvarid != CDI_UNDEFID) ? ncvar->zvarid : ncdims[zdimid].ncvarid;
                  zsize  = ncdims[zdimid].len;
                }
            }

	  if ( CDI_Debug ) Message("nlevs = %zu", zsize);

	  double *zvar = NULL;
          char **zcvals = NULL;
          size_t zclength = 0;

	  int zaxisType = CDI_UNDEFID;
	  if ( zvarid != CDI_UNDEFID ) zaxisType = ncvars[zvarid].zaxistype;
	  if ( zaxisType == CDI_UNDEFID ) zaxisType = ZAXIS_GENERIC;

	  int zdatatype = CDI_DATATYPE_FLT64;
	  double *restrict lbounds = NULL;
	  double *restrict ubounds = NULL;

	  if ( zvarid != CDI_UNDEFID )
	    {
	      positive  = ncvars[zvarid].positive;
	      pname     = ncvars[zvarid].name;
	      plongname = ncvars[zvarid].longname;
	      punits    = ncvars[zvarid].units;
	      if ( ncvars[zvarid].xtype == NC_FLOAT ) zdatatype = CDI_DATATYPE_FLT32;
	      /* don't change the name !!! */
	      /*
	      if ( (len = strlen(pname)) > 2 )
		if ( pname[len-2] == '_' && isdigit((int) pname[len-1]) )
		  pname[len-2] = 0;
	      */
              if ( zaxisType == ZAXIS_CHAR )
                {
                  if ( ncvars[zvarid].ndims == 2 )
                    {
                      zdatatype = CDI_DATATYPE_UINT8;
                      zclength = ncdims[ncvars[zvarid].dimids[1]].len;
                      cdf_load_cvals(zsize*zclength, zvarid, ncvar, &zcvals, zsize);
                    }
                }

              if ( zaxisType == ZAXIS_HYBRID && ncvars[zvarid].vct )
                {
                  vct = ncvars[zvarid].vct;
                  vctsize = ncvars[zvarid].vctsize;

                  if ( ncvars[zvarid].psvarid != -1 ) psvarid = ncvars[zvarid].psvarid;
                  if ( ncvars[zvarid].p0varid != -1 ) p0varid = ncvars[zvarid].p0varid;
                }

              if ( zaxisType != ZAXIS_CHAR )
                {
                  zvar = (double*) Malloc(zsize*sizeof(double));
                  cdf_get_var_double(ncvars[zvarid].ncid, zvarid, zvar);
                }

	      if ( ncvars[zvarid].bounds != CDI_UNDEFID )
		{
		  int nbdims = ncvars[ncvars[zvarid].bounds].ndims;
		  if ( nbdims == 2 || is_scalar )
		    {
		      size_t nlevel  = is_scalar ? 1 : (int)ncdims[ncvars[ncvars[zvarid].bounds].dimids[0]].len;
		      int nvertex = (int)ncdims[ncvars[ncvars[zvarid].bounds].dimids[1-is_scalar]].len;
		      if ( nlevel == zsize && nvertex == 2 )
			{
			  with_bounds = true;
			  lbounds = (double *) Malloc(4 * nlevel*sizeof(double));
			  ubounds = lbounds + nlevel;
			  double *restrict zbounds = lbounds + 2 * nlevel;
			  cdf_get_var_double(ncvars[zvarid].ncid, ncvars[zvarid].bounds, zbounds);
			  for ( size_t i = 0; i < nlevel; ++i )
			    {
			      lbounds[i] = zbounds[i*2];
			      ubounds[i] = zbounds[i*2+1];
			    }
			}
		    }
		}
	    }
	  else
	    {
              pname     = (zdimid != CDI_UNDEFID) ? ncdims[zdimid].name : NULL;
	      plongname = NULL;
	      punits    = NULL;

	      if ( zsize == 1 && zdimid == CDI_UNDEFID )
		{
                  zaxisType = (ncvar->zaxistype != CDI_UNDEFID) ? ncvar->zaxistype : ZAXIS_SURFACE;
                  // if ( pname )
                    {
                      zvar = (double*) Malloc(sizeof(double));
                      zvar[0] = 0;
                    }
                }
	    }

          if ( zsize > INT_MAX )
            {
              Warning("Size limit exceeded for z-axis dimension (limit=%d)!", INT_MAX);
              return CDI_EDIMSIZE;
            }

      	  ncvar->zaxisID = varDefZaxis(vlistID, zaxisType, (int) zsize, zvar, (const char **)zcvals, zclength, with_bounds, lbounds, ubounds,
                                       (int)vctsize, vct, pname, plongname, punits, zdatatype, 1, 0);

          int zaxisID = ncvar->zaxisID;

          if ( CDI_cmor_mode && zsize == 1 && zaxisType != ZAXIS_HYBRID ) zaxisDefScalar(zaxisID);

          if ( uuidOfVGrid[0] != 0 )
            zaxisDefUUID(zaxisID, uuidOfVGrid);

          if ( zaxisType == ZAXIS_HYBRID )
            {
              if ( psvarid != -1 )
                cdiZaxisDefKeyStr(zaxisID, CDI_KEY_PSNAME, (int)(strlen(ncvars[psvarid].name)+1), ncvars[psvarid].name);
              if ( p0varid != -1 )
                {
                  double px = 1;
                  cdf_get_var_double(ncvars[p0varid].ncid, p0varid, &px);
                  cdiZaxisDefKeyFlt(zaxisID, CDI_KEY_P0VALUE, px);
                  cdiZaxisDefKeyStr(zaxisID, CDI_KEY_P0NAME, (int)(strlen(ncvars[p0varid].name)+1), ncvars[p0varid].name);
                }
            }

          if ( positive > 0 ) zaxisDefPositive(zaxisID, positive);
          if ( is_scalar ) zaxisDefScalar(zaxisID);

          if ( zdimid != CDI_UNDEFID )
            cdiZaxisDefKeyStr(zaxisID, CDI_KEY_DIMNAME, (int)(strlen(ncdims[zdimid].name)+1), ncdims[zdimid].name);
          /*
          if ( vdimid != -1 )
            cdiZaxisDefKeyStr(zaxisID, CDI_KEY_VDIMNAME, strlen(ncdims[vdimid].name)+1, ncdims[vdimid].name);
          */
	  if ( zvar    ) Free(zvar);
	  if ( zcvals  )
            {
              for ( size_t i = 0; i < zsize; i++ )
                Free(zcvals[i]);
              Free(zcvals);
            }
	  if ( lbounds ) Free(lbounds);

          if ( zvarid != CDI_UNDEFID )
            {
              int ncid = ncvars[zvarid].ncid;
              int nvatts = ncvars[zvarid].natts;
              for ( int iatt = 0; iatt < nvatts; ++iatt )
                {
                  int attnum = ncvars[zvarid].atts[iatt];
                  cdf_set_cdi_attr(ncid, zvarid, attnum, zaxisID, CDI_GLOBAL);
                }
            }

          int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
	  streamptr->zaxisID[zaxisindex] = zdimid;

	  if ( CDI_Debug )
	    Message("zaxisID %d %d %s", zaxisID, ncvarid, ncvar->name);

	  for ( int ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
	    if ( ncvars[ncvarid2].isvar == TRUE && ncvars[ncvarid2].zaxisID == CDI_UNDEFID /*&& ncvars[ncvarid2].zaxistype == CDI_UNDEFID*/ )
	      {
                int zvarid2 = CDI_UNDEFID;
                if ( ncvars[ncvarid2].zvarid != CDI_UNDEFID && ncvars[ncvars[ncvarid2].zvarid].ndims == 0 )
                  zvarid2 = ncvars[ncvarid2].zvarid;

		int zdimid2 = CDI_UNDEFID;
		ndims = ncvars[ncvarid2].ndims;
		for ( int i = 0; i < ndims; i++ )
		  {
		    if ( ncvars[ncvarid2].dimtype[i] == Z_AXIS )
		      zdimid2 = ncvars[ncvarid2].dimids[i];
		  }

		if ( zdimid == zdimid2 /* && zvarid == zvarid2 */)
		  {
                    if ( (zdimid != CDI_UNDEFID && ncvars[ncvarid2].zaxistype == CDI_UNDEFID) ||
                         (zdimid == CDI_UNDEFID && zvarid != CDI_UNDEFID && zvarid == zvarid2) ||
                         (zdimid == CDI_UNDEFID && zaxisType == ncvars[ncvarid2].zaxistype) ||
                         (zdimid == CDI_UNDEFID && zvarid2 == CDI_UNDEFID && ncvars[ncvarid2].zaxistype == CDI_UNDEFID) )
                      {
                        if ( CDI_Debug )
                          Message("zaxisID %d %d %s", zaxisID, ncvarid2, ncvars[ncvarid2].name);
                        ncvars[ncvarid2].zaxisID = zaxisID;
                      }
                  }
	      }
	}
    }

  return 0;
}


struct cdf_varinfo
{
  int        varid;
  const char *name;
};

static
int cdf_cmp_varname(const void *s1, const void *s2)
{
  const struct cdf_varinfo *x = (const struct cdf_varinfo *)s1,
                           *y = (const struct cdf_varinfo *)s2;
  return strcmp(x->name, y->name);
}

/* define all input data variables */
static
void cdf_define_all_vars(stream_t *streamptr, int vlistID, int instID, int modelID, int *varids, int nvars, int num_ncvars, ncvar_t *ncvars)
{
  if ( CDI_Debug )
    for ( int i = 0; i < nvars; i++ ) Message("varids[%d] = %d", i, varids[i]);

  if ( streamptr->sortname )
    {
      struct cdf_varinfo *varInfo
        = (struct cdf_varinfo *) Malloc((size_t)nvars * sizeof(struct cdf_varinfo));

      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int ncvarid = varids[varID];
	  varInfo[varID].varid = ncvarid;
	  varInfo[varID].name = ncvars[ncvarid].name;
	}
      qsort(varInfo, (size_t)nvars, sizeof(varInfo[0]), cdf_cmp_varname);
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  varids[varID] = varInfo[varID].varid;
	}
      Free(varInfo);
      if ( CDI_Debug )
        for ( int i = 0; i < nvars; i++ ) Message("sorted varids[%d] = %d", i, varids[i]);
    }

  for ( int varID1 = 0; varID1 < nvars; varID1++ )
    {
      int ncvarid = varids[varID1];
      int gridID  = ncvars[ncvarid].gridID;
      int zaxisID = ncvars[ncvarid].zaxisID;

      stream_new_var(streamptr, gridID, zaxisID, CDI_UNDEFID);
      int varID = vlistDefVar(vlistID, gridID, zaxisID, ncvars[ncvarid].timetype);

#if  defined  (HAVE_NETCDF4)
      if ( ncvars[ncvarid].deflate )
	vlistDefVarCompType(vlistID, varID, CDI_COMPRESS_ZIP);

      if ( ncvars[ncvarid].chunked && ncvars[ncvarid].chunktype != CDI_UNDEFID )
        vlistDefVarChunkType(vlistID, varID, ncvars[ncvarid].chunktype);
#endif

      streamptr->vars[varID1].defmiss = false;
      streamptr->vars[varID1].ncvarid = ncvarid;

      vlistDefVarName(vlistID, varID, ncvars[ncvarid].name);
      if ( ncvars[ncvarid].param != CDI_UNDEFID ) vlistDefVarParam(vlistID, varID, ncvars[ncvarid].param);
      if ( ncvars[ncvarid].code != CDI_UNDEFID )  vlistDefVarCode(vlistID, varID, ncvars[ncvarid].code);
      if ( ncvars[ncvarid].code != CDI_UNDEFID )
	{
	  int param = cdiEncodeParam(ncvars[ncvarid].code, ncvars[ncvarid].tabnum, 255);
	  vlistDefVarParam(vlistID, varID, param);
	}
      if ( ncvars[ncvarid].longname[0] )  vlistDefVarLongname(vlistID, varID, ncvars[ncvarid].longname);
      if ( ncvars[ncvarid].stdname[0] )   vlistDefVarStdname(vlistID, varID, ncvars[ncvarid].stdname);
      if ( ncvars[ncvarid].units[0] )     vlistDefVarUnits(vlistID, varID, ncvars[ncvarid].units);

      if ( ncvars[ncvarid].lvalidrange )
        vlistDefVarValidrange(vlistID, varID, ncvars[ncvarid].validrange);

      if ( IS_NOT_EQUAL(ncvars[ncvarid].addoffset, 0) )
	vlistDefVarAddoffset(vlistID, varID, ncvars[ncvarid].addoffset);
      if ( IS_NOT_EQUAL(ncvars[ncvarid].scalefactor, 1) )
	vlistDefVarScalefactor(vlistID, varID, ncvars[ncvarid].scalefactor);

      vlistDefVarDatatype(vlistID, varID, cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned));

      vlistDefVarInstitut(vlistID, varID, instID);
      vlistDefVarModel(vlistID, varID, modelID);
      if ( ncvars[ncvarid].tableID != CDI_UNDEFID )
	vlistDefVarTable(vlistID, varID, ncvars[ncvarid].tableID);

      if ( ncvars[ncvarid].deffillval == false && ncvars[ncvarid].defmissval )
        {
          ncvars[ncvarid].deffillval = true;
          ncvars[ncvarid].fillval    = ncvars[ncvarid].missval;
        }

      if ( ncvars[ncvarid].deffillval )
        vlistDefVarMissval(vlistID, varID, ncvars[ncvarid].fillval);

      if ( CDI_Debug )
	Message("varID = %d  gridID = %d  zaxisID = %d", varID,
		vlistInqVarGrid(vlistID, varID), vlistInqVarZaxis(vlistID, varID));

      int gridindex = vlistGridIndex(vlistID, gridID);
      int xdimid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      int ydimid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];

      int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
      int zdimid = streamptr->zaxisID[zaxisindex];

      int ndims = ncvars[ncvarid].ndims;
      int iodim = 0;
      int ixyz = 0;
      static const int ipow10[4] = {1, 10, 100, 1000};

      if ( ncvars[ncvarid].timetype != TIME_CONSTANT ) iodim++;

      const int *dimids = ncvars[ncvarid].dimids;

      if ( gridInqType(gridID) == GRID_UNSTRUCTURED && ndims-iodim <= 2 && ydimid == xdimid )
        {
          ixyz = (xdimid == dimids[ndims-1]) ? 321 : 213;
        }
      else
        {
          for ( int idim = iodim; idim < ndims; idim++ )
            {
              if      ( xdimid == dimids[idim] ) ixyz += 1*ipow10[ndims-idim-1];
              else if ( ydimid == dimids[idim] ) ixyz += 2*ipow10[ndims-idim-1];
              else if ( zdimid == dimids[idim] ) ixyz += 3*ipow10[ndims-idim-1];
            }
        }

      vlistDefVarXYZ(vlistID, varID, ixyz);
      /*
      printf("ixyz %d\n", ixyz);
      printf("ndims %d\n", ncvars[ncvarid].ndims);
      for ( int i = 0; i < ncvars[ncvarid].ndims; ++i )
        printf("dimids: %d %d\n", i, dimids[i]);
      printf("xdimid, ydimid %d %d\n", xdimid, ydimid);
      */
      if ( ncvars[ncvarid].ensdata != NULL )
        {
          vlistDefVarEnsemble( vlistID, varID, ncvars[ncvarid].ensdata->ens_index,
                               ncvars[ncvarid].ensdata->ens_count,
                               ncvars[ncvarid].ensdata->forecast_init_type );
          Free(ncvars[ncvarid].ensdata);
          ncvars[ncvarid].ensdata = NULL;
        }

      if ( ncvars[ncvarid].extra[0] != 0 )
        {
          vlistDefVarExtra(vlistID, varID, ncvars[ncvarid].extra);
        }
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      int ncvarid = varids[varID];
      int ncid = ncvars[ncvarid].ncid;

      int nvatts = ncvars[ncvarid].natts;
      for ( int iatt = 0; iatt < nvatts; ++iatt )
        {
          int attnum = ncvars[ncvarid].atts[iatt];
          cdf_set_cdi_attr(ncid, ncvarid, attnum, vlistID, varID);
        }

      if ( ncvars[ncvarid].atts )
        {
          Free(ncvars[ncvarid].atts);
          ncvars[ncvarid].atts = NULL;
        }

      if ( ncvars[ncvarid].vct )
        {
          Free(ncvars[ncvarid].vct);
          ncvars[ncvarid].vct = NULL;
        }
    }

  /* release mem of not freed attributes */
  for ( int ncvarid = 0; ncvarid < num_ncvars; ncvarid++ )
    if ( ncvars[ncvarid].atts ) Free(ncvars[ncvarid].atts);

  if ( varids ) Free(varids);

  for ( int varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarCode(vlistID, varID) == -varID-1 )
	{
          char name[CDI_MAX_NAME]; name[0] = 0;
	  vlistInqVarName(vlistID, varID, name);
	  size_t len = strlen(name);
	  if ( len > 3 && isdigit((int) name[3]) )
	    {
	      if ( str_is_equal(name, "var") )
		{
		  vlistDefVarCode(vlistID, varID, atoi(name+3));
                  // vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 4 && isdigit((int) name[4]) )
	    {
	      if ( str_is_equal(name, "code") )
		{
		  vlistDefVarCode(vlistID, varID, atoi(name+4));
		  // vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 5 && isdigit((int) name[5]) )
	    {
	      if ( str_is_equal(name, "param") )
		{
		  int pnum = -1, pcat = 255, pdis = 255;
		  sscanf(name+5, "%d.%d.%d", &pnum, &pcat, &pdis);
		  vlistDefVarParam(vlistID, varID, cdiEncodeParam(pnum, pcat, pdis));
                  // vlistDestroyVarName(vlistID, varID);
		}
	    }
	}
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      int varInstID  = vlistInqVarInstitut(vlistID, varID);
      int varModelID = vlistInqVarModel(vlistID, varID);
      int varTableID = vlistInqVarTable(vlistID, varID);
      int code = vlistInqVarCode(vlistID, varID);
      if ( cdiDefaultTableID != CDI_UNDEFID )
	{
          char name[CDI_MAX_NAME]; name[0] = 0;
          char longname[CDI_MAX_NAME]; longname[0] = 0;
          char units[CDI_MAX_NAME]; units[0] = 0;
          tableInqEntry(cdiDefaultTableID, code, -1, name, longname, units);
	  if ( name[0] )
	    {
	      vlistDestroyVarName(vlistID, varID);
	      vlistDestroyVarLongname(vlistID, varID);
	      vlistDestroyVarUnits(vlistID, varID);

	      if ( varTableID != CDI_UNDEFID )
		{
		  vlistDefVarName(vlistID, varID, name);
		  if ( longname[0] ) vlistDefVarLongname(vlistID, varID, longname);
		  if ( units[0] ) vlistDefVarUnits(vlistID, varID, units);
		}
	      else
		{
		  varTableID = cdiDefaultTableID;
		}
	    }

	  if ( cdiDefaultModelID != CDI_UNDEFID ) varModelID = cdiDefaultModelID;
	  if ( cdiDefaultInstID  != CDI_UNDEFID ) varInstID  = cdiDefaultInstID;
	}
      if ( varInstID  != CDI_UNDEFID ) vlistDefVarInstitut(vlistID, varID, varInstID);
      if ( varModelID != CDI_UNDEFID ) vlistDefVarModel(vlistID, varID, varModelID);
      if ( varTableID != CDI_UNDEFID ) vlistDefVarTable(vlistID, varID, varTableID);
    }
}

static
void cdf_scan_global_attr(int fileID, int vlistID, stream_t *streamptr, int ngatts, int *instID, int *modelID, bool *ucla_les, unsigned char *uuidOfHGrid, unsigned char *uuidOfVGrid, char *gridfile, int *number_of_grid_used)
{
  nc_type xtype;
  size_t attlen;
  char attname[CDI_MAX_NAME];
  char attstring[65636];

  for ( int iatt = 0; iatt < ngatts; iatt++ )
    {
      cdf_inq_attname(fileID, NC_GLOBAL, iatt, attname);
      cdf_inq_atttype(fileID, NC_GLOBAL, attname, &xtype);
      cdf_inq_attlen(fileID, NC_GLOBAL, attname, &attlen);

      if ( xtypeIsText(xtype) )
	{
	  cdfGetAttText(fileID, NC_GLOBAL, attname, sizeof(attstring), attstring);

          size_t attstrlen = strlen(attstring);

	  if ( attlen > 0 && attstring[0] != 0 )
	    {
	      if ( strcmp(attname, "history") == 0 )
		{
		  streamptr->historyID = iatt;
		}
	      else if ( strcmp(attname, "institution") == 0 )
		{
		  *instID = institutInq(0, 0, NULL, attstring);
		  if ( *instID == CDI_UNDEFID )
		    *instID = institutDef(0, 0, NULL, attstring);
		}
	      else if ( strcmp(attname, "source") == 0 )
		{
		  *modelID = modelInq(-1, 0, attstring);
		  if ( *modelID == CDI_UNDEFID )
		    *modelID = modelDef(-1, 0, attstring);
		}
	      else if ( strcmp(attname, "Source") == 0 )
		{
		  if ( strncmp(attstring, "UCLA-LES", 8) == 0 )
		    *ucla_les = true;
		}
	      /*
	      else if ( strcmp(attname, "Conventions") == 0 )
		{
		}
	      */
	      else if ( strcmp(attname, "CDI") == 0 )
		{
		}
	      else if ( strcmp(attname, "CDO") == 0 )
		{
		}
              /*
	      else if ( strcmp(attname, "forecast_reference_time") == 0 )
		{
                  memcpy(fcreftime, attstring, attstrlen+1);
		}
              */
	      else if ( strcmp(attname, "grid_file_uri") == 0 )
		{
                  memcpy(gridfile, attstring, attstrlen+1);
		}
	      else if ( strcmp(attname, "uuidOfHGrid") == 0 && attstrlen == 36 )
		{
                  attstring[36] = 0;
                  cdiStr2UUID(attstring, uuidOfHGrid);
                  //   printf("uuid: %d %s\n", attlen, attstring);
		}
	      else if ( strcmp(attname, "uuidOfVGrid") == 0 && attstrlen == 36 )
		{
                  attstring[36] = 0;
                  cdiStr2UUID(attstring, uuidOfVGrid);
		}
	      else
		{
                  if ( strcmp(attname, "ICON_grid_file_uri") == 0 && gridfile[0] == 0 )
                    memcpy(gridfile, attstring, attstrlen+1);

		  cdiDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)attstrlen, attstring);
		}
	    }
	}
      else if ( xtype == NC_SHORT || xtype == NC_INT )
	{
	  if ( strcmp(attname, "number_of_grid_used") == 0 )
	    {
	      (*number_of_grid_used) = CDI_UNDEFID;
	      cdfGetAttInt(fileID, NC_GLOBAL, attname, 1, number_of_grid_used);
	    }
 	  else
            {
              int attint[attlen];
              cdfGetAttInt(fileID, NC_GLOBAL, attname, attlen, attint);
              int datatype = (xtype == NC_SHORT) ? CDI_DATATYPE_INT16 : CDI_DATATYPE_INT32;
              cdiDefAttInt(vlistID, CDI_GLOBAL, attname, datatype, (int)attlen, attint);
            }
        }
      else if ( xtype == NC_FLOAT || xtype == NC_DOUBLE )
	{
	  double attflt[attlen];
	  cdfGetAttDouble(fileID, NC_GLOBAL, attname, attlen, attflt);
          int datatype = (xtype == NC_FLOAT) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
          cdiDefAttFlt(vlistID, CDI_GLOBAL, attname, datatype, (int)attlen, attflt);
	}
    }
}

static
int find_leadtime(int nvars, ncvar_t *ncvars)
{
  int leadtime_id = CDI_UNDEFID;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].ndims == 1 )
        if ( ncvars[ncvarid].stdname[0] && strcmp(ncvars[ncvarid].stdname, "forecast_period") == 0 )
          {
            leadtime_id = ncvarid;
            break;
          }
    }

  return leadtime_id;
}

static
void find_time_vars(int nvars, ncvar_t *ncvars, ncdim_t *ncdims, int timedimid, stream_t *streamptr,
                    bool *time_has_units, bool *time_has_bounds, bool *time_climatology)
{
  int ncvarid;

  if ( timedimid == CDI_UNDEFID )
    {
      char timeunits[CDI_MAX_NAME];

      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
        {
          if ( ncvars[ncvarid].ndims == 0 && strcmp(ncvars[ncvarid].name, "time") == 0 )
            {
              if ( ncvars[ncvarid].units[0] )
                {
                  strcpy(timeunits, ncvars[ncvarid].units);
                  str_tolower(timeunits);

                  if ( is_time_units(timeunits) )
                    {
                      streamptr->basetime.ncvarid = ncvarid;
                      break;
                    }
                }
            }
        }
    }
  else
    {
      bool ltimevar = false;

      if ( ncdims[timedimid].ncvarid != CDI_UNDEFID )
        {
          streamptr->basetime.ncvarid = ncdims[timedimid].ncvarid;
          ltimevar = true;
        }

      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
        if ( ncvarid != streamptr->basetime.ncvarid &&
             ncvars[ncvarid].ndims == 1 &&
             timedimid == ncvars[ncvarid].dimids[0] &&
             !xtypeIsText(ncvars[ncvarid].xtype) &&
             is_timeaxis_units(ncvars[ncvarid].units) )
          {
            ncvars[ncvarid].isvar = FALSE;

            if ( !ltimevar )
              {
                streamptr->basetime.ncvarid = ncvarid;
                ltimevar = true;
                if ( CDI_Debug )
                  fprintf(stderr, "timevar %s\n", ncvars[ncvarid].name);
              }
            else
              {
                Warning("Found more than one time variable, skipped variable %s!", ncvars[ncvarid].name);
              }
          }

      if ( ltimevar == false ) /* search for WRF time description */
        {
          for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
            if ( ncvarid != streamptr->basetime.ncvarid &&
                 ncvars[ncvarid].ndims == 2 &&
                 timedimid == ncvars[ncvarid].dimids[0] &&
                 xtypeIsText(ncvars[ncvarid].xtype) &&
                 ncdims[ncvars[ncvarid].dimids[1]].len == 19 )
              {
                streamptr->basetime.ncvarid = ncvarid;
                streamptr->basetime.lwrf    = true;
                break;
              }
        }

      /* time varID */
      ncvarid = streamptr->basetime.ncvarid;

      if ( ncvarid == CDI_UNDEFID )
        {
          Warning("Time variable >%s< not found!", ncdims[timedimid].name);
        }
    }

  /* time varID */
  ncvarid = streamptr->basetime.ncvarid;

  if ( ncvarid != CDI_UNDEFID && streamptr->basetime.lwrf == false )
    {
      if ( ncvars[ncvarid].units[0] != 0 ) *time_has_units = true;

      if ( ncvars[ncvarid].bounds != CDI_UNDEFID )
        {
          int nbdims = ncvars[ncvars[ncvarid].bounds].ndims;
          if ( nbdims == 2 )
            {
              int len = (int) ncdims[ncvars[ncvars[ncvarid].bounds].dimids[nbdims-1]].len;
              if ( len == 2 && timedimid == ncvars[ncvars[ncvarid].bounds].dimids[0] )
                {
                  *time_has_bounds = true;
                  streamptr->basetime.ncvarboundsid = ncvars[ncvarid].bounds;
                  if ( ncvars[ncvarid].climatology ) *time_climatology = true;
                }
            }
        }
    }
}

static
void read_vct_echam(int fileID, int nvars, ncvar_t *ncvars, ncdim_t *ncdims, double **vct, size_t *pvctsize)
{
  /* find ECHAM VCT */
  int nvcth_id = CDI_UNDEFID, vcta_id = CDI_UNDEFID, vctb_id = CDI_UNDEFID;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].ndims == 1 )
        {
          size_t len = strlen(ncvars[ncvarid].name);
          if ( len == 4 && ncvars[ncvarid].name[0] == 'h' && ncvars[ncvarid].name[1] == 'y' )
            {
              if ( ncvars[ncvarid].name[2] == 'a' && ncvars[ncvarid].name[3] == 'i' ) // hyai
                {
                  vcta_id = ncvarid;
                  nvcth_id = ncvars[ncvarid].dimids[0];
                  ncvars[ncvarid].isvar = FALSE;
                }
              else if ( ncvars[ncvarid].name[2] == 'b' && ncvars[ncvarid].name[3] == 'i' ) //hybi
                {
                  vctb_id = ncvarid;
                  nvcth_id = ncvars[ncvarid].dimids[0];
                  ncvars[ncvarid].isvar = FALSE;
                }
              else if ( (ncvars[ncvarid].name[2] == 'a' || ncvars[ncvarid].name[2] == 'b') && ncvars[ncvarid].name[3] == 'm' )
                {
                  ncvars[ncvarid].isvar = FALSE; // hyam or hybm
                }
            }
	}
    }

  /* read VCT */
  if ( nvcth_id != CDI_UNDEFID && vcta_id != CDI_UNDEFID && vctb_id != CDI_UNDEFID )
    {
      size_t vctsize = ncdims[nvcth_id].len;
      vctsize *= 2;
      *vct = (double *) Malloc(vctsize*sizeof(double));
      cdf_get_var_double(fileID, vcta_id, *vct);
      cdf_get_var_double(fileID, vctb_id, *vct+vctsize/2);
      *pvctsize = vctsize;
    }
}

static
void cdf_set_ucla_dimtype(int ndims, ncdim_t *ncdims, ncvar_t *ncvars)
{
  for ( int ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      int ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
        {
          if ( ncdims[ncdimid].dimtype == CDI_UNDEFID && ncvars[ncvarid].units[0] == 'm' )
            {
              if      ( ncvars[ncvarid].name[0] == 'x' ) ncdims[ncdimid].dimtype = X_AXIS;
              else if ( ncvars[ncvarid].name[0] == 'y' ) ncdims[ncdimid].dimtype = Y_AXIS;
              else if ( ncvars[ncvarid].name[0] == 'z' ) ncdims[ncdimid].dimtype = Z_AXIS;
            }
        }
    }
}

static
int cdf_check_vars(int nvars, ncvar_t *ncvars, size_t ntsteps, int timedimid)
{
  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( timedimid != CDI_UNDEFID )
	if ( ncvars[ncvarid].isvar == -1 &&
	     ncvars[ncvarid].ndims > 1   &&
	     timedimid == ncvars[ncvarid].dimids[0] )
	  cdf_set_var(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims == 0 )
	cdf_set_var(ncvars, ncvarid, FALSE);

      //if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims > 1 )
      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims >= 1 )
	cdf_set_var(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Variable %s has an unknown type, skipped!", ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].ndims > 4 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("%d dimensional variables are not supported, skipped variable %s!",
		ncvars[ncvarid].ndims, ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].ndims == 4 && timedimid == CDI_UNDEFID )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("%d dimensional variables without time dimension are not supported, skipped variable %s!",
		ncvars[ncvarid].ndims, ncvars[ncvarid].name);
	  continue;
	}

      if ( xtypeIsText(ncvars[ncvarid].xtype) )
	{
	  ncvars[ncvarid].isvar = 0;
	  continue;
	}

      if ( cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned) == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Unsupported data type, skipped variable %s!", ncvars[ncvarid].name);
	  continue;
	}

      if ( timedimid != CDI_UNDEFID && ntsteps == 0 && ncvars[ncvarid].ndims > 0 )
	{
	  if ( timedimid == ncvars[ncvarid].dimids[0] )
	    {
	      ncvars[ncvarid].isvar = 0;
	      Warning("Number of time steps undefined, skipped variable %s!", ncvars[ncvarid].name);
	      continue;
	    }
	}
    }

  return timedimid;
}


int cdfInqContents(stream_t *streamptr)
{
  int ndims, nvars, ngatts, unlimdimid;
  int ncvarid;
  int ncdimid;
  int *varids;
  int nvarids;
  bool time_has_units = false;
  bool time_has_bounds = false;
  bool time_climatology = false;
  int leadtime_id = CDI_UNDEFID;
  int nvars_data;
  int instID  = CDI_UNDEFID;
  int modelID = CDI_UNDEFID;
  int calendar = CDI_UNDEFID;
  int format = 0;
  bool ucla_les = false;
  char gridfile[8912];
  char fcreftime[CDI_MAX_NAME];
  int number_of_grid_used = CDI_UNDEFID;

  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  unsigned char uuidOfVGrid[CDI_UUID_SIZE];
  memset(uuidOfHGrid, 0, CDI_UUID_SIZE);
  memset(uuidOfVGrid, 0, CDI_UUID_SIZE);
  gridfile[0] = 0;
  fcreftime[0] = 0;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  if ( CDI_Debug ) Message("streamID = %d, fileID = %d", streamptr->self, fileID);

#if  defined  (HAVE_NETCDF4)
  nc_inq_format(fileID, &format);
#endif

  cdf_inq(fileID, &ndims , &nvars, &ngatts, &unlimdimid);

  if ( CDI_Debug )
    Message("root: ndims %d, nvars %d, ngatts %d", ndims, nvars, ngatts);
  /*
  if ( ndims == 0 )
    {
      Warning("No dimensions found!");
      return CDI_EUFSTRUCT;
    }
  */
  // alloc ncdims
  ncdim_t *ncdims = ndims ? (ncdim_t *) Malloc((size_t)ndims * sizeof(ncdim_t)) : NULL;
  init_ncdims(ndims, ncdims);

#if  defined  (TEST_GROUPS)
#if  defined  (HAVE_NETCDF4)
  if ( format == NC_FORMAT_NETCDF4 )
    {
      int ncid;
      int numgrps;
      int ncids[NC_MAX_VARS];
      char name1[CDI_MAX_NAME];
      int gndims, gnvars, gngatts, gunlimdimid;
      nc_inq_grps(fileID, &numgrps, ncids);
      for ( int i = 0; i < numgrps; ++i )
        {
          ncid = ncids[i];
          nc_inq_grpname(ncid, name1);
          cdf_inq(ncid, &gndims , &gnvars, &gngatts, &gunlimdimid);

          if ( CDI_Debug )
            Message("%s: ndims %d, nvars %d, ngatts %d", name1, gndims, gnvars, gngatts);

          if ( gndims == 0 )
            {
            }
        }
    }
#endif
#endif

  if ( nvars == 0 )
    {
      Warning("No arrays found!");
      return CDI_EUFSTRUCT;
    }

  // alloc ncvars
  ncvar_t *ncvars = nvars ? (ncvar_t *) Malloc((size_t)nvars * sizeof (ncvar_t)) : NULL;
  init_ncvars(nvars, ncvars);

  for ( ncvarid = 0; ncvarid < nvars; ++ncvarid ) ncvars[ncvarid].ncid = fileID;


  // scan global attributes
  cdf_scan_global_attr(fileID, vlistID, streamptr, ngatts, &instID, &modelID, &ucla_les,
                       uuidOfHGrid, uuidOfVGrid, gridfile, &number_of_grid_used);

  // find time dim
  int timedimid = (unlimdimid >= 0) ? unlimdimid : cdf_time_dimid(fileID, ndims, nvars);

  streamptr->basetime.ncdimid = timedimid;

  size_t ntsteps = 0;
  if ( timedimid != CDI_UNDEFID ) cdf_inq_dimlen(fileID, timedimid, &ntsteps);
  if ( ntsteps > INT_MAX )
    {
      Warning("Size limit exceeded for time dimension (limit=%d)!", INT_MAX);
      return CDI_EDIMSIZE;
    }

  if ( CDI_Debug ) Message("Number of timesteps = %zu", ntsteps);
  if ( CDI_Debug ) Message("Time dimid = %d", streamptr->basetime.ncdimid);

  // read ncdims
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      cdf_inq_dimlen(fileID, ncdimid, &ncdims[ncdimid].len);
      cdf_inq_dimname(fileID, ncdimid, ncdims[ncdimid].name);
      if ( timedimid == ncdimid )
	ncdims[ncdimid].dimtype = T_AXIS;
    }

  if ( CDI_Debug ) cdf_print_vars(ncvars, nvars, "cdf_scan_var_attr");

  // scan attributes of all variables
  cdf_scan_var_attr(nvars, ncvars, ncdims, timedimid, modelID, format);


  if ( CDI_Debug ) cdf_print_vars(ncvars, nvars, "find coordinate vars");

  // find coordinate vars
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
	{
	  if ( ncvars[ncvarid].ndims == 1 )
	    {
	      if ( timedimid != CDI_UNDEFID && timedimid == ncvars[ncvarid].dimids[0] )
		{
		  if ( ncvars[ncvarid].isvar != FALSE ) cdf_set_var(ncvars, ncvarid, TRUE);
		}
	      else
		{
                  //  if ( ncvars[ncvarid].isvar != TRUE ) cdf_set_var(ncvars, ncvarid, FALSE);
		}
	      // if ( ncvars[ncvarid].isvar != TRUE ) cdf_set_var(ncvars, ncvarid, FALSE);

	      if ( ncdimid == ncvars[ncvarid].dimids[0] && ncdims[ncdimid].ncvarid == CDI_UNDEFID )
		if ( strcmp(ncvars[ncvarid].name, ncdims[ncdimid].name) == 0 )
		  {
		    ncdims[ncdimid].ncvarid = ncvarid;
		    ncvars[ncvarid].isvar = FALSE;
		  }
	    }
	}
    }

  // find time vars
  find_time_vars(nvars, ncvars, ncdims, timedimid, streamptr, &time_has_units, &time_has_bounds, &time_climatology);

  leadtime_id = find_leadtime(nvars, ncvars);
  if ( leadtime_id != CDI_UNDEFID ) ncvars[leadtime_id].isvar = FALSE;

  // check ncvars
  timedimid = cdf_check_vars(nvars, ncvars, ntsteps, timedimid);

  // verify coordinate vars - first scan (dimname == varname)
  bool lhybrid_cf = false;
  verify_coordinate_vars_1(fileID, ndims, ncdims, ncvars, timedimid, &lhybrid_cf);

  // verify coordinate vars - second scan (all other variables)
  verify_coordinate_vars_2(nvars, ncvars);

  if ( CDI_Debug ) cdf_print_vars(ncvars, nvars, "verify_coordinate_vars");

  if ( ucla_les ) cdf_set_ucla_dimtype(ndims, ncdims, ncvars);

  /*
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  printf("coord var %d %s %s\n", ncvarid, ncvars[ncvarid].name, ncvars[ncvarid].units);
	  if ( ncdims[ncdimid].dimtype == X_AXIS )
	    printf("coord var %d %s is x dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Y_AXIS )
	    printf("coord var %d %s is y dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Z_AXIS )
	    printf("coord var %d %s is z dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == T_AXIS )
	    printf("coord var %d %s is t dim\n", ncvarid, ncvars[ncvarid].name);

	  if ( ncvars[ncvarid].islon )
	    printf("coord var %d %s is lon\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islat )
	    printf("coord var %d %s is lat\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islev )
	    printf("coord var %d %s is lev\n", ncvarid, ncvars[ncvarid].name);
	}
    }
  */

  // Set coordinate varids (att: associate)
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ncvar_t *ncvar = &ncvars[ncvarid];
      if ( ncvar->isvar == TRUE && ncvar->ncoordvars )
	{
	  int ncoordvars = ncvar->ncoordvars;
	  for ( int i = 0; i < ncoordvars; i++ )
	    {
	      if      ( ncvars[ncvar->coordvarids[i]].islon ||
                        ncvars[ncvar->coordvarids[i]].isx )   ncvar->xvarid = ncvar->coordvarids[i];
	      else if ( ncvars[ncvar->coordvarids[i]].islat ||
                        ncvars[ncvar->coordvarids[i]].isy )   ncvar->yvarid = ncvar->coordvarids[i];
	      else if ( ncvars[ncvar->coordvarids[i]].islev ) ncvar->zvarid = ncvar->coordvarids[i];
	      else if ( ncvars[ncvar->coordvarids[i]].isc )   ncvar->cvarids[i] = ncvar->coordvarids[i];
	    }
	}
    }

  // set dim type
  cdf_set_dimtype(nvars, ncvars, ncdims);

  // read ECHAM VCT if present
  size_t vctsize = 0;
  double *vct = NULL;
  if ( !lhybrid_cf ) read_vct_echam(fileID, nvars, ncvars, ncdims, &vct, &vctsize);


  if ( CDI_Debug ) cdf_print_vars(ncvars, nvars, "cdf_define_all_grids");

  // define all grids
  int status;
  status = cdf_define_all_grids(streamptr->ncgrid, vlistID, ncdims, nvars, ncvars, timedimid, uuidOfHGrid, gridfile, number_of_grid_used);
  if ( status < 0 ) return status;

  // define all zaxes
  status = cdf_define_all_zaxes(streamptr, vlistID, ncdims, nvars, ncvars, vctsize, vct, uuidOfVGrid);
  if ( vct ) Free(vct);
  if ( status < 0 ) return status;


  // select vars
  varids = (int *) Malloc((size_t)nvars * sizeof (int));
  nvarids = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    if ( ncvars[ncvarid].isvar == TRUE ) varids[nvarids++] = ncvarid;

  nvars_data = nvarids;

  if ( CDI_Debug ) Message("time varid = %d", streamptr->basetime.ncvarid);
  if ( CDI_Debug ) Message("ntsteps = %zu", ntsteps);
  if ( CDI_Debug ) Message("nvars_data = %d", nvars_data);


  if ( nvars_data == 0 )
    {
      streamptr->ntsteps = 0;
      return CDI_EUFSTRUCT;
    }

  if ( ntsteps == 0 && streamptr->basetime.ncdimid == CDI_UNDEFID && streamptr->basetime.ncvarid != CDI_UNDEFID )
    ntsteps = 1;

  streamptr->ntsteps = (long)ntsteps;

  // define all data variables
  cdf_define_all_vars(streamptr, vlistID, instID, modelID, varids, nvars_data, nvars, ncvars);


  cdiCreateTimesteps(streamptr);

  // time varID
  int nctimevarid = streamptr->basetime.ncvarid;

  if ( time_has_units )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;

      if ( setBaseTime(ncvars[nctimevarid].units, taxis) == 1 )
        {
          nctimevarid = CDI_UNDEFID;
          streamptr->basetime.ncvarid = CDI_UNDEFID;
        }

      if ( leadtime_id != CDI_UNDEFID && taxis->type == TAXIS_RELATIVE )
        {
          streamptr->basetime.leadtimeid = leadtime_id;
          taxis->type = TAXIS_FORECAST;

          int timeunit = -1;
          if ( ncvars[leadtime_id].units[0] != 0 ) timeunit = scanTimeUnit(ncvars[leadtime_id].units);
          if ( timeunit == -1 ) timeunit = taxis->unit;
          taxis->fc_unit = timeunit;

          setForecastTime(fcreftime, taxis);
        }
    }

  if ( time_has_bounds )
    {
      streamptr->tsteps[0].taxis.has_bounds = true;
      if ( time_climatology ) streamptr->tsteps[0].taxis.climatology = true;
    }

  if ( nctimevarid != CDI_UNDEFID )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;
      ptaxisDefName(taxis, ncvars[nctimevarid].name);

      if ( ncvars[nctimevarid].longname[0] )
        ptaxisDefLongname(taxis, ncvars[nctimevarid].longname);

      if ( ncvars[nctimevarid].units[0] )
        ptaxisDefUnits(taxis, ncvars[nctimevarid].units);

      int datatype = (ncvars[nctimevarid].xtype == NC_FLOAT) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
      ptaxisDefDatatype(taxis, datatype);
    }

  if ( nctimevarid != CDI_UNDEFID )
    if ( ncvars[nctimevarid].calendar == true )
      {
        char attstring[1024];
	cdfGetAttText(fileID, nctimevarid, "calendar", sizeof(attstring), attstring);
	str_tolower(attstring);
        set_calendar(attstring, &calendar);
      }

  int taxisID;
  if ( streamptr->tsteps[0].taxis.type == TAXIS_FORECAST )
    {
      taxisID = taxisCreate(TAXIS_FORECAST);
    }
  else if ( streamptr->tsteps[0].taxis.type == TAXIS_RELATIVE )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      if ( !time_has_units )
	{
	  taxisDefTunit(taxisID, TUNIT_DAY);
	  streamptr->tsteps[0].taxis.unit = TUNIT_DAY;
	}
    }


  if ( calendar == CDI_UNDEFID && streamptr->tsteps[0].taxis.type != TAXIS_ABSOLUTE )
    {
      calendar = CALENDAR_STANDARD;
    }

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wstrict-overflow"
#endif
  if ( calendar != CDI_UNDEFID )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;
      taxis->calendar = calendar;
      taxisDefCalendar(taxisID, calendar);
    }
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif

  vlistDefTaxis(vlistID, taxisID);

  streamptr->curTsID = 0;
  streamptr->rtsteps = 1;

  (void) cdfInqTimestep(streamptr, 0);

  cdfCreateRecords(streamptr, 0);

  // free ncdims
  if ( ncdims ) Free(ncdims);

  // free ncvars
  if ( ncvars ) Free(ncvars);

  return 0;
}

static
void wrf_read_timestep(int fileID, int nctimevarid, int tsID, taxis_t *taxis)
{
  size_t start[2], count[2];
  char stvalue[32];
  start[0] = (size_t) tsID; start[1] = 0;
  count[0] = 1; count[1] = 19;
  stvalue[0] = 0;
  cdf_get_vara_text(fileID, nctimevarid, start, count, stvalue);
  stvalue[19] = 0;
  {
    int year = 1, month = 1, day = 1 , hour = 0, minute = 0, second = 0;
    if ( strlen(stvalue) == 19 )
      sscanf(stvalue, "%d-%d-%d_%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
    taxis->vdate = cdiEncodeDate(year, month, day);
    taxis->vtime = cdiEncodeTime(hour, minute, second);
    taxis->type = TAXIS_ABSOLUTE;
  }
}

static
double get_timevalue(int fileID, int nctimevarid, int tsID, timecache_t *tcache)
{
  double timevalue = 0;

  if ( tcache )
    {
      if ( tcache->size == 0 || (tsID < tcache->startid || tsID > (tcache->startid+tcache->size-1)) )
        {
          int maxvals = MAX_TIMECACHE_SIZE;
          tcache->startid = (tsID/MAX_TIMECACHE_SIZE)*MAX_TIMECACHE_SIZE;
          if ( (tcache->startid + maxvals) > tcache->maxvals ) maxvals = (tcache->maxvals)%MAX_TIMECACHE_SIZE;
          tcache->size = maxvals;
          size_t index = (size_t) tcache->startid;
          // fprintf(stderr, "fill time cache: %d %d %d %d %d\n", tcache->maxvals, tsID, tcache->startid, tcache->startid+maxvals-1, maxvals);
          for ( int ival = 0; ival < maxvals; ++ival )
            {
              cdf_get_var1_double(fileID, nctimevarid, &index, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;
              tcache->cache[ival] = timevalue;
              index++;
            }
        }

      timevalue = tcache->cache[tsID%MAX_TIMECACHE_SIZE];
    }
  else
    {
      size_t index = (size_t) tsID;
      cdf_get_var1_double(fileID, nctimevarid, &index, &timevalue);
      if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;
    }

  return timevalue;
}


int cdfInqTimestep(stream_t * streamptr, int tsID)
{
  if ( CDI_Debug ) Message("streamID = %d  tsID = %d", streamptr->self, tsID);

  if ( tsID < 0 ) Error("unexpected tsID = %d", tsID);

  if ( tsID < streamptr->ntsteps && streamptr->ntsteps > 0 )
    {
      cdfCreateRecords(streamptr, tsID);

      taxis_t *taxis = &streamptr->tsteps[tsID].taxis;
      if ( tsID > 0 )
	ptaxisCopy(taxis, &streamptr->tsteps[0].taxis);

      double timevalue = tsID;

      int nctimevarid = streamptr->basetime.ncvarid;
      if ( nctimevarid != CDI_UNDEFID )
	{
	  int fileID = streamptr->fileID;
	  size_t index = (size_t)tsID;

	  if ( streamptr->basetime.lwrf )
	    {
              wrf_read_timestep(fileID, nctimevarid, tsID, taxis);
	    }
	  else
	    {
#if defined (USE_TIMECACHE)
              if ( streamptr->basetime.timevar_cache == NULL )
                {
                  streamptr->basetime.timevar_cache = (timecache_t *) Malloc(MAX_TIMECACHE_SIZE*sizeof(timecache_t));
                  streamptr->basetime.timevar_cache->size = 0;
                  streamptr->basetime.timevar_cache->maxvals = streamptr->ntsteps;
                }
#endif
              timevalue = get_timevalue(fileID, nctimevarid, tsID, streamptr->basetime.timevar_cache);
	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate, &taxis->vtime);
	    }

	  int nctimeboundsid = streamptr->basetime.ncvarboundsid;
	  if ( nctimeboundsid != CDI_UNDEFID )
	    {
	      size_t start[2], count[2];
              start[0] = index; count[0] = 1; start[1] = 0; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_lb, &taxis->vtime_lb);

              start[0] = index; count[0] = 1; start[1] = 1; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_ub, &taxis->vtime_ub);
	    }

          int leadtimeid = streamptr->basetime.leadtimeid;
          if ( leadtimeid != CDI_UNDEFID )
            {
              timevalue = get_timevalue(fileID, leadtimeid, tsID, NULL);
              cdiSetForecastPeriod(timevalue, taxis);
            }
	}
    }

  streamptr->curTsID = tsID;
  long nrecs = streamptr->tsteps[tsID].nrecs;

  return (int) nrecs;
}


int cdfInqHistorySize(stream_t *streamptr)
{
  size_t size = 0;
  int ncid = streamptr->fileID;
  if ( streamptr->historyID != CDI_UNDEFID )
    cdf_inq_attlen(ncid, NC_GLOBAL, "history", &size);

  return (int) size;
}


void cdfInqHistoryString(stream_t *streamptr, char *history)
{
  int ncid = streamptr->fileID;
  if ( streamptr->historyID != CDI_UNDEFID )
    {
      nc_type atttype;
      cdf_inq_atttype(ncid, NC_GLOBAL, "history", &atttype);

      if ( atttype == NC_CHAR )
        {
          cdf_get_att_text(ncid, NC_GLOBAL, "history", history);
        }
#if  defined  (HAVE_NETCDF4)
      else if ( atttype == NC_STRING )
        {
          // ToDo
          Warning("History attribute with type NC_STRING unsupported!");
        }
#endif
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
