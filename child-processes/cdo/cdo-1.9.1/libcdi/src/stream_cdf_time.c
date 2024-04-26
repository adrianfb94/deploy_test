#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include <stdio.h>
#include <string.h>

#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"
#include "cdf_int.h"

static
int cdfDefTimeBounds(int fileID, int nctimevarid, int nctimedimid, const char *taxis_name, taxis_t* taxis)
{
  int time_bndsid = -1;
  int dims[2];

  dims[0] = nctimedimid;

  /* fprintf(stderr, "time has bounds\n"); */
  static const char bndsName[] = "bnds";
  if ( nc_inq_dimid(fileID, bndsName, &dims[1]) != NC_NOERR )
    cdf_def_dim(fileID, bndsName, 2, &dims[1]);

  const char *bndsAttName, *bndsAttVal;
  size_t bndsAttValLen;
  char tmpstr[CDI_MAX_NAME];
  if ( taxis->climatology )
    {
      static const char climatology_bndsName[] = "climatology_bnds",
        climatology_bndsAttName[] = "climatology";
      bndsAttName = climatology_bndsAttName;
      bndsAttValLen = sizeof (climatology_bndsName) - 1;
      bndsAttVal = climatology_bndsName;
    }
  else
    {
      size_t taxisnameLen = strlen(taxis_name);
      memcpy(tmpstr, taxis_name, taxisnameLen);
      tmpstr[taxisnameLen] = '_';
      memcpy(tmpstr + taxisnameLen + 1, bndsName, sizeof (bndsName));
      size_t tmpstrLen = taxisnameLen + sizeof (bndsName);
      static const char generic_bndsAttName[] = "bounds";
      bndsAttName = generic_bndsAttName;
      bndsAttValLen = tmpstrLen;
      bndsAttVal = tmpstr;
    }
  cdf_def_var(fileID, bndsAttVal, NC_DOUBLE, 2, dims, &time_bndsid);
  cdf_put_att_text(fileID, nctimevarid, bndsAttName, bndsAttValLen, bndsAttVal);

  return time_bndsid;
}

static
void cdfDefTimeUnits(char *unitstr, taxis_t *taxis0, taxis_t *taxis)
{
  if ( taxis->units && taxis->units[0] )
    {
      strcpy(unitstr, taxis->units);
    }
  else
    {
      unitstr[0] = 0;

      if ( taxis0->type == TAXIS_ABSOLUTE )
        {
          static const char *const unitstrfmt[3]
            = { "year as %Y.%f",
                "month as %Y%m.%f",
                "day as %Y%m%d.%f" };
          size_t fmtidx = (taxis0->unit == TUNIT_YEAR ? 0
                           : (taxis0->unit == TUNIT_MONTH ? 1
                              : 2));
          strcpy(unitstr, unitstrfmt[fmtidx]);
        }
      else
        {
          int year, month, day, hour, minute, second;
          cdiDecodeDate(taxis->rdate, &year, &month, &day);
          cdiDecodeTime(taxis->rtime, &hour, &minute, &second);

          int timeunit = taxis->unit  != -1 ? taxis->unit  : TUNIT_HOUR;
          if      ( timeunit == TUNIT_QUARTER   ) timeunit = TUNIT_MINUTE;
          else if ( timeunit == TUNIT_30MINUTES ) timeunit = TUNIT_MINUTE;
          else if (    timeunit == TUNIT_3HOURS
                    || timeunit == TUNIT_6HOURS
                    || timeunit == TUNIT_12HOURS ) timeunit = TUNIT_HOUR;

          sprintf(unitstr, "%s since %d-%d-%d %02d:%02d:%02d",
                  tunitNamePtr(timeunit), year, month, day, hour, minute, second);
        }
    }
}

static
void cdfDefForecastTimeUnits(char *unitstr, int timeunit)
{
  unitstr[0] = 0;

  if ( timeunit == -1 ) timeunit = TUNIT_HOUR;
  else if ( timeunit == TUNIT_QUARTER   ) timeunit = TUNIT_MINUTE;
  else if ( timeunit == TUNIT_30MINUTES ) timeunit = TUNIT_MINUTE;
  else if (    timeunit == TUNIT_3HOURS
            || timeunit == TUNIT_6HOURS
            || timeunit == TUNIT_12HOURS ) timeunit = TUNIT_HOUR;

  strcpy(unitstr, tunitNamePtr(timeunit));
}

static
void cdfDefCalendar(int fileID, int ncvarid, int calendar)
{
  static const struct { int calCode; const char *calStr; } calTab[] = {
    { CALENDAR_STANDARD, "standard" },
    { CALENDAR_GREGORIAN, "gregorian" },
    { CALENDAR_PROLEPTIC, "proleptic_gregorian" },
    { CALENDAR_NONE, "none" },
    { CALENDAR_360DAYS, "360_day" },
    { CALENDAR_365DAYS, "365_day" },
    { CALENDAR_366DAYS, "366_day" },
  };
  enum { calTabSize = sizeof calTab / sizeof calTab[0] };

  for ( size_t i = 0; i < calTabSize; ++i )
    if ( calTab[i].calCode == calendar )
      {
        const char *calstr = calTab[i].calStr;
        size_t len = strlen(calstr);
        cdf_put_att_text(fileID, ncvarid, "calendar", len, calstr);
        break;
      }
}


void cdfDefTime(stream_t* streamptr)
{
  int time_varid;
  int time_dimid;
  int time_bndsid = -1;
  static const char default_name[] = "time";

  if ( streamptr->basetime.ncvarid != CDI_UNDEFID ) return;

  int fileID = streamptr->fileID;

  if ( streamptr->ncmode == 0 ) streamptr->ncmode = 1;
  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  taxis_t *taxis = &streamptr->tsteps[0].taxis;

  const char *taxis_name = (taxis->name && taxis->name[0]) ? taxis->name : default_name ;

  cdf_def_dim(fileID, taxis_name, NC_UNLIMITED, &time_dimid);
  streamptr->basetime.ncdimid = time_dimid;

  nc_type xtype = (taxis->datatype == CDI_DATATYPE_FLT32) ? NC_FLOAT : NC_DOUBLE;

  cdf_def_var(fileID, taxis_name, xtype, 1, &time_dimid, &time_varid);

  streamptr->basetime.ncvarid = time_varid;

#if  defined  (HAVE_NETCDF4)
  if ( streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C )
    {
      size_t chunk = 512;
      cdf_def_var_chunking(fileID, time_varid, NC_CHUNKED, &chunk);
    }
#endif

  {
    static const char timeStr[] = "time";
    cdf_put_att_text(fileID, time_varid, "standard_name", sizeof(timeStr) - 1, timeStr);
  }

  if ( taxis->longname && taxis->longname[0] )
    cdf_put_att_text(fileID, time_varid, "long_name", strlen(taxis->longname), taxis->longname);

  if ( taxis->has_bounds )
    {
      time_bndsid = cdfDefTimeBounds(fileID, time_varid, time_dimid, taxis_name, taxis);
      streamptr->basetime.ncvarboundsid = time_bndsid;
    }

  {
    char unitstr[CDI_MAX_NAME];
    cdfDefTimeUnits(unitstr, &streamptr->tsteps[0].taxis, taxis);
    size_t len = strlen(unitstr);
    if ( len )
      {
        cdf_put_att_text(fileID, time_varid, "units", len, unitstr);
        /*
          if ( taxis->has_bounds )
          cdf_put_att_text(fileID, time_bndsid, "units", len, unitstr);
        */
      }
  }

  if ( taxis->calendar != -1 )
    {
      cdfDefCalendar(fileID, time_varid, taxis->calendar);
      /*
      if ( taxis->has_bounds )
        cdfDefCalendar(fileID, time_bndsid, taxis->calendar);
      */
    }

  if ( taxis->type == TAXIS_FORECAST )
    {
      int leadtimeid;
      cdf_def_var(fileID, "leadtime", xtype, 1, &time_dimid, &leadtimeid);
      streamptr->basetime.leadtimeid = leadtimeid;

      {
        static const char stdname[] = "forecast_period";
        cdf_put_att_text(fileID, leadtimeid, "standard_name", sizeof(stdname) - 1, stdname);
      }

      {
        static const char lname[] = "Time elapsed since the start of the forecast";
        cdf_put_att_text(fileID, leadtimeid, "long_name", sizeof(lname) - 1, lname);
      }

      {
          char unitstr[CDI_MAX_NAME];
          cdfDefForecastTimeUnits(unitstr, taxis->fc_unit);
          size_t len = strlen(unitstr);
          if ( len )
            cdf_put_att_text(fileID, leadtimeid, "units", len, unitstr);
      }
    }

  cdf_put_att_text(fileID, time_varid, "axis", 1, "T");

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
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
