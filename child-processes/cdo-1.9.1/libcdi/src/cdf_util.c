#include <string.h>
#include <ctype.h>
#include "dmemory.h"
#include "cdi.h"
#include "cdf_util.h"
#include "error.h"


void str_tolower(char *str)
{
  if ( str )
    for ( size_t i = 0; str[i]; ++i )
      str[i] = (char)tolower((int)str[i]);
}


bool str_is_equal(const char *vstr, const char *cstr)
{
  bool is_equal = false;
  size_t clen = (cstr != NULL) ? strlen(cstr) : 0;

  if ( vstr && *vstr ) is_equal = (memcmp(vstr, cstr, clen) == 0);

  return is_equal;
}

int get_timeunit(size_t len, const char *ptu)
{
  int timeunit = -1;

  while ( isspace(*ptu) && len ) { ptu++; len--; }

  if ( len > 2 )
    {
      if      ( str_is_equal(ptu, "sec") )            timeunit = TUNIT_SECOND;
      else if ( str_is_equal(ptu, "minute") )         timeunit = TUNIT_MINUTE;
      else if ( str_is_equal(ptu, "hour") )           timeunit = TUNIT_HOUR;
      else if ( str_is_equal(ptu, "day") )            timeunit = TUNIT_DAY;
      else if ( str_is_equal(ptu, "month") )          timeunit = TUNIT_MONTH;
      else if ( str_is_equal(ptu, "calendar_month") ) timeunit = TUNIT_MONTH;
      else if ( str_is_equal(ptu, "year") )           timeunit = TUNIT_YEAR;
    }
  else if ( len == 1 && ptu[0] == 's' )  timeunit = TUNIT_SECOND;

  return timeunit;
}


bool is_time_units(const char *timeunits)
{
  while ( isspace(*timeunits) ) timeunits++; 

  bool status = str_is_equal(timeunits, "sec")
             || str_is_equal(timeunits, "minute")
             || str_is_equal(timeunits, "hour")
             || str_is_equal(timeunits, "day")
             || str_is_equal(timeunits, "month")
             || str_is_equal(timeunits, "calendar_month")
             || str_is_equal(timeunits, "year");

  return status;
}


bool is_timeaxis_units(const char *timeunits)
{
  bool status = false;

  size_t len = strlen(timeunits);
  char *tu = (char *) Malloc((len+1)*sizeof(char));
  memcpy(tu, timeunits, (len+1) * sizeof(char));
  char *ptu = tu;

  for ( size_t i = 0; i < len; i++ ) ptu[i] = (char)tolower((int)ptu[i]);

  int timeunit = get_timeunit(len, ptu);
  if ( timeunit != -1 )
    {
      while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
      if ( *ptu )
        {
          while ( isspace(*ptu) ) ptu++;

          int timetype = str_is_equal(ptu, "as") ? TAXIS_ABSOLUTE :
                         str_is_equal(ptu, "since") ? TAXIS_RELATIVE : -1;

          status = timetype != -1;
        }
    }

  Free(tu);

  return status;
}


bool is_height_units(const char *units)
{
  int u0 = units[0];

  bool status
    = (u0=='m' && (!units[1] || strncmp(units, "meter", 5) == 0))
    || (!units[2] && units[1]=='m' && (u0=='c' || u0=='d' || u0=='k'));

  return status;
}


bool is_pressure_units(const char *units)
{
  bool status = false;

  if ( strncmp(units, "millibar", 8) == 0 ||
       strncmp(units, "mb", 2)       == 0 ||
       strncmp(units, "hectopas", 8) == 0 ||
       strncmp(units, "hPa", 3)      == 0 ||
       strncmp(units, "Pa", 2)       == 0 )
    {
      status = true;
    }

  return status;
}


bool is_DBL_axis(/*const char *units,*/ const char *longname)
{
  bool status = false;

  if ( strcmp(longname, "depth below land")         == 0 ||
       strcmp(longname, "depth_below_land")         == 0 ||
       strcmp(longname, "levels below the surface") == 0 )
    {
      /*
      if ( strcmp(ncvars[ncvarid].units, "cm") == 0 ||
           strcmp(ncvars[ncvarid].units, "dm") == 0 ||
           strcmp(ncvars[ncvarid].units, "m")  == 0 )
      */
        status = true;
    }

  return status;
}


bool is_depth_axis(const char *stdname, const char *longname)
{
  bool status = false;

  if ( strcmp(stdname, "depth") == 0 )
    status = true;
  else
    if ( strcmp(longname, "depth_below_sea") == 0 ||
         strcmp(longname, "depth below sea") == 0 )
      {
        status = true;
      }

  return status;
}


bool is_height_axis(const char *stdname, const char *longname)
{
  bool status = false;

  if ( strcmp(stdname, "height") == 0 )
    status = true;
  else
    if ( strcmp(longname, "height") == 0 ||
         strcmp(longname, "height above the surface") == 0 )
      {
        status = true;
      }

  return status;
}


bool is_lon_axis(const char *units, const char *stdname)
{
  bool status = false;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  str_tolower(lc_units);

  if ( (str_is_equal(lc_units, "degree") || str_is_equal(lc_units, "radian")) &&
       (str_is_equal(stdname, "grid_longitude") || str_is_equal(stdname, "longitude")) )
    {
      status = true;
    }
  else if ( str_is_equal(lc_units, "degree")
            && !str_is_equal(stdname, "grid_latitude")
            && !str_is_equal(stdname, "latitude") )
    {
      int ioff = 6;
      if ( lc_units[ioff] == 's' ) ioff++;
      if ( lc_units[ioff] == '_' ) ioff++;
      if ( lc_units[ioff] == 'e' ) status = true;
    }

  return status;
}


bool is_lat_axis(const char *units, const char *stdname)
{
  bool status = false;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  str_tolower(lc_units);

  if ( (str_is_equal(lc_units, "degree") || str_is_equal(lc_units, "radian")) &&
        (str_is_equal(stdname, "grid_latitude") || str_is_equal(stdname, "latitude")) )
    {
      status = true;
    }
  else if ( str_is_equal(lc_units, "degree")
            && !str_is_equal(stdname, "grid_longitude")
            && !str_is_equal(stdname, "longitude") )
    {
      int ioff = 6;
      if ( lc_units[ioff] == 's' ) ioff++;
      if ( lc_units[ioff] == '_' ) ioff++;
      if ( lc_units[ioff] == 'n' || lc_units[ioff] == 's' ) status = true;
    }

  return status;
}


bool is_x_axis(const char *units, const char *stdname)
{
  (void)units;
  return (strcmp(stdname, "projection_x_coordinate") == 0);
}


bool is_y_axis(const char *units, const char *stdname)
{
  (void)units;
  return (strcmp(stdname, "projection_y_coordinate") == 0);
}


void set_gridtype(const char *attstring, int *gridtype)
{
  if      ( strcmp(attstring, "gaussian reduced") == 0 )
    *gridtype = GRID_GAUSSIAN_REDUCED;
  else if ( strcmp(attstring, "gaussian") == 0 )
    *gridtype = GRID_GAUSSIAN;
  else if ( strncmp(attstring, "spectral", 8) == 0 )
    *gridtype = GRID_SPECTRAL;
  else if ( strncmp(attstring, "fourier", 7) == 0 )
    *gridtype = GRID_FOURIER;
  else if ( strcmp(attstring, "trajectory") == 0 )
    *gridtype = GRID_TRAJECTORY;
  else if ( strcmp(attstring, "generic") == 0 )
    *gridtype = GRID_GENERIC;
  else if ( strcmp(attstring, "cell") == 0 )
    *gridtype = GRID_UNSTRUCTURED;
  else if ( strcmp(attstring, "unstructured") == 0 )
    *gridtype = GRID_UNSTRUCTURED;
  else if ( strcmp(attstring, "curvilinear") == 0 )
    *gridtype = GRID_CURVILINEAR;
  else if ( strcmp(attstring, "characterxy") == 0 )
    *gridtype = GRID_CHARXY;
  else if ( strcmp(attstring, "sinusoidal") == 0 )
    ;
  else if ( strcmp(attstring, "laea") == 0 )
    ;
  else if ( strcmp(attstring, "lcc2") == 0 )
    ;
  else if ( strcmp(attstring, "linear") == 0 ) // ignore grid type linear
    ;
  else
    {
      static bool warn = true;
      if ( warn )
        {
          warn = false;
          Warning("NetCDF attribute grid_type='%s' unsupported!", attstring);
        }
    }
}


void set_zaxistype(const char *attstring, int *zaxistype)
{
  if      ( strcmp(attstring, "toa") == 0 ) *zaxistype = ZAXIS_TOA;
  else if ( strcmp(attstring, "cloudbase") == 0 ) *zaxistype = ZAXIS_CLOUD_BASE;
  else if ( strcmp(attstring, "cloudtop") == 0 ) *zaxistype = ZAXIS_CLOUD_TOP;
  else if ( strcmp(attstring, "isotherm0") == 0 ) *zaxistype = ZAXIS_ISOTHERM_ZERO;
  else if ( strcmp(attstring, "seabottom") == 0 ) *zaxistype = ZAXIS_SEA_BOTTOM;
  else if ( strcmp(attstring, "lakebottom") == 0 ) *zaxistype = ZAXIS_LAKE_BOTTOM;
  else if ( strcmp(attstring, "sedimentbottom") == 0 ) *zaxistype = ZAXIS_SEDIMENT_BOTTOM;
  else if ( strcmp(attstring, "sedimentbottomta") == 0 ) *zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;
  else if ( strcmp(attstring, "sedimentbottomtw") == 0 ) *zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;
  else if ( strcmp(attstring, "mixlayer") == 0 ) *zaxistype = ZAXIS_MIX_LAYER;
  else if ( strcmp(attstring, "atmosphere") == 0 ) *zaxistype = ZAXIS_ATMOSPHERE;
  else
    {
      static bool warn = true;
      if ( warn )
        {
          warn = false;
          Warning("NetCDF attribute level_type='%s' unsupported!", attstring);
        }
    }  
}


void set_calendar(const char *attstring, int *calendar)
{
  if ( str_is_equal(attstring, "standard") )
    *calendar = CALENDAR_STANDARD;
  else if ( str_is_equal(attstring, "gregorian") )
    *calendar = CALENDAR_GREGORIAN;
  else if ( str_is_equal(attstring, "none") )
    *calendar = CALENDAR_NONE;
  else if ( str_is_equal(attstring, "proleptic") )
    *calendar = CALENDAR_PROLEPTIC;
  else if ( str_is_equal(attstring, "360") )
    *calendar = CALENDAR_360DAYS;
  else if ( str_is_equal(attstring, "365") ||
            str_is_equal(attstring, "noleap") )
    *calendar = CALENDAR_365DAYS;
  else if ( str_is_equal(attstring, "366") ||
            str_is_equal(attstring, "all_leap") )
    *calendar = CALENDAR_366DAYS;
  else
    Warning("calendar >%s< unsupported!", attstring);
}
