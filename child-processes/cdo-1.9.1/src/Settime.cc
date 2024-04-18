/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Settime    setdate         Set date
      Settime    settime         Set time
      Settime    setday          Set day
      Settime    setmon          Set month
      Settime    setyear         Set year
      Settime    settunits       Set time units
      Settime    settaxis        Set time axis
      Settime    setreftime      Set reference time
      Settime    setcalendar     Set calendar
      Settime    shifttime       Shift timesteps
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "calendar.h"
#include "pstream.h"


int get_tunits(const char *unit, int *incperiod, int *incunit, int *tunit)
{
  size_t len = strlen(unit);
  
  if      ( memcmp(unit, "seconds", len) == 0 ) { *incunit =     1; *tunit = TUNIT_SECOND;  }
  else if ( memcmp(unit, "minutes", len) == 0 ) { *incunit =    60; *tunit = TUNIT_MINUTE;  }
  else if ( memcmp(unit, "hours", len)   == 0 ) { *incunit =  3600; *tunit = TUNIT_HOUR;    }
  else if ( memcmp(unit, "3hours", len)  == 0 ) { *incunit = 10800; *tunit = TUNIT_3HOURS;  }
  else if ( memcmp(unit, "6hours", len)  == 0 ) { *incunit = 21600; *tunit = TUNIT_6HOURS;  }
  else if ( memcmp(unit, "12hours", len) == 0 ) { *incunit = 43200; *tunit = TUNIT_12HOURS; }
  else if ( memcmp(unit, "days", len)    == 0 ) { *incunit = 86400; *tunit = TUNIT_DAY;     }
  else if ( memcmp(unit, "months", len)  == 0 ) { *incunit =     1; *tunit = TUNIT_MONTH;   }
  else if ( memcmp(unit, "years", len)   == 0 ) { *incunit =    12; *tunit = TUNIT_YEAR;    }
  else cdoAbort("Time unit >%s< unsupported!", unit);

  if ( *tunit == TUNIT_HOUR )
    {
      if      ( *incperiod ==  3 ) { *incperiod = 1; *incunit = 10800; *tunit = TUNIT_3HOURS;  }
      else if ( *incperiod ==  6 ) { *incperiod = 1; *incunit = 21600; *tunit = TUNIT_6HOURS;  }
      else if ( *incperiod == 12 ) { *incperiod = 1; *incunit = 43200; *tunit = TUNIT_12HOURS; }
    }

  return 0;
}

static
void shifttime(int calendar, int tunit, int ijulinc, int *pdate, int *ptime)
{
  int vdate = *pdate;
  int vtime = *ptime;

  if ( tunit == TUNIT_MONTH || tunit == TUNIT_YEAR )
    {
      int year, month, day;
      cdiDecodeDate(vdate, &year, &month, &day);
	      
      month += ijulinc;

      while ( month > 12 ) { month -= 12; year++; }
      while ( month <  1 ) { month += 12; year--; }

      vdate = cdiEncodeDate(year, month, day);

      *pdate = vdate;
    }
  else
    {
      juldate_t juldate;
      juldate = juldate_encode(calendar, vdate, vtime);
      juldate = juldate_add_seconds(ijulinc, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      *pdate = vdate;
      *ptime = vtime;

      if ( cdoVerbose )
	cdoPrint("juldate, ijulinc, vdate, vtime: %g %d %d %d",
		 juldate_to_seconds(juldate), ijulinc, vdate, vtime);
    }
}

static
void gen_bounds(int calendar, int tunit, int incperiod, int vdate, int vtime, int *vdateb, int *vtimeb)
{
  juldate_t juldate;
  int year, month, day;
  
  vdateb[0] = vdate;
  vdateb[1] = vdate;
  vtimeb[0] = 0;
  vtimeb[1] = 0;
          
  cdiDecodeDate(vdate, &year, &month, &day);
  if ( tunit == TUNIT_MONTH )
    {
      vdateb[0] = cdiEncodeDate(year, month, 1);
      month++;
      if ( month > 12 ) { month = 1; year++; }
      vdateb[1] = cdiEncodeDate(year, month, 1);
    }
  else if ( tunit == TUNIT_YEAR )
    {
      vdateb[0] = cdiEncodeDate(year,   1, 1);
      vdateb[1] = cdiEncodeDate(year+1, 1, 1);
    }
  else if ( tunit == TUNIT_DAY )
    {
      vdateb[0] = vdate;
      juldate = juldate_encode(calendar, vdateb[0], vtimeb[0]);
      juldate = juldate_add_seconds(86400, juldate);
      juldate_decode(calendar, juldate, &vdateb[1], &vtimeb[1]);
    }
  else if ( tunit == TUNIT_HOUR || tunit == TUNIT_3HOURS || tunit == TUNIT_6HOURS || tunit == TUNIT_12HOURS )
    {
      if ( incperiod == 0 ) incperiod = 1;
      if ( incperiod > 24 ) cdoAbort("Time period must be less equal 24!");

      if      ( tunit == TUNIT_3HOURS  ) incperiod = 3;
      else if ( tunit == TUNIT_6HOURS  ) incperiod = 6;
      else if ( tunit == TUNIT_12HOURS ) incperiod = 12;

      int hour, minute, second;
      cdiDecodeTime(vtime, &hour, &minute, &second);
      int h0 = (hour/incperiod) * incperiod;
      vtimeb[0] = cdiEncodeTime(h0, 0, 0);
      int h1 = h0 + incperiod;
      if ( h1 >= 24 )
        {
          vdateb[1] = cdiEncodeDate(year, month, day+1);
          juldate = juldate_encode(calendar, vdateb[0], vtimeb[0]);
          juldate = juldate_add_seconds(incperiod*3600, juldate);
          juldate_decode(calendar, juldate, &vdateb[1], &vtimeb[1]);
        }
      else
        vtimeb[1] = cdiEncodeTime(h1, 0, 0);
    }
}


void *Settime(void *argument)
{
  int nrecs, newval = 0;
  int varID, levelID;
  int vdateb[2], vtimeb[2];
  int sdate = 0, stime = 0;
  int taxisID2 = CDI_UNDEFID;
  int nmiss;
  int gridsize;
  int tunit = TUNIT_DAY;
  int ijulinc = 0, incperiod = 1, incunit = 86400;
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  int day0;
  bool copy_timestep = false;
  int newcalendar = CALENDAR_STANDARD;
  // int nargs;
  juldate_t juldate;

  cdoInitialize(argument);

  // clang-format off
  int SETYEAR     = cdoOperatorAdd("setyear",      0,  1, "year");
  int SETMON      = cdoOperatorAdd("setmon",       0,  1, "month");
  int SETDAY      = cdoOperatorAdd("setday",       0,  1, "day");
  int SETDATE     = cdoOperatorAdd("setdate",      0,  1, "date (format: YYYY-MM-DD)");
  int SETTIME     = cdoOperatorAdd("settime",      0,  1, "time (format: hh:mm:ss)");
  int SETTUNITS   = cdoOperatorAdd("settunits",    0,  1, "time units (seconds, minutes, hours, days, months, years)");
  int SETTAXIS    = cdoOperatorAdd("settaxis",     0, -2, "date<,time<,increment>> (format YYYY-MM-DD,hh:mm:ss)");
  int SETTBOUNDS  = cdoOperatorAdd("settbounds",   0,  1, "frequency (day, month, year)");
  int SETREFTIME  = cdoOperatorAdd("setreftime",   0, -2, "date<,time<,units>> (format YYYY-MM-DD,hh:mm:ss)");
  int SETCALENDAR = cdoOperatorAdd("setcalendar",  0,  1, "calendar (standard, proleptic_gregorian, 360_day, 365_day, 366_day)");
  int SHIFTTIME   = cdoOperatorAdd("shifttime",    0,  1, "shift value");
  // clang-format on

  int operatorID = cdoOperatorID();
  // nargs = cdoOperatorF2(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));

  //  if ( operatorArgc()

  if ( operatorID == SETTAXIS || operatorID == SETREFTIME )
    {
      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 3 ) cdoAbort("Too many arguments!");

      const char *datestr = operatorArgv()[0];
      const char *timestr = operatorArgv()[1];

      if ( strchr(datestr+1, '-') )
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  sdate = cdiEncodeDate(year, month, day);
	}
      else
	{
	  sdate = *datestr ? parameter2int(datestr) : 10101;
	}

      if ( operatorArgc() > 1 )
        {
          if ( strchr(timestr, ':') )
            {
              sscanf(timestr, "%d:%d:%d", &hour, &minute, &second);
              stime = cdiEncodeTime(hour, minute, second);
            }
          else
            {
              stime = parameter2int(timestr);
            }

          if ( operatorArgc() == 3 )
            {
              const char *timeunits = operatorArgv()[2];
              int ich = timeunits[0];
              if ( ich == '-' || ich == '+' || isdigit(ich) )
                {
                  incperiod = (int)strtol(timeunits, NULL, 10);
                  if ( ich == '-' || ich == '+' ) timeunits++;
                  while ( isdigit((int) *timeunits) ) timeunits++;
                }
              get_tunits(timeunits, &incperiod, &incunit, &tunit);
            }
        }

      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETDATE )
    {
      operatorCheckArgc(1);
      const char *datestr = operatorArgv()[0];
      if ( strchr(datestr, '-') )
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  newval = cdiEncodeDate(year, month, day);
	}
      else
	{
	  newval = parameter2int(datestr);
	}
    }
  else if ( operatorID == SETTIME )
    {
      operatorCheckArgc(1);
      const char *timestr = operatorArgv()[0];

      if ( strchr(timestr, ':') )
	{
	  sscanf(timestr, "%d:%d:%d", &hour, &minute, &second);
	  newval = cdiEncodeTime(hour, minute, second);
	}
      else
	{
	  newval = parameter2int(timestr);
	}
    }
  else if ( operatorID == SHIFTTIME )
    {
      operatorCheckArgc(1);
      const char *timeunits = operatorArgv()[0];
      incperiod = (int)strtol(timeunits, NULL, 10);
      if ( timeunits[0] == '-' || timeunits[0] == '+' ) timeunits++;
      while ( isdigit((int) *timeunits) ) timeunits++;

      get_tunits(timeunits, &incperiod, &incunit, &tunit);

      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETTUNITS || operatorID == SETTBOUNDS )
    {
      operatorCheckArgc(1);
      const char *timeunits = operatorArgv()[0];
      incperiod = (int)strtol(timeunits, NULL, 10);
      if ( timeunits[0] == '-' || timeunits[0] == '+' ) timeunits++;
      while ( isdigit((int) *timeunits) ) timeunits++;

      get_tunits(timeunits, &incperiod, &incunit, &tunit);

      if ( operatorID == SETTBOUNDS &&
           !(tunit == TUNIT_HOUR || tunit == TUNIT_3HOURS || tunit == TUNIT_6HOURS || tunit == TUNIT_12HOURS ||
             tunit == TUNIT_DAY || tunit == TUNIT_MONTH || tunit == TUNIT_YEAR) )
        cdoAbort("Unsupported frequency %s! Use hour, 3hours, 6hours, day, month or year.", timeunits);
    }
  else if ( operatorID == SETCALENDAR )
    {
      operatorCheckArgc(1);
      char *cname = operatorArgv()[0];
      strtolower(cname);
      size_t len = strlen(cname);
      if ( len < 3 ) len = 7;
      if      ( strcmp(cname, "standard")  == 0 ) newcalendar = CALENDAR_STANDARD;
      else if ( strcmp(cname, "gregorian") == 0 ) newcalendar = CALENDAR_GREGORIAN;
      else if ( strcmp(cname, "proleptic") == 0 ) newcalendar = CALENDAR_PROLEPTIC;
      else if ( strcmp(cname, "proleptic_gregorian") == 0 ) newcalendar = CALENDAR_PROLEPTIC;
      else if ( strncmp(cname, "360days", len) == 0 ) newcalendar = CALENDAR_360DAYS;
      else if ( strncmp(cname, "360_day", len) == 0 ) newcalendar = CALENDAR_360DAYS;
      else if ( strncmp(cname, "365days", len) == 0 ) newcalendar = CALENDAR_365DAYS;
      else if ( strncmp(cname, "365_day", len) == 0 ) newcalendar = CALENDAR_365DAYS;
      else if ( strncmp(cname, "366days", len) == 0 ) newcalendar = CALENDAR_366DAYS;
      else if ( strncmp(cname, "366_day", len) == 0 ) newcalendar = CALENDAR_366DAYS;
      else cdoAbort("Calendar >%s< unsupported! Available %s", cname, cdoOperatorEnter(operatorID));
    }
  else
    {
      operatorCheckArgc(1);
      newval = parameter2int(operatorArgv()[0]);
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  bool taxis_has_bounds = taxisHasBounds(taxisID1) > 0;
  int ntsteps  = vlistNtsteps(vlistID1);
  int nvars    = vlistNvars(vlistID1);

  if ( ntsteps == 1 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	if ( vlistInqVarTimetype(vlistID1, varID) != TIME_CONSTANT ) break;

      if ( varID == nvars ) ntsteps = 0;
    }

  if ( ntsteps == 0 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
    }

  int calendar = taxisInqCalendar(taxisID1);

  if ( cdoVerbose ) cdoPrint("calendar = %d", calendar);

  if ( operatorID == SETREFTIME )
    {
      copy_timestep = true;

      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{
	  cdoPrint("Changing absolute to relative time axis!");

	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);

      if ( operatorArgc() != 3 ) tunit = taxisInqTunit(taxisID1);
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETTUNITS )
    {
      copy_timestep = true;

      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{
	  cdoPrint("Changing absolute to relative time axis!");

	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	  taxisDefTunit(taxisID2, tunit);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);
    }
  else if ( operatorID == SETCALENDAR )
    {
      copy_timestep = true;
      /*
      if ( ((char *)argument)[0] == '-' )
	cdoAbort("This operator does not work with pipes!");
      */
      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{/*
	  if ( cdoFiletype() != CDI_FILETYPE_NC )
	    cdoAbort("This operator does not work on an absolute time axis!");
	 */
	  cdoPrint("Changing absolute to relative time axis!");
	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);
    }
  else
    taxisID2 = taxisDuplicate(taxisID1);

  if ( operatorID == SETTAXIS )
    {
      taxisDefTunit(taxisID2, tunit);
      taxisDefRdate(taxisID2, sdate);
      taxisDefRtime(taxisID2, stime);
      juldate = juldate_encode(calendar, sdate, stime);
    }
  else if ( operatorID == SETTUNITS )
    {
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETCALENDAR )
    {
      taxisDefCalendar(taxisID2, newcalendar);
    }
  else if ( operatorID == SETTBOUNDS )
    {
      taxisWithBounds(taxisID2);
    }

  if ( operatorID != SHIFTTIME )
    if ( taxis_has_bounds && !copy_timestep )
      {
	cdoWarning("Time bounds unsupported by this operator, removed!");
	taxisDeleteBounds(taxisID2);
	taxis_has_bounds = false;
      }

  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID1 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      if ( operatorID == SETTAXIS )
	{
	  if ( tunit == TUNIT_MONTH || tunit == TUNIT_YEAR )
	    {
	      vtime = stime;
	      if ( tsID1 == 0 )
		{
		  vdate = sdate;
		  cdiDecodeDate(vdate, &year, &month, &day0);
		}
	      else
		{	      
		  month += ijulinc;

		  while ( month > 12 ) { month -= 12; year++; }
		  while ( month <  1 ) { month += 12; year--; }

                  day = (day0 == 31) ? days_per_month(calendar, year, month) : day0;

		  vdate = cdiEncodeDate(year, month, day);
		}
	    }
	  else
	    {
	      juldate_decode(calendar, juldate, &vdate, &vtime);
	      juldate = juldate_add_seconds(ijulinc, juldate);
	    }
	}
      else if ( operatorID == SETTBOUNDS )
	{
          gen_bounds(calendar, tunit, incperiod, vdate, vtime, vdateb, vtimeb);
          
          if ( CDO_CMOR_Mode )
            {
              juldate_t juldate1 = juldate_encode(calendar, vdateb[0], vtimeb[0]);
              juldate_t juldate2 = juldate_encode(calendar, vdateb[1], vtimeb[1]);
              double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
              juldate_t juldatem = juldate_add_seconds((int)lround(seconds), juldate1);
              juldate_decode(calendar, juldatem, &vdate, &vtime);
            }
	}
      else if ( operatorID == SHIFTTIME )
	{
	  shifttime(calendar, tunit, ijulinc, &vdate, &vtime);
	  if ( taxis_has_bounds )
	    {
	      taxisInqVdateBounds(taxisID1, &vdateb[0], &vdateb[1]);
	      taxisInqVtimeBounds(taxisID1, &vtimeb[0], &vtimeb[1]);	      
	      shifttime(calendar, tunit, ijulinc, &vdateb[0], &vtimeb[0]);
	      shifttime(calendar, tunit, ijulinc, &vdateb[1], &vtimeb[1]);
	    }
	}
      else if ( operatorID == SETREFTIME || operatorID == SETCALENDAR || operatorID == SETTUNITS )
	{
	  ;
	}
      else
	{
	  cdiDecodeDate(vdate, &year, &month, &day);

	  if ( operatorID == SETYEAR ) year  = newval;
	  if ( operatorID == SETMON  ) month = newval;
	  if ( operatorID == SETMON && (month < 0 || month > 16) ) cdoAbort("parameter month=%d out of range!", month);
	  if ( operatorID == SETDAY  ) day   = newval;
	  if ( operatorID == SETDAY && (day < 0 || day > 31) ) cdoAbort("parameter day=%d %d out of range!", day);
      
	  vdate = cdiEncodeDate(year, month, day);

	  if ( operatorID == SETDATE  ) vdate = newval;
	  if ( operatorID == SETTIME  ) vtime = newval;
	}

      if ( copy_timestep )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);
	  if ( operatorID == SETREFTIME )
	    {
	      taxisDefRdate(taxisID2, sdate);
	      taxisDefRtime(taxisID2, stime);
	    }
	}
      else
	{
	  int numavg = taxisInqNumavg(taxisID1);
	  taxisDefNumavg(taxisID2, numavg);

	  taxisDefVdate(taxisID2, vdate);
	  taxisDefVtime(taxisID2, vtime);

	  if ( taxis_has_bounds || operatorID == SETTBOUNDS )
	    {
	      taxisDefVdateBounds(taxisID2, vdateb[0], vdateb[1]);
	      taxisDefVtimeBounds(taxisID2, vtimeb[0], vtimeb[1]);
	    }
	}

      pstreamDefTimestep(streamID2, tsID1);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);
	  
	  pstreamReadRecord(streamID1, array, &nmiss);
	  pstreamWriteRecord(streamID2, array, nmiss);
	}
      
      tsID1++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}
