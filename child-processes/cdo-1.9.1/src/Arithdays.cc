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

      Arithdays  muldpm          Multiply with days per month
      Arithdays  divdpm          Divide by days per month
      Arithdays  muldpy          Multiply with days per year
      Arithdays  divdpy          Divide by days per year
      Arithdays  muldoy          Multiply with day of year
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "calendar.h"
#include "pstream.h"


static
double dayofyear(int calendar, int vdate, int vtime)
{
  int month_360[12] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
  int month_365[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int month_366[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int *dpm;
  int year, month, day;
  int hour, minute, second;
  double doy = 0;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  int dpy = days_per_year(calendar, year);

  for ( int im = 1; im < month; ++im )
    {
      if      ( dpy == 360 ) dpm = month_360;
      else if ( dpy == 365 ) dpm = month_365;
      else                   dpm = month_366;

      if ( im >= 1 && im <= 12 ) doy += dpm[im-1];
    }

  doy += (day-1);
  doy += (second+minute*60+hour*3600)/86400.;

  if ( cdoVerbose )
    cdoPrint("vdate, vtime, dpy, doy: %d %d %d %g", vdate, vtime, dpy, doy);

  return doy;
}


void *Arithdays(void *argument)
{
  int nrecs;
  int varID, levelID;
  int year, month, day;
  int nmiss;
  double rconst;

  cdoInitialize(argument);

  // clang-format off
               cdoOperatorAdd("muldpm", func_mul, func_month, NULL);
               cdoOperatorAdd("divdpm", func_div, func_month, NULL);
               cdoOperatorAdd("muldpy", func_mul, func_year,  NULL);
               cdoOperatorAdd("divdpy", func_div, func_year,  NULL);
  int MULDOY = cdoOperatorAdd("muldoy", func_mul,         0,  NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int operfunc2 = cdoOperatorF2(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int calendar = taxisInqCalendar(taxisID1);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr    = (double*) Malloc(gridsize*sizeof(double));
  field.weight = NULL;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( operatorID == MULDOY )
	{
	  rconst = dayofyear(calendar, vdate, vtime);
	}
      else
	{
	  if ( operfunc2 == func_month )
	    rconst = days_per_month(calendar, year, month);
	  else
	    rconst = days_per_year(calendar, year);
	}

      if ( cdoVerbose )
	cdoPrint("calendar %d  year %d  month %d  result %g", calendar, year, month, rconst);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field.ptr, &nmiss);

          field.nmiss   = (size_t)nmiss;
	  field.grid    = vlistInqVarGrid(vlistID1, varID);
	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  farcfun(&field, rconst, operfunc);

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, field.ptr, (int)field.nmiss);
	}
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( field.ptr ) Free(field.ptr);

  cdoFinish();

  return 0;
}
