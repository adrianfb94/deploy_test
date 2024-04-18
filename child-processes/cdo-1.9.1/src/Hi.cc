/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Hi      hi           Compute the humidity index
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static const char HI_NAME[]     = "hum_index";
static const char HI_LONGNAME[] = "Humindex describes empirically in units of temperature how the temperature and humidity influence the wellness of a human being. HI = T + 5/9 * (A - 10) with A = e * (6.112 * 10 ** ((7.5 * T)/(237.7 + T)) * R), T  = air temperature in degree Celsius, R = relative humidity, e = vapour pressure. Humindex is only defined for temperatures of at least 26 degree Celsius and relative humidity of at least 40 percent.";
static const char HI_UNITS[]    = "Celsius";

static const int FIRST_VAR = 0;


static double humidityIndex(double t, double e, double r, double missval)
{
  static const double tmin = 26.0;
  static const double rmin = 40.0;
  
  if ( t < tmin || r < rmin )
    return missval;
    
  return t + (5.0 / 9.0) * ((0.01 * r * e * 6.112 * pow(10.0, (7.5 * t) / (237.7 + t))) - 10.0);
}


static void farexpr(field_type *field1, field_type field2, field_type field3, double (*expression)(double, double, double, double))
{
  const int     grid1    = field1->grid;
  const int     nmiss1   = field1->nmiss;
  const double  missval1 = field1->missval;
  double       *array1   = field1->ptr;
  const int     grid2    = field2.grid;
  const int     nmiss2   = field2.nmiss;
  const double  missval2 = field2.missval;
  const double *array2   = field2.ptr;
  const int     grid3    = field3.grid;
  const int     nmiss3   = field3.nmiss;
  const double  missval3 = field3.missval;
  const double *array3   = field3.ptr;

  int len = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) || len != gridInqSize(grid3) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 || nmiss3 > 0 )
    {
      for ( int i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) || DBL_IS_EQUAL(array3[i], missval3))  
	  array1[i] = missval1;
	else
	  array1[i] = expression(array1[i], array2[i], array3[i], missval1);
    }
  else
    {
      for ( int i = 0; i < len; i++ )
        array1[i] = expression(array1[i], array2[i], array3[i], missval1);  
    }

  field1->nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}

   
void *Hi(void *argument)
{
  int nmiss;
  int nrecs;
  int varID1, varID2, varID3;
  int levelID1, levelID2, levelID3;

  cdoInitialize(argument);
  cdoOperatorAdd("hi", 0, 0, NULL);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));
  int streamID3 = pstreamOpenRead(cdoStreamName(2));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = pstreamInqVlist(streamID3);

  int taxisID1 = vlistInqTaxis(vlistID1);
  //taxisID2 = vlistInqTaxis(vlistID2);
  //taxisID3 = vlistInqTaxis(vlistID3);

  vlistCompare(vlistID1, vlistID2, CMP_DIM);
  vlistCompare(vlistID1, vlistID3, CMP_DIM);
  
  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field2, field3;
  field_init(&field1);
  field_init(&field2);
  field_init(&field3);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));
  field3.ptr = (double*) Malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d, file3 %d",
	     vlistNtsteps(vlistID1), vlistNtsteps(vlistID2), vlistNtsteps(vlistID3));

  int vlistID4 = vlistCreate();
  int gridID   = vlistInqVarGrid(vlistID1, FIRST_VAR);
  int zaxisID  = vlistInqVarZaxis(vlistID1, FIRST_VAR);
  int varID4   = vlistDefVar(vlistID4, gridID, zaxisID, TIME_VARYING);

  int taxisID4 = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(taxisID4, TUNIT_MINUTE);
  taxisDefCalendar(taxisID4, CALENDAR_STANDARD);
  taxisDefRdate(taxisID4, 19550101);
  taxisDefRtime(taxisID4, 0);
  vlistDefTaxis(vlistID4, taxisID4);

  vlistDefVarName(vlistID4, varID4, HI_NAME);
  vlistDefVarLongname(vlistID4, varID4, HI_LONGNAME);
  vlistDefVarUnits(vlistID4, varID4, HI_UNITS);

  int streamID4 = pstreamOpenWrite(cdoStreamName(3), cdoFiletype());

  pstreamDefVlist(streamID4, vlistID4);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int nrecs2 = pstreamInqTimestep(streamID2, tsID);
      int nrecs3 = pstreamInqTimestep(streamID3, tsID);
      if ( nrecs2 == 0 || nrecs3 == 0 )
        cdoAbort("Input streams have different number of timesteps!");

      taxisCopyTimestep(taxisID4, taxisID1);
      pstreamDefTimestep(streamID4, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID1, &levelID1);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = (size_t) nmiss;
          
	  pstreamInqRecord(streamID2, &varID2, &levelID2);
	  pstreamReadRecord(streamID2, field2.ptr, &nmiss);
          field2.nmiss = (size_t) nmiss;
	  
	  pstreamInqRecord(streamID3, &varID3, &levelID3);
	  pstreamReadRecord(streamID3, field3.ptr, &nmiss);
          field3.nmiss = (size_t) nmiss;
	  
	  if ( varID1 != varID2 || varID1 != varID3 || levelID1 != levelID2 || levelID1 != levelID3 )
	    cdoAbort("Input streams have different structure!");
	    
          if ( varID1 != FIRST_VAR ) continue;
            
	  field1.grid    = vlistInqVarGrid(vlistID1, varID1);
	  field1.missval = vlistInqVarMissval(vlistID1, varID1);

	  field2.grid    = vlistInqVarGrid(vlistID2, varID2);
	  field2.missval = vlistInqVarMissval(vlistID2, varID2);

	  field3.grid    = vlistInqVarGrid(vlistID3, varID3);
	  field3.missval = vlistInqVarMissval(vlistID3, varID3);

	  farexpr(&field1, field2, field3, humidityIndex);
	  
	  pstreamDefRecord(streamID4, varID4, levelID1);
	  pstreamWriteRecord(streamID4, field1.ptr, (int)field1.nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);
  if ( field3.ptr ) Free(field3.ptr);

  cdoFinish();

  return 0;
}
