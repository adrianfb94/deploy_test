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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *EstFreq(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridsize, nmiss;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int ntsteps = vlistNtsteps(vlistID1);

  double *array = NULL;
  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  int fyear, lyear, fmonth, lmonth, lymonth, dummy;
  int step_per_year = 0, step_per_month, currentyear, currentmon;

  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID == 0 )
        {
          cdiDecodeDate(taxisInqVdate(taxisID1), &fyear, &fmonth, &dummy);
          currentyear = fyear;
          currentmon = fmonth;
        }
      else
        cdiDecodeDate(taxisInqVdate(taxisID1), &currentyear, &currentmon, &dummy);
      if ( currentyear == fyear )
        {
          lymonth = currentmon;
          step_per_year++;
        }
      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  if ( lcopy )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}

      tsID++;
    }

  char frequency[CDI_MAX_NAME];

  if ( ntsteps > 2 )
    {
      int reclast = streamInqTimestep(streamID2, ntsteps);    
      cdiDecodeDate(taxisInqVdate(taxisID2), &lyear, &lmonth, &dummy);
/* First, estimation by maximal number of time steps divided by covered years between last and first time step */
      if ( cdoVerbose )
        printf("Frequency is calculated by dividing the number of time steps '%d' included in the time axis by the covered years of the time axis\ncomputed by the difference of the year of the last time stamp '%d' and the year of the first time stamp '%d'.\n", ntsteps, lyear, fyear);
      double covered_years = lyear-fyear + 1.0;
      if ( DBL_IS_EQUAL(ntsteps / covered_years, 1.) )
        strcpy(frequency, "yr");
      else if ( DBL_IS_EQUAL(ntsteps / covered_years, 12.) )
        strcpy(frequency, "mon");
      else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365.) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 365.25) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 366.) )
        strcpy(frequency, "day");
      else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365.*4) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 365.25*4) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 366.*4) )
        strcpy(frequency, "6hr");
      else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365.*8) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 365.25*8) ||
                DBL_IS_EQUAL(ntsteps / covered_years, 366.*8) )
        strcpy(frequency, "3hr");
      else 
        {
          int covered_months = lmonth-fmonth+1;
          if ( cdoVerbose )
            printf("The fraction ntsteps / covered_years = '%f' is neither 1, 12, 365, 365.25, 366 nor a multiple of 365 which would correspond to frequencies yearly, monthly, daily or subdaily respectively.\n Next try:\n\nFrequency is calculated by dividing the number of time steps '%d' in year '%d' by the covered months in that year '%d'.\n", ntsteps/covered_years, step_per_year, fyear, covered_months);
          if ( step_per_year > 366*8 )
            cdoAbort("Step per year '%d' in year '%d' is bigger than 366*8 which corresponds to a frequency of sub-3hourly! This is not yet enabled.", step_per_year, fyear);
          else
            {
              if ( (double)step_per_year / (double)covered_months > 31*8 )
                cdoAbort("Frequency is sub-3hourly! Not yet enabled.");
              else if ( (double)step_per_year / (double)covered_months > 31*4 )
                strcpy(frequency, "3hr");
              else if ( (double)step_per_year / (double)covered_months > 31 )
                strcpy(frequency, "6hr");
              else if ( (double)step_per_year / (double)covered_months > 1 )
                strcpy(frequency, "day");
              else
                strcpy(frequency, "mon");
            }
        }
    }
  else
    cdoAbort("For %d found timesteps no frequency can be computed - at least 3 timesteps are required.", ntsteps);
  if ( cdoVerbose )
    printf("Your file indicates a frequency of '%s'.\n", frequency);
  cdiDefAttTxt(vlistID2, CDI_GLOBAL, "frequency", 3, frequency);

  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}
