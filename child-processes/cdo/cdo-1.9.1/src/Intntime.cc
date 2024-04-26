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

      Intntime   intntime        Time interpolation
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


void *Intntime(void *argument)
{
  int nlevel;
  int varID, levelID;
  int vdate, vtime;
  int offset;
  double *single1, *single2;
  double *vardatap;

  cdoInitialize(argument);

  operatorInputArg("number of timesteps between 2 timesteps");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");

  int numts = parameter2int(operatorArgv()[0]);
  if ( numts < 2 ) cdoAbort("parameter must be greater than 1!");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int **nmiss1   = (int **) Malloc(nvars*sizeof(int *));
  int **nmiss2   = (int **) Malloc(nvars*sizeof(int *));
  double **vardata1 = (double **) Malloc(nvars*sizeof(double *));
  double **vardata2 = (double **) Malloc(nvars*sizeof(double *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      nmiss1[varID]   = (int*) Malloc(nlevel*sizeof(int));
      nmiss2[varID]   = (int*) Malloc(nlevel*sizeof(int));
      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int calendar = taxisInqCalendar(taxisID1);

  int tsID = 0;
  int tsIDo = 0;
  int nrecs = pstreamInqTimestep(streamID1, tsID++);
  int vdate1 = taxisInqVdate(taxisID1);
  int vtime1 = taxisInqVtime(taxisID1);
  juldate_t juldate1 = juldate_encode(calendar, vdate1, vtime1);

  taxisCopyTimestep(taxisID2, taxisID1);
  pstreamDefTimestep(streamID2, tsIDo++);
  for ( int recID = 0; recID < nrecs; recID++ )
    {
      pstreamInqRecord(streamID1, &varID, &levelID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      offset   = gridsize*levelID;
      single1  = vardata1[varID] + offset;
      pstreamReadRecord(streamID1, single1, &nmiss1[varID][levelID]);

      pstreamDefRecord(streamID2, varID, levelID);
      pstreamWriteRecord(streamID2, single1, nmiss1[varID][levelID]);
    }

  while ( (nrecs = pstreamInqTimestep(streamID1, tsID++)) )
    {
      int vdate2 = taxisInqVdate(taxisID1);
      int vtime2 = taxisInqVtime(taxisID1);
      juldate_t juldate2 = juldate_encode(calendar, vdate2, vtime2);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single2  = vardata2[varID] + offset;
	  pstreamReadRecord(streamID1, single2, &nmiss2[varID][levelID]);
	}

      for ( int it = 1; it < numts; it++ )
	{
	  double seconds = it * juldate_to_seconds(juldate_sub(juldate2, juldate1)) / numts;
	  juldate_t juldate = juldate_add_seconds((int)lround(seconds), juldate1);

	  juldate_decode(calendar, juldate, &vdate, &vtime);

	  if ( cdoVerbose )
	    {
	      char vdatestr[32], vtimestr[32];	  
	      /*
		cdoPrint("juldate1 %f", juldate_to_seconds(juldate1));
		cdoPrint("juldate  %f", juldate_to_seconds(juldate));
		cdoPrint("juldate2 %f", juldate_to_seconds(juldate2));
	      */
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      cdoPrint("%s %s", vdatestr, vtimestr);
	    }

	  taxisDefVdate(taxisID2, vdate);
	  taxisDefVtime(taxisID2, vtime);
	  pstreamDefTimestep(streamID2, tsIDo++);

	  double fac1 = juldate_to_seconds(juldate_sub(juldate2, juldate)) / 
	                juldate_to_seconds(juldate_sub(juldate2, juldate1));
	  double fac2 = juldate_to_seconds(juldate_sub(juldate,  juldate1)) / 
	                juldate_to_seconds(juldate_sub(juldate2, juldate1));

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      varID    = recVarID[recID];
	      levelID  = recLevelID[recID];
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single1  = vardata1[varID] + offset;
	      single2  = vardata2[varID] + offset;

	      int nmiss3 = 0;

	      if ( nmiss1[varID][levelID] > 0 || nmiss2[varID][levelID] > 0 )
		{
		  double missval1 = vlistInqVarMissval(vlistID1, varID);
		  double missval2 = vlistInqVarMissval(vlistID2, varID);

		  for ( int i = 0; i < gridsize; i++ )
		    {
		      if ( !DBL_IS_EQUAL(single1[i], missval1) &&
			   !DBL_IS_EQUAL(single2[i], missval2) )
			array[i] = single1[i]*fac1 + single2[i]*fac2;
		      else if ( DBL_IS_EQUAL(single1[i], missval1) &&
				!DBL_IS_EQUAL(single2[i], missval2) && fac2 >= 0.5 )
			array[i] = single2[i];
		      else if ( DBL_IS_EQUAL(single2[i], missval2) &&
				!DBL_IS_EQUAL(single1[i], missval1) && fac1 >= 0.5 )
			array[i] = single1[i];
		      else
			{
			  array[i] = missval1;
			  nmiss3++;
			}
		    }
		}
	      else
		{
		  for ( int i = 0; i < gridsize; i++ )
		    array[i] = single1[i]*fac1 + single2[i]*fac2;
		}

	      pstreamDefRecord(streamID2, varID, levelID);
	      pstreamWriteRecord(streamID2, array, nmiss3);
	    }
	}

      taxisDefVdate(taxisID2, vdate2);
      taxisDefVtime(taxisID2, vtime2);
      pstreamDefTimestep(streamID2, tsIDo++);
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  single2  = vardata2[varID] + offset;

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, single2, nmiss2[varID][levelID]);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  vardatap        = vardata1[varID];
	  vardata1[varID] = vardata2[varID];
	  vardata2[varID] = vardatap;
	}

      //vdate1 = vdate2;
      //vtime1 = vtime2;
      juldate1 = juldate2;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(nmiss1[varID]);
      Free(nmiss2[varID]);
      Free(vardata1[varID]);
      Free(vardata2[varID]);
    }

  Free(nmiss1);
  Free(nmiss2);
  Free(vardata1);
  Free(vardata2);

  if ( array )  Free(array);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
