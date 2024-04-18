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

     Splityear  splityear       Split in years
     Splityear  splityearmon    Split in years and month
*/

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define MAX_YEARS 99999

void *Splityear(void *argument)
{
  int streamID2 = -1;
  int varID;
  int nrecs;
  int levelID;
  int day;
  int year1, year2;
  int mon1, mon2;
  int gridsize;
  int ic = 0;
  int cyear[MAX_YEARS];
  int nmiss;
  int gridID;
  int nlevel;
  char filesuffix[32];
  char filename[8192];
  double *array = NULL;
  field_type **vars = NULL;

  cdoInitialize(argument);

  if ( processSelf().m_ID != 0 ) cdoAbort("This operator can't be combined with other operators!");
  
  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int SPLITYEAR    = cdoOperatorAdd("splityear",     0, 10000, NULL);
  int SPLITYEARMON = cdoOperatorAdd("splityearmon",  0,   100, NULL);
  // clang-format on
  
  int operatorID = cdoOperatorID();
  int operintval = cdoOperatorF2(operatorID);

  memset(cyear, 0, MAX_YEARS*sizeof(int));

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  const char *refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  // if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) nconst++;

  if ( nconst )
    {
      vars = (field_type **) Malloc(nvars*sizeof(field_type *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      gridsize = gridInqSize(gridID);
		  
	      vars[varID] = (field_type*) Malloc(nlevel*sizeof(field_type));

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  field_init(&vars[varID][levelID]);
		  vars[varID][levelID].grid = gridID;
		  vars[varID][levelID].ptr  = (double*) Malloc(gridsize*sizeof(double));
		}
	    }
	}
    }

  int index1 = -INT_MAX;
  int index2;
  year1 = -1;
  mon1  = -1;
  int tsID  = 0;
  int tsID2 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      cdiDecodeDate(vdate, &year2, &mon2, &day);

      if ( operatorID == SPLITYEAR )
	{
	  if ( tsID == 0 || year1 != year2 || (year1 == year2 && mon1 > mon2) )
	    {
	      tsID2 = 0;

	      if ( year1 != year2 ) ic = 0;
	      else                  ic++;

	      if ( year2 >= 0 && year2 < MAX_YEARS )
		{
		  ic = cyear[year2];
		  cyear[year2]++;
		}

	      year1 = year2;

	      if ( streamID2 >= 0 ) pstreamClose(streamID2);

	      sprintf(filename+nchars, "%04d", year1);
	      if ( ic > 0 ) sprintf(filename+strlen(filename), "_%d", ic+1);
	      if ( filesuffix[0] )
		sprintf(filename+strlen(filename), "%s", filesuffix);
	  
	      if ( cdoVerbose ) cdoPrint("create file %s", filename);

	      argument_t *fileargument = file_argument_new(filename);
	      streamID2 = pstreamOpenWrite(fileargument, cdoFiletype());
	      file_argument_free(fileargument);

	      pstreamDefVlist(streamID2, vlistID2);
	    }
	  mon1 = mon2;
	}
      else if ( operatorID == SPLITYEARMON )
	{
	  index2 = (vdate/operintval);
	  
	  if ( tsID == 0 || index1 != index2 )
	    {
	      tsID2 = 0;

	      index1 = index2;

	      if ( streamID2 >= 0 ) pstreamClose(streamID2);

	      sprintf(filename+nchars, "%04d", index1);
	      //if ( ic > 0 ) sprintf(filename+strlen(filename), "_%d", ic+1);
	      if ( filesuffix[0] )
		sprintf(filename+strlen(filename), "%s", filesuffix);
	  
	      if ( cdoVerbose ) cdoPrint("create file %s", filename);

	      argument_t *fileargument = file_argument_new(filename);
	      streamID2 = pstreamOpenWrite(fileargument, cdoFiletype());
	      file_argument_free(fileargument);

	      pstreamDefVlist(streamID2, vlistID2);
	    }
	}
      
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID2);

      if ( tsID > 0 && tsID2 == 0 && nconst )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
		{
		  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      pstreamDefRecord(streamID2, varID, levelID);
		      nmiss = vars[varID][levelID].nmiss;
		      pstreamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);

	  if ( lcopy && !(tsID == 0 && nconst) )
	    {
	      pstreamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      pstreamReadRecord(streamID1, array, &nmiss);
	      pstreamWriteRecord(streamID2, array, nmiss);

	      if ( tsID == 0 && nconst )
		{
		  if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
		    {
		      gridID  = vlistInqVarGrid(vlistID1, varID);
		      gridsize = gridInqSize(gridID);
		      memcpy(vars[varID][levelID].ptr, array, gridsize*sizeof(double));
		      vars[varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}

      tsID2++;
      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);
 
  if ( array ) Free(array);

  if ( nconst )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTimetype(vlistID2, varID) == TIME_CONSTANT )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( vars[varID][levelID].ptr )
		  Free(vars[varID][levelID].ptr);

	      Free(vars[varID]);
	    }
	}

      if ( vars ) Free(vars);
    }

  cdoFinish();

  return 0;
}
