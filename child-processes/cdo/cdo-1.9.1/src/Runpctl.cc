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

      Runpctl    runpctl         Running percentiles
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles.h"


void *Runpctl(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int varID;
  int levelID;
  int nmiss;

  cdoInitialize(argument);

  cdoOperatorAdd("runpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number, number of timesteps");
  operatorCheckArgc(2);
  double pn  = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);
  int ndates = parameter2int(operatorArgv()[1]);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  field_type ***vars1 = (field_type ***) Malloc((ndates+1)*sizeof(field_type **));
  double *array = (double*) Malloc(ndates*sizeof(double));
  
  for ( int its = 0; its < ndates; its++ )
    vars1[its] = field_malloc(vlistID1, FIELD_PTR);

  int tsID;
  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      int nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) cdoAbort("File has less than %d timesteps!", ndates);

      dtlist_taxisInqTimestep(dtlist, taxisID1, tsID);
        
      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if ( tsID == 0 )
            {
              recVarID[recID]   = varID;
              recLevelID[recID] = levelID;
            }
          
          pstreamReadRecord(streamID1, vars1[tsID][varID][levelID].ptr, &nmiss);
          vars1[tsID][varID][levelID].nmiss = nmiss;
        }
    }

  int otsID = 0;
  while ( TRUE )
    {
      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;
          
          int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          double missval  = vlistInqVarMissval(vlistID1, varID);
          int nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              nmiss = 0;  
              for ( int i = 0; i < gridsize; i++ )
                {
                  int j = 0;
                  for ( int inp = 0, j = 0; inp < ndates; inp++ )
                    {
                      double val = vars1[inp][varID][levelID].ptr[i];
                      if ( !DBL_IS_EQUAL(val, missval) ) array[j++] = val;
                    }
                  
                  if ( j > 0 )
                    {
                      vars1[0][varID][levelID].ptr[i] = percentile(array, j, pn);
                    }
                  else
                    {
                      vars1[0][varID][levelID].ptr[i] = missval;
                      nmiss++;
                    }
                }
              vars1[0][varID][levelID].nmiss = nmiss;  
            }
        }

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, ndates);
      pstreamDefTimestep(streamID2, otsID);

      for ( int recID = 0; recID < nrecords; recID++ )
        {
          varID    = recVarID[recID];
          levelID  = recLevelID[recID];

	  if ( otsID && vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
        }

      otsID++;

      dtlist_shift(dtlist);

      vars1[ndates] = vars1[0];
      for ( int inp = 0; inp < ndates; inp++ ) vars1[inp] = vars1[inp+1];

      int nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      dtlist_taxisInqTimestep(dtlist, taxisID1, ndates-1);

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
          vars1[ndates-1][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  for ( int its = 0; its < ndates; its++ )
    field_free(vars1[its], vlistID1);

  Free(vars1);
  Free(array);
  
  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  dtlist_delete(dtlist);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
