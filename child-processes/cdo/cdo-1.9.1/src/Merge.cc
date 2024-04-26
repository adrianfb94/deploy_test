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

      Merge      merge           Merge datasets with different fields
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


static
void checkDupEntry(int vlistID1, int vlistID2, const char *filename)
{
  char vname1[CDI_MAX_NAME], vname2[CDI_MAX_NAME];
  int k;
  int mlev1 = 0, mlev2 = 0;
  double *lev1 = NULL, *lev2 = NULL;

  int nvars1 = vlistNvars(vlistID1);
  int nvars2 = vlistNvars(vlistID2);

  for ( int varID1 = 0; varID1 < nvars1; ++varID1 )
    {
      vlistInqVarName(vlistID1, varID1, vname1);
      int param1   = vlistInqVarParam(vlistID1, varID1);
      int gridID1  = vlistInqVarGrid(vlistID1, varID1);
      int zaxisID1 = vlistInqVarZaxis(vlistID1, varID1);
      int gtype1   = gridInqType(gridID1);
      int gsize1   = gridInqSize(gridID1);
      int ztype1   = zaxisInqType(zaxisID1);
      int nlev1    = zaxisInqSize(zaxisID1);
      if ( nlev1 > mlev1 )
	{
	  mlev1 = nlev1;
	  lev1 = (double*) Realloc(lev1, mlev1*sizeof(double));
	}
      cdoZaxisInqLevels(zaxisID1, lev1);

      for ( int varID2 = 0; varID2 < nvars2; ++varID2 )
	{
	  vlistInqVarName(vlistID2, varID2, vname2);
	  int param2   = vlistInqVarParam(vlistID2, varID2);
	  int gridID2  = vlistInqVarGrid(vlistID2, varID2);
	  int zaxisID2 = vlistInqVarZaxis(vlistID2, varID2);
	  int gtype2   = gridInqType(gridID2);
	  int gsize2   = gridInqSize(gridID2);
	  int ztype2   = zaxisInqType(zaxisID2);
	  int nlev2    = zaxisInqSize(zaxisID2);
	  if ( gtype1 == gtype2 && gsize1 == gsize2 && ztype1 == ztype2 && nlev1 == nlev2 )
	    {
	      if ( nlev2 > mlev2 )
		{
		  mlev2 = nlev2;
		  lev2 = (double*) Realloc(lev2, mlev2*sizeof(double));
		}
	      cdoZaxisInqLevels(zaxisID2, lev2);

              if ( zaxisInqLevels(zaxisID1, NULL) && zaxisInqLevels(zaxisID2, NULL) )
                {
                  for ( k = 0; k < nlev2; ++k )
                    if ( !IS_EQUAL(lev1[k], lev2[k]) ) break;

                  if ( k == nlev2 )
                    {
                      if ( param1 < 0 || param2 < 0 )
                        {
                          if ( strcmp(vname1, vname2) == 0 )
                            cdoWarning("Duplicate entry of parameter %s in %s!", vname2, filename);
                        }
                      else
                        {
                          if ( param1 == param2 )
                            {
                              char paramstr[32];
                              cdiParamToString(param2, paramstr, sizeof(paramstr));
                              cdoWarning("Duplicate entry of parameter %s in %s!", paramstr, filename);
                            }
                        }
                    }
		}
	    }
	}
    }

  if ( lev1 ) Free(lev1);
  if ( lev2 ) Free(lev2);
}
/*
static
int vlistConstVars(int vlistID)
{
  int nvars = vlistNvars(vlistID);

  for ( int varID = 0; varID < nvars; ++varID )
    if ( vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT ) return 0;

  return 1;
}
*/

void *Merge(void *argument)
{
  int streamID1 = -1;
  int varID, varID2;
  int nrecs = 0;
  int levelID, levelID2;
  int index;
  int gridsize;
  int nmiss;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  int streamCnt = cdoStreamCnt();
  int nmerge    = streamCnt - 1;

  const char *ofilename = cdoStreamName(streamCnt-1)->args;

  if ( !cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename) )
    cdoAbort("Outputfile %s already exists!", ofilename);

  int *streamIDs = (int*) Malloc(nmerge*sizeof(int));
  int *vlistIDs  = (int*) Malloc(nmerge*sizeof(int));
  int *numrecs   = (int*) Malloc(nmerge*sizeof(int));
  int *numsteps  = (int*) Malloc(nmerge*sizeof(int));

  for ( index = 0; index < nmerge; index++ )
    {
      streamIDs[index] = pstreamOpenRead(cdoStreamName(index));
      vlistIDs[index]  = pstreamInqVlist(streamIDs[index]);
    }

  int taxisindex = 0;
  int vlistID1 = vlistIDs[0];
  int taxisID1 = vlistInqTaxis(vlistID1);
  if ( vlistNtsteps(vlistID1) == 0 )
    for ( index = 1; index < nmerge; index++ )
      {
        vlistID1 = vlistIDs[index];
        if ( vlistNtsteps(vlistID1) != 0 )
          {
            taxisindex = index;
            taxisID1 = vlistInqTaxis(vlistID1);
            break;
          }
      }

  int taxisID2 = taxisDuplicate(taxisID1);

  int vlistID2 = vlistCreate();
  vlistCopy(vlistID2, vlistIDs[0]);
  for ( index = 1; index < nmerge; index++ )
    {
      checkDupEntry(vlistID2, vlistIDs[index], cdoStreamName(index)->args);
      /* vlistCat(vlistID2, vlistIDs[index]); */
      vlistMerge(vlistID2, vlistIDs[index]);
    }

  for ( index = 0; index < nmerge; index++ ) numsteps[index] = -1;

  int numconst = 0;
  for ( index = 0; index < nmerge; index++ )
    {
      vlistID1  = vlistIDs[index];
      numsteps[index] = vlistNtsteps(vlistID1);
      if ( numsteps[index] == 0 ) numsteps[index] = 1;
      if ( numsteps[index] == 1 ) numconst++; 
    }

  if ( numconst > 0 && numconst < nmerge )
    for ( index = 0; index < nmerge; index++ )
      {
	if ( numsteps[index] == 1 ) 
	  {
	    vlistID1  = vlistIDs[index];
	    int nvars = vlistNvars(vlistID1);
	    for ( int varID = 0; varID < nvars; ++varID )
	      {
		varID2 = vlistMergedVar(vlistID1, varID);
		vlistDefVarTimetype(vlistID2, varID2, TIME_CONSTANT);
	      }
	  }
      }

  if ( cdoVerbose ) 
    {
      for ( index = 0; index < nmerge; index++ ) vlistPrint(vlistIDs[index]);
      vlistPrint(vlistID2);
    }
       
  int streamID2 = pstreamOpenWrite(cdoStreamName(streamCnt-1), cdoFiletype());

  vlistDefTaxis(vlistID2, taxisID2);
  pstreamDefVlist(streamID2, vlistID2);

  double *array = NULL;
  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID2);
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( tsID >= 0 )
    {
      for ( index = 0; index < nmerge; index++ )
	{
	  streamID1 = streamIDs[index];
	  vlistID1  = vlistIDs[index];
	  if ( vlistID1 == -1 ) continue;

	  numrecs[index] = pstreamInqTimestep(streamID1, tsID);
	}

      for ( index = 0; index < nmerge; index++ ) if ( numrecs[index] != 0 ) break;
      if ( index == nmerge ) break; // EOF on all input streams
      
      if ( tsID == 1 )
	{
	  for ( index = 0; index < nmerge; index++ )
	    if ( numrecs[index] == 0 && numsteps[index] == 1 ) vlistIDs[index] = -1;
	  /*
	  for ( index = 0; index < nmerge; index++ )
	    if ( vlistIDs[index] != -1 )
	      {
		firstindex = index;
		break;
	      }
	  */
	}
      /*
      for ( index = 0; index < nmerge; index++ )
	printf("tsID %d   %d sID %d vID %d nrecs %d\n", tsID, index, streamIDs[index], vlistIDs[index], numrecs[index]);
      */
      if ( numrecs[taxisindex] == 0 )
	{
	  for ( index = 1; index < nmerge; index++ )
	    if ( vlistIDs[index] != -1 && numrecs[index] != 0 )
	      cdoWarning("Input stream %d has %d timestep%s. Stream %d has more timesteps, skipped!", taxisindex+1, tsID, tsID==1?"":"s", index+1);
	  break;
	}
      else
	{
	  for ( index = 1; index < nmerge; index++ )
	    if ( vlistIDs[index] != -1 && numrecs[index] == 0 )
	      {
		cdoWarning("Input stream %d has %d timestep%s. Stream %d has more timesteps, skipped!", index+1, tsID, tsID==1?"":"s", taxisindex+1);
		break;
	      }
	  if ( index < nmerge ) break;
	}

      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( index = 0; index < nmerge; index++ )
	{
	  streamID1 = streamIDs[index];
	  vlistID1  = vlistIDs[index];
	  nrecs = numrecs[index];

	  if ( vlistID1 == -1 ) continue;

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);

	      varID2   = vlistMergedVar(vlistID1, varID);
	      levelID2 = vlistMergedLevel(vlistID1, varID, levelID);

	      if ( cdoVerbose )	cdoPrint("var %d %d %d %d", varID, levelID, varID2, levelID2);

	      pstreamDefRecord(streamID2, varID2, levelID2);
	      if ( lcopy )
		{
		  pstreamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  pstreamReadRecord(streamID1, array, &nmiss);
		  pstreamWriteRecord(streamID2, array, nmiss);
		}
	    }
	}

      tsID++;
   }

  for ( index = 0; index < nmerge; index++ )
    pstreamClose(streamIDs[index]);

  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( streamIDs ) Free(streamIDs);
  if ( vlistIDs  ) Free(vlistIDs);
  if ( numrecs   ) Free(numrecs);
  if ( numsteps  ) Free(numsteps);
 
  if ( ! lcopy )
    if ( array ) Free(array);

  cdoFinish();

  return 0;
}
