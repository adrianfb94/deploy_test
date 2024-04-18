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

      Split      splitcode       Split codes
      Split      splitparam      Split parameters
      Split      splitname       Split variables
      Split      splitlevel      Split levels
      Split      splitgrid       Split grids
      Split      splitzaxis      Split zaxis
      Split      splittabnum     Split table numbers
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
void gen_filename(char *filename, bool swap_obase, const char *obase, const char *suffix)
{
  if ( swap_obase ) strcat(filename, obase);
  if ( suffix[0] ) strcat(filename, suffix);
}


void *Split(void *argument)
{
  int nchars = 0;
  int varID;
  int levelID, levID;
  int varID2, levelID2;
  int vlistID2;
  int *vlistIDs = NULL, *streamIDs = NULL;
  int  itmp[999];
  double ftmp[999];
  char filesuffix[32];
  char filename[8192];
  int nsplit = 0;
  int nmiss;
  bool swap_obase = false;
  const char *uuid_attribute = NULL;

  cdoInitialize(argument);

  if ( processSelf().m_ID != 0 ) cdoAbort("This operator can't be combined with other operators!");

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int SPLITCODE   = cdoOperatorAdd("splitcode",   0, 0, NULL);
  int SPLITPARAM  = cdoOperatorAdd("splitparam",  0, 0, NULL);
  int SPLITNAME   = cdoOperatorAdd("splitname",   0, 0, NULL);
  int SPLITLEVEL  = cdoOperatorAdd("splitlevel",  0, 0, NULL);
  int SPLITGRID   = cdoOperatorAdd("splitgrid",   0, 0, NULL);
  int SPLITZAXIS  = cdoOperatorAdd("splitzaxis",  0, 0, NULL);
  int SPLITTABNUM = cdoOperatorAdd("splittabnum", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  for( int i = 0; i < operatorArgc(); ++i )
    {
      if ( strcmp("swap", operatorArgv()[i]) == 0 )
          swap_obase = true;
      else if ( strncmp("uuid=", operatorArgv()[i], 5 ) == 0 )
          uuid_attribute = operatorArgv()[i] + 5;
      else cdoAbort("Unknown parameter: >%s<", operatorArgv()[0]); 
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  int nvars  = vlistNvars(vlistID1);

  if ( !swap_obase )
    {
      strcpy(filename, cdoStreamName(1)->args);
      nchars = strlen(filename);
    }

  const char *refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);
  
  if ( operatorID == SPLITCODE )
    {
      nsplit = 0;
      for ( varID = 0; varID < nvars; varID++ )
	{
	  int code = vlistInqVarCode(vlistID1, varID);
          int index;
	  for ( index = 0; index < varID; index++ )
	    if ( code == vlistInqVarCode(vlistID1, index) ) break;

	  if ( index == varID )
	    {
	      itmp[nsplit] = code;
	      nsplit++;
	    }
	}

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<int> codes(nsplit);
      for ( int index = 0; index < nsplit; ++index ) codes[index] = itmp[index];

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int code = vlistInqVarCode(vlistID1, varID);
	      if ( codes[index] == code )
		{
                  int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }

	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  if ( codes[index] > 9999 )
	    {
	      sprintf(filename+nchars, "%05d", codes[index]);
	      gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);
	    }
	  else if ( codes[index] > 999 )
	    {
	      sprintf(filename+nchars, "%04d", codes[index]);
	      gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);
	    }
	  else
	    {
	      sprintf(filename+nchars, "%03d", codes[index]);
	      gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);
	    }

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITPARAM )
    {
      char paramstr[32];
      nsplit = 0;
      for ( varID = 0; varID < nvars; varID++ )
	{
	  int param = vlistInqVarParam(vlistID1, varID);
          int index;
	  for ( index = 0; index < varID; index++ )
	    if ( param == vlistInqVarParam(vlistID1, index) ) break;

	  if ( index == varID )
	    {
	      itmp[nsplit] = param;
	      nsplit++;
	    }
	}

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<int> params(nsplit);
      for ( int index = 0; index < nsplit; ++index ) params[index] = itmp[index];

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int param   = vlistInqVarParam(vlistID1, varID);
	      if ( params[index] == param )
		{
                  int nlevs   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }

	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  cdiParamToString(params[index], paramstr, sizeof(paramstr));

	  filename[nchars] = '\0';
	  strcat(filename, paramstr);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITTABNUM )
    {
      nsplit = 0;
      for ( varID = 0; varID < nvars; varID++ )
	{
          int tabnum  = tableInqNum(vlistInqVarTable(vlistID1, varID));
          int index;
	  for ( index = 0; index < varID; index++ )
	    if ( tabnum == tableInqNum(vlistInqVarTable(vlistID1, index)) ) break;

	  if ( index == varID )
	    {
	      itmp[nsplit] = tabnum;
	      nsplit++;
	    }
	}

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<int> tabnums(nsplit);
      for ( int index = 0; index < nsplit; ++index ) tabnums[index] = itmp[index];

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int tabnum  = tableInqNum(vlistInqVarTable(vlistID1, varID));
	      if ( tabnums[index] == tabnum )
		{
                  int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%03d", tabnums[index]);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITNAME )
    {
      char varname[CDI_MAX_NAME];
      nsplit = nvars;

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));

      for ( int index = 0; index < nsplit; index++ )
	{
	  vlistClearFlag(vlistID1);
	  varID = index;
	  int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levID = 0; levID < nlevs; levID++ )
	    {
	      vlistDefIndex(vlistID1, varID, levID, index);
	      vlistDefFlag(vlistID1, varID, levID, TRUE);
	    }

	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  filename[nchars] = '\0';
	  vlistInqVarName(vlistID1, varID, varname);
	  strcat(filename, varname);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITLEVEL )
    {
      double level;
      int nzaxis = vlistNzaxis(vlistID1);
      nsplit = 0;
      for ( int index = 0; index < nzaxis; index++ )
	{
          int zaxisID = vlistZaxis(vlistID1, index);
	  int nlevs = zaxisInqSize(zaxisID);
	  for ( levID = 0; levID < nlevs; levID++ )
	    {
	      level = cdoZaxisInqLevel(zaxisID, levID);
              int i;
	      for ( i = 0; i < nsplit; i++ )
		if ( IS_EQUAL(level, ftmp[i]) ) break;
	      if ( i == nsplit )
		ftmp[nsplit++] = level;
	    }
	}

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<double> levels(nsplit);
      for ( int index = 0; index < nsplit; ++index ) levels[index] = ftmp[index];

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
              int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs = zaxisInqSize(zaxisID);
	      for ( levID = 0; levID < nlevs; levID++ )
		{
                  level = cdoZaxisInqLevel(zaxisID, levID);
		  if ( IS_EQUAL(levels[index], level) )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%06g", levels[index]);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);
   
	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITGRID )
    {
      int gridID;

      nsplit = vlistNgrids(vlistID1);

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<int> gridIDs(nsplit);
      for ( int index = 0; index < nsplit; ++index )
	gridIDs[index] = vlistGrid(vlistID1, index);

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      if ( gridIDs[index] == gridID )
		{
                  int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%02d", vlistGridIndex(vlistID1, gridIDs[index])+1);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else if ( operatorID == SPLITZAXIS )
    {
      nsplit = vlistNzaxis(vlistID1);

      vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
      streamIDs = (int*) Malloc(nsplit*sizeof(int));
      std::vector<int> zaxisIDs(nsplit);
      for ( int index = 0; index < nsplit; ++index )
	zaxisIDs[index] = vlistZaxis(vlistID1, index);

      for ( int index = 0; index < nsplit; ++index )
	{
	  vlistClearFlag(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      if ( zaxisIDs[index] == zaxisID )
		{
                  int nlevs = zaxisInqSize(zaxisID);
                  for ( levID = 0; levID < nlevs; levID++ )
		    {
		      vlistDefIndex(vlistID1, varID, levID, index);
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	    }
	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID1);
	  vlistIDs[index] = vlistID2;

	  sprintf(filename+nchars, "%02d", vlistZaxisIndex(vlistID1, zaxisIDs[index])+1);
	  gen_filename(filename, swap_obase, cdoStreamName(1)->args, filesuffix);

	  argument_t *fileargument = file_argument_new(filename);
	  streamIDs[index] = pstreamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);
	}
    }
  else
    {
      cdoAbort("not implemented!");
    }

  for ( int  index = 0; index < nsplit; index++ )
    {
      if ( uuid_attribute ) cdo_def_tracking_id(vlistIDs[index], uuid_attribute);

      pstreamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  double *array = NULL;
  if ( ! lcopy )
    {
      int gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double *) Malloc(gridsize*sizeof(double));
    }

  int nrecs;
  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      for ( int index = 0; index < nsplit; index++ )
	pstreamDefTimestep(streamIDs[index], tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  int index = vlistInqIndex(vlistID1, varID, levelID);
	  vlistID2 = vlistIDs[index];
	  varID2   = vlistFindVar(vlistID2, varID);
	  levelID2 = vlistFindLevel(vlistID2, varID, levelID);
	  /*
	    printf("%d %d %d %d %d %d\n", index, vlistID2, varID, levelID, varID2, levelID2);
	  */
	  pstreamDefRecord(streamIDs[index], varID2, levelID2);
	  if ( lcopy )
	    {
	      pstreamCopyRecord(streamIDs[index], streamID1);
	    }
	  else
	    {
	      pstreamReadRecord(streamID1, array, &nmiss);
	      pstreamWriteRecord(streamIDs[index], array, nmiss);
	    }
	}

      tsID++;
    }

  pstreamClose(streamID1);

  for ( int index = 0; index < nsplit; index++ )
    {
      pstreamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }
 
  if ( ! lcopy )
    if ( array ) Free(array);

  if ( vlistIDs  ) Free(vlistIDs);
  if ( streamIDs ) Free(streamIDs);

  cdoFinish();

  return 0;
}
