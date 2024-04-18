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

      Replace    replace         Replace variables
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_VARS  1024

void *Replace(void *argument)
{
  int varID, varID1, varID2;
  int nrecs = 0;
  int levelID, levelID2;
  int nrecs2;
  int nchvars = 0;
  int idx;
  char varname1[CDI_MAX_NAME], varname2[CDI_MAX_NAME];
  int nmiss;
  int varlist1[MAX_VARS], varlist2[MAX_VARS];
  int **varlevel = NULL;
  int **varnmiss2 = NULL;
  double **vardata2 = NULL;
  double *parray;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID2 = pstreamInqVlist(streamID2);

  /* compare all variables in vlistID2 */

  int nvars1 = vlistNvars(vlistID1);
  int nvars2 = vlistNvars(vlistID2);

  for ( varID2 = 0; varID2 < nvars2; varID2++ )
    {
      vlistInqVarName(vlistID2, varID2, varname2);

      for ( varID1 = 0; varID1 < nvars1; varID1++ )
	{
	  vlistInqVarName(vlistID1, varID1, varname1);
	  if ( strcmp(varname1, varname2) == 0 ) break;
	}

      if ( varID1 < nvars1 )
	{
	  int gridsize1 = gridInqSize(vlistInqVarGrid(vlistID1, varID1));
	  int nlevel1   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

	  int gridsize2 = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  int nlevel2   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));

	  if ( gridsize1 != gridsize2 )
	    cdoAbort("Variables have different gridsize!");

	  if ( nlevel1 < nlevel2 )
	    cdoAbort("Variables have different number of levels!");

	  if ( cdoVerbose ) cdoPrint("Variable %s replaced.", varname1);

	  varlist1[nchvars] = varID1;
	  varlist2[nchvars] = varID2;
	  nchvars++;
	  if ( nchvars > MAX_VARS ) cdoAbort("Internal problem - too many variables!");
	}
      else
	{
	  cdoPrint("Variable %s not found!", varname2);
	}
    }

  if ( nchvars )
    {
      vardata2  = (double **) Malloc(nchvars*sizeof(double *));
      varnmiss2 = (int **) Malloc(nchvars*sizeof(int *));
      varlevel  = (int **) Malloc(nchvars*sizeof(int *));
      for ( idx = 0; idx < nchvars; idx++ )
	{
	  varID1 = varlist1[idx];
	  varID2 = varlist2[idx];
	  int nlevel1  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));
	  int nlevel2  = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
	  int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  vardata2[idx]  = (double*) Malloc(nlevel2*gridsize*sizeof(double));
	  varnmiss2[idx] = (int*) Malloc(nlevel2*sizeof(int));
	  varlevel[idx] = (int*) Malloc(nlevel1*sizeof(int));
	  /*
	  for ( levelID = 0; levelID < nlevel1; levelID++ )
	    varlevel[idx][levelID] = levelID;
	  */
	  if ( nlevel2 <= nlevel1 )
	    {
	      double *level1 = (double*) Malloc(nlevel1*sizeof(double));
	      double *level2 = (double*) Malloc(nlevel2*sizeof(double));
	      cdoZaxisInqLevels(vlistInqVarZaxis(vlistID1, varID1), level1);
	      cdoZaxisInqLevels(vlistInqVarZaxis(vlistID2, varID2), level2);

	      for ( levelID = 0; levelID < nlevel1; levelID++ )
		varlevel[idx][levelID] = -1;
	      
	      for ( int l2 = 0; l2 < nlevel2; l2++ )
		{
		  int l1;
		  for ( l1 = 0; l1 < nlevel1; l1++ )
		    if ( IS_EQUAL(level2[l2], level1[l1]) )
		      {
			varlevel[idx][l1] = l2;
			break;
		      }

		  if ( l1 == nlevel1 ) cdoWarning("Level %g not found!", level2[l2]);
		}

	      Free(level1);
	      Free(level2);
	    }
	}
    }

  int vlistID3 = vlistDuplicate(vlistID1);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistDefTaxis(vlistID3, taxisID3);
  pstreamDefVlist(streamID3, vlistID3);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int nts2 = vlistNtsteps(vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      if ( tsID == 0 || (nts2 != 0 && nts2 != 1) )
	{
	  nrecs2 = pstreamInqTimestep(streamID2, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");

	  for ( int recID = 0; recID < nrecs2; recID++ )
	    {
	      pstreamInqRecord(streamID2, &varID, &levelID);
	      
	      for ( idx = 0; idx < nchvars; idx++ )
		if ( varlist2[idx] == varID )
		  {
		    int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		    int offset   = gridsize*levelID;
		    parray = vardata2[idx]+offset;
		    pstreamReadRecord(streamID2, parray, &nmiss);
		    varnmiss2[idx][levelID] = nmiss;
		    break;
		  }
	    }
	}

      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  parray = array;

	  for ( idx = 0; idx < nchvars; idx++ )
	    if ( varlist1[idx] == varID )
	      {
		levelID2 = varlevel[idx][levelID];
		if ( levelID2 != -1 )
		  {
		    int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		    int offset   = gridsize*levelID2;
		    parray   = vardata2[idx]+offset;
		    nmiss    = varnmiss2[idx][levelID2];
		    break;
		  }
	      }

	  if ( idx == nchvars ) pstreamReadRecord(streamID1, parray, &nmiss);

	  pstreamDefRecord(streamID3, varID, levelID);
	  pstreamWriteRecord(streamID3, parray, nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
 
  if ( vardata2 )
    {
      for ( idx = 0; idx < nchvars; idx++ )
	{
	  Free(vardata2[idx]);
	  Free(varnmiss2[idx]);
	  Free(varlevel[idx]);
	}

      Free(vardata2);
      Free(varnmiss2);
      Free(varlevel);
    }

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}
