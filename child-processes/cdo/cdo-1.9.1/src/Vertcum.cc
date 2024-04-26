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

      Vertcum    vertcum         Vertical cumulative
      Vertcum    vertcumhl       Vertical cumulative on hybrid sigma half level
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"


#define IS_SURFACE_LEVEL(zaxisID)  (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

static
void add_vars_mv(int gridsize, double missval, const double *restrict var1, const double *restrict var2, double *restrict var3)
{
  double missval1 = missval;
  double missval2 = missval;
  /*
  for ( int i = 0; i < gridsize; ++i )
    var3[i] = ADDMN(var2[i], var1[i]);
  */
  for ( int i = 0; i < gridsize; ++i )
    {
      var3[i] = var2[i]; 
      if ( !DBL_IS_EQUAL(var1[i], missval1) )
        {
          if ( !DBL_IS_EQUAL(var2[i], missval2) )
            var3[i] += var1[i];
          else
            var3[i] = var1[i];
        }      
    }
}


void *Vertcum(void *argument)
{
  int nrecs;
  int i, nlevshl = 0;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  // clang-format off
                  cdoOperatorAdd("vertcum",    0,  0, NULL);
  int VERTCUMHL = cdoOperatorAdd("vertcumhl",  0,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID1);
  int **varnmiss = (int**) Malloc(nvars*sizeof(int*));
  double ***vardata1 = (double***) Malloc(nvars*sizeof(double**));
  double ***vardata2 = (double***) Malloc(nvars*sizeof(double**));

  int zaxisIDhl = -1;

  if ( operatorID == VERTCUMHL )
    {
      double *vct = NULL;
      int lhybrid = FALSE;
      int nzaxis  = vlistNzaxis(vlistID1);
      for ( int i = 0; i < nzaxis; ++i )
        {
          int zaxisID = vlistZaxis(vlistID1, i);
          int nlevs   = zaxisInqSize(zaxisID);

          if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevs > 1 )
            {
              int nvct = zaxisInqVctSize(zaxisID);
              if ( nlevs == (nvct/2 - 1) )
                {
                  if ( lhybrid == FALSE )
                    {
                      lhybrid = TRUE;
                      nlevshl = nlevs+1;

                      double *vct = (double*) Malloc(nvct*sizeof(double));
                      zaxisInqVct(zaxisID, vct);

                      zaxisIDhl = zaxisCreate(ZAXIS_HYBRID_HALF, nlevshl);
                      double *levels = (double*) Malloc(nlevshl*sizeof(double));
                      for ( levelID = 0; levelID < nlevshl; ++levelID ) levels[levelID] = levelID+1;
                      zaxisDefLevels(zaxisIDhl, levels);
                      Free(levels);
                      zaxisDefVct(zaxisIDhl, nvct, vct);
                      vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
                    }
                  else if ( vct )
                    {
                      if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
                        vlistChangeZaxisIndex(vlistID2, i, zaxisIDhl);
                    }
                }
            }
        }
      if ( vct ) Free(vct);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      int nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      int nlevs2   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

      varnmiss[varID] = (int*) Malloc(nlevs*sizeof(int));
      vardata1[varID] = (double**) Malloc(nlevs*sizeof(double*));
      vardata2[varID] = (double**) Malloc(nlevs2*sizeof(double*));
      for ( levelID = 0; levelID < nlevs; ++levelID )
        vardata1[varID][levelID] = (double*) Malloc(gridsize*sizeof(double));
      for ( levelID = 0; levelID < nlevs2; ++levelID )
        vardata2[varID][levelID] = (double*) Malloc(gridsize*sizeof(double));
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, vardata1[varID][levelID], &varnmiss[varID][levelID]);
        }

      for ( varID = 0; varID < nvars; ++varID )
	{
          double missval  = vlistInqVarMissval(vlistID2, varID);
          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          int nlevs2   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

          if ( operatorID == VERTCUMHL && nlevs2 == nlevshl )
            {
              for ( i = 0; i < gridsize; ++i ) vardata2[varID][0][i] = 0;
            }
          else
            {
              for ( i = 0; i < gridsize; ++i ) vardata2[varID][0][i] = vardata1[varID][0][i];
            }

          for ( levelID = 1; levelID < nlevs2; ++levelID )
            {
              if ( operatorID == VERTCUMHL && nlevs2 == nlevshl )
                add_vars_mv(gridsize, missval, vardata1[varID][levelID-1], vardata2[varID][levelID-1], vardata2[varID][levelID]);
              else
                add_vars_mv(gridsize, missval, vardata1[varID][levelID], vardata2[varID][levelID-1], vardata2[varID][levelID]);
            }

          if ( operatorID == VERTCUMHL && nlevs2 == nlevshl )
            {
              double *var1 = vardata2[varID][nlevs2-1];
              for ( levelID = 0; levelID < nlevs2; ++levelID )
                {
                  double *var2 = vardata2[varID][levelID];
                  for ( i = 0; i < gridsize; ++i )
                    {
                      if ( IS_NOT_EQUAL(var1[i], 0) )
                        var2[i] /= var1[i];
                      else
                        var2[i] = 0;
                    }
                }
            }
	}

      for ( varID = 0; varID < nvars; ++varID )
	{
          double missval  = vlistInqVarMissval(vlistID2, varID);
          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          int nlevs2   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
          for ( levelID = 0; levelID < nlevs2; ++levelID )
	    {
              double *single = vardata2[varID][levelID];

              nmiss = 0;
              for ( int i = 0; i < gridsize; ++i )
                if ( DBL_IS_EQUAL(single[i], missval) ) nmiss++;

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, single, nmiss);
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; ++varID )
    {
      int nlevs  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      int nlevs2 = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
      for ( levelID = 0; levelID < nlevs; ++levelID ) Free(vardata1[varID][levelID]);
      for ( levelID = 0; levelID < nlevs2; ++levelID ) Free(vardata2[varID][levelID]);
      Free(vardata1[varID]);
      Free(vardata2[varID]);
      Free(varnmiss[varID]);
    }

  Free(vardata1);
  Free(vardata2);
  Free(varnmiss);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
