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

      Fldstat2    fldcor         Correlation of two fields
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


/* routine corr copied from PINGO */
/* correclation in space */
static
double correlation_s(const double * restrict in0, const double * restrict in1,
		     const double * restrict weight, double missval1, double missval2, long gridsize)
{
  double sum0, sum1, sum00, sum01, sum11, wsum0;
  sum0 = sum1 = sum00 = sum01 = sum11 = 0;
  wsum0 = 0;
	
  for ( long i = 0; i < gridsize; ++i )
    {
      if ( IS_NOT_EQUAL(weight[i], missval1) && IS_NOT_EQUAL(in0[i], missval1) && IS_NOT_EQUAL(in1[i], missval2) )
	    {
	      sum0  += weight[i] * in0[i];
	      sum1  += weight[i] * in1[i];
	      sum00 += weight[i] * in0[i] * in0[i];
	      sum01 += weight[i] * in0[i] * in1[i];
	      sum11 += weight[i] * in1[i] * in1[i];
	      wsum0 += weight[i];
	    }
    }

  double out = IS_NOT_EQUAL(wsum0, 0) ?
               DIVMN((sum01 * wsum0 - sum0 * sum1),
	            SQRTMN((sum00 * wsum0 - sum0 * sum0) *
	                 (sum11 * wsum0 - sum1 * sum1))) : missval1;

  return out;
}

/* covariance in space */
static
double covariance_s(const double * restrict in0, const double * restrict in1,
		    const double * restrict weight, double missval1, double missval2, long gridsize)
{
  double sum0, sum1, sum01, wsum0, wsum00;
  sum0 = sum1 = sum01 = 0;
  wsum0 = wsum00 = 0;

  for ( long i = 0; i < gridsize; ++i )
    {
      if ( IS_NOT_EQUAL(weight[i], missval1) && IS_NOT_EQUAL(in0[i], missval1) && IS_NOT_EQUAL(in1[i], missval2) )
	{
	  sum0   += weight[i] * in0[i];
	  sum1   += weight[i] * in1[i];
	  sum01  += weight[i] * in0[i] * in1[i];
	  wsum0  += weight[i];
	  wsum00 += weight[i] * weight[i];
	}
    }

  double out = IS_NOT_EQUAL(wsum0, 0) ?
               (sum01 * wsum0 - sum0 * sum1) / (wsum0 * wsum0) : missval1;

  return out;
}


void *Fldstat2(void *argument)
{
  int gridID, lastgridID = -1;
  int gridID3;
  int nrecs, nrecs2;
  int varID, levelID;
  int nmiss1, nmiss2;
  bool wstatus = false;
  bool needWeights = true;
  double sglval = 0;
  char varname[CDI_MAX_NAME];

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("fldcor",   func_cor,   0, NULL);
  cdoOperatorAdd("fldcovar", func_covar, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  double slon = 0;
  double slat = 0;
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  int ngrids = vlistNgrids(vlistID1);

  for ( int index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID3, index, gridID3);

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamDefVlist(streamID3, vlistID3);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));
  double *weight = needWeights ? (double*) Malloc(gridsize*sizeof(double)) : NULL;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      nrecs2 = pstreamInqTimestep(streamID2, tsID);

      if ( nrecs2 == 0 )
	{
	  cdoWarning("Input streams have different number of time steps!");
	  break;
	}

      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamInqRecord(streamID2, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss1);
	  pstreamReadRecord(streamID2, array2, &nmiss2);

	  gridID = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  if ( needWeights && gridID != lastgridID )
	    {
	      lastgridID = gridID;
	      wstatus = gridWeights(gridID, weight) != 0;
	    }
	  if ( wstatus && tsID == 0 && levelID == 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      cdoWarning("Using constant grid cell area weights for variable %s!", varname);
	    }

	  double missval1 = vlistInqVarMissval(vlistID1, varID);
	  double missval2 = vlistInqVarMissval(vlistID2, varID);

	  if ( operfunc == func_cor )
	    {
	      sglval = correlation_s(array1, array2, weight, missval1, missval2, gridsize);
	    }
	  else if ( operfunc == func_covar )
	    {
	      sglval = covariance_s(array1, array2, weight, missval1, missval2, gridsize);
	    }

          int nmiss3 = DBL_IS_EQUAL(sglval, missval1) ? 1 : 0;

	  pstreamDefRecord(streamID3, varID,  levelID);
	  pstreamWriteRecord(streamID3, &sglval, nmiss3);
	}

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);
  if ( weight ) Free(weight);

  cdoFinish();

  return 0;
}
