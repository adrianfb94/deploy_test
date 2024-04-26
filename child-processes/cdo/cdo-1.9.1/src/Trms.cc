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

*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void trms(field_type field1, field_type field2, double *dp, field_type *field3)
{
  int i, k, nlev, len, rnmiss = 0;
  int    zaxis    = field1.zaxis;
  int    grid1    = field1.grid;
  double *array1  = field1.ptr;
  int    grid2    = field2.grid;
  double *array2  = field2.ptr;
  double missval1 = field1.missval;
  double missval2 = field2.missval;
  double *w       = field1.weight;
  double rsum = 0, rsumw = 0, ravg = 0, wp;

  nlev   = zaxisInqSize(zaxis);
  len    = gridInqSize(grid1);
  if ( len != gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  for ( k = 0; k < nlev; k++ ) 
    for ( i = 0; i < len; i++ ) 
      {
	wp = w[i]*dp[k*len+i];
	rsum  = ADDMN(rsum, MULMN(wp, MULMN( SUBMN(array2[k*len+i], array1[k*len+i]),
				      SUBMN(array2[k*len+i], array1[k*len+i]))));
	rsumw = ADDMN(rsumw, wp);
      }

  ravg = SQRTMN( DIVMN(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}

void *Trms(void *argument)
{
  int gridID1, gridID3, lastgrid = -1;
  int code = 0, oldcode = 0;
  int zaxisID;
  int nrecs;
  int nmiss;
  int varID, levelID;
  int pcode = 152, pvarID = -1;
  long offset;
  size_t vctsize = 0;
  const double *va = NULL, *vb = NULL;
  double *single;
  double sglval;

  cdoInitialize(argument);

  bool needWeights = true;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);

  double slon = 0;
  double slat = 0;
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  vlistClearFlag(vlistID1);
  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarCode(vlistID1, varID) == pcode )
	pvarID = varID;
      else
	vlistDefFlag(vlistID1, varID, 0, TRUE);
    }

  if ( pvarID == -1 ) cdoAbort("pressure variable missing!");

  int vlistID3 = vlistCreate();
  cdoVlistCopyFlag(vlistID3, vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int ngrids = vlistNgrids(vlistID1);
  int index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  
  if ( needWeights &&
       gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  vlistChangeGridIndex(vlistID3, index, gridID3);
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");

  int nzaxis = vlistNzaxis(vlistID1);
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID1, index);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	{
	  vctsize = zaxisInqVctSize(zaxisID);
	  const double *vct = zaxisInqVctPtr(zaxisID);
	  va = vct;
	  vb = vct + vctsize/2;
	  /*
	  for ( i = 0; i < vctsize/2; i++ )
	    fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
	  for ( i = 0; i < vctsize/2-1; i++ )
	    fprintf(stdout, "%5d %25.17f %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i],
		    (va[i+1] + vb[i+1]*101300) - (va[i] + vb[i]*101300));
	  */

	  break;
	}
    }

  if ( vctsize == 0 ) cdoAbort("VCT missing!");

  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  double **vardata1 = (double**) Malloc(nvars*sizeof(double*));
  double **vardata2 = (double**) Malloc(nvars*sizeof(double*));

  int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, pvarID));
  int nlevel   = vctsize/2 - 1;
  double *dp = (double*) Malloc(gridsize*nlevel*sizeof(double));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

  field_type field1, field2, field3;
  field_init(&field1);
  field_init(&field2);
  field_init(&field3);

  int lim = vlistGridsizeMax(vlistID1);
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) Malloc(lim*sizeof(double));

  field2.weight = NULL;

  field3.ptr  = &sglval;
  field3.grid = gridID3;

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      nrecs = pstreamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata1[varID] + offset;

	  pstreamReadRecord(streamID1, single, &nmiss);
	  if ( nmiss ) cdoAbort("Missing values unsupported for this operator!");

	  pstreamInqRecord(streamID2, &varID, &levelID);

	  single   = vardata2[varID] + offset;
	  pstreamReadRecord(streamID2, single, &nmiss);
	  if ( nmiss ) cdoAbort("Missing values unsupported for this operator!");
	}

      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, pvarID));
      for ( int i = 0; i < gridsize; i++ )
	{
	  vardata1[pvarID][i] = exp(vardata1[pvarID][i]);
	  vardata2[pvarID][i] = exp(vardata2[pvarID][i]);
	}

      nlevel = vctsize/2 - 1;
      for ( int k = 0; k < nlevel; k++ )
	{
	  offset = gridsize*k;
	  for ( int i = 0; i < gridsize; i++ )
	    {
	      double dp1 = (va[k+1] + vb[k+1]*vardata1[pvarID][i]) - (va[k] + vb[k]*vardata1[pvarID][i]);
	      double dp2 = (va[k+1] + vb[k+1]*vardata2[pvarID][i]) - (va[k] + vb[k]*vardata2[pvarID][i]);

	      dp[offset+i] = 0.5 * (dp1 + dp2);
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  field1.ptr = vardata1[varID];
	  field2.ptr = vardata2[varID];

	  field1.zaxis   = vlistInqVarZaxis(vlistID1, varID);
	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
          bool wstatus = false;
	  if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      wstatus = gridWeights(field1.grid, field1.weight);
	    }
	  code = vlistInqVarCode(vlistID1, varID);
	  if ( wstatus != 0 && tsID == 0 && code != oldcode )
	    cdoWarning("Using constant area weights for code %d!", oldcode=code);

	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID2, varID);
	  field3.missval = vlistInqVarMissval(vlistID3, varID);

	  trms(field1, field2, dp, &field3);

	  pstreamDefRecord(streamID3, varID, 0);
	  pstreamWriteRecord(streamID3, &sglval, field3.nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID3);

  if ( field1.weight ) Free(field1.weight);

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(vardata1[varID]);
      Free(vardata2[varID]);
    }

  Free(vardata1);
  Free(vardata2);
  Free(dp);

  cdoFinish();

  return 0;
}
