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

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>
#include "cdo_int.h"
#include "pstream.h"


static
void invertLevDes(int vlistID)
{
  int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; index++ )
    {
      int zaxisID1 = vlistZaxis(vlistID, index);
      int zaxisID2 = zaxisDuplicate(zaxisID1);
      int zaxistype = zaxisInqType(zaxisID1);

      int nlev = zaxisInqSize(zaxisID1);
      if ( nlev <= 1 ) continue;

      if ( zaxisInqLevels(zaxisID1, NULL) )
	{
          double *yv1 = (double*) Malloc(nlev*sizeof(double));
          double *yv2 = (double*) Malloc(nlev*sizeof(double));
	  zaxisInqLevels(zaxisID1, yv1);
	  for ( int ilev = 0; ilev < nlev; ++ilev ) yv2[nlev-ilev-1] = yv1[ilev];
	  zaxisDefLevels(zaxisID2, yv2);
          Free(yv1);
          Free(yv2);
	}

      if ( zaxisInqLbounds(zaxisID1, NULL) && zaxisInqUbounds(zaxisID1, NULL) )
	{
          double *yb1 = (double*) Malloc(nlev*sizeof(double));
          double *yb2 = (double*) Malloc(nlev*sizeof(double));
	  zaxisInqLbounds(zaxisID1, yb1);
	  for ( int ilev = 0; ilev < nlev; ++ilev ) yb2[nlev-ilev-1] = yb1[ilev];
	  zaxisDefLbounds(zaxisID2, yb2);

	  zaxisInqUbounds(zaxisID1, yb1);
	  for ( int ilev = 0; ilev < nlev; ++ilev ) yb2[nlev-ilev-1] = yb1[ilev];
	  zaxisDefUbounds(zaxisID2, yb2);
          Free(yb1);
          Free(yb2);
	}

      if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
        {
          int vctsize = zaxisInqVctSize(zaxisID1);		
          if ( vctsize && vctsize%2 == 0 )
            {
              double *vct1 = (double*) Malloc(vctsize*sizeof(double));
              double *vct2 = (double*) Malloc(vctsize*sizeof(double));
              zaxisInqVct(zaxisID1, vct1);
              for ( int i = 0; i < vctsize/2; ++i )
                {
                  vct2[vctsize/2-1-i] = vct1[i];
                  vct2[vctsize-1-i]   = vct1[vctsize/2+i];
                }
              zaxisDefVct(zaxisID2, vctsize, vct2);
              Free(vct1);
              Free(vct2);
            }
        }

      vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
    }
}


void *Invertlev(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int nlev, nlevel;
  int gridID, zaxisID, offset;
  bool linvert = false;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  cdoOperatorAdd("invertlev", func_all, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc == func_all || operfunc == func_hrd ) invertLevDes(vlistID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array = (double*) Malloc(gridsize*sizeof(double));

  int nvars = vlistNvars(vlistID1);

  double **vardata  = (double**) Malloc(nvars*sizeof(double*));
  int **varnmiss = (int**) Malloc(nvars*sizeof(int*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID    = vlistInqVarGrid(vlistID1, varID);
      zaxisID   = vlistInqVarZaxis(vlistID1, varID);
      gridsize  = gridInqSize(gridID);
      nlev      = zaxisInqSize(zaxisID);

      if ( nlev <= 1 )
	{
	  vardata[varID]  = NULL;
	  varnmiss[varID] = NULL;
	}
      else
	{
	  linvert = true;
	  vardata[varID]  = (double*) Malloc(gridsize*nlev*sizeof(double));
	  varnmiss[varID] = (int*) Malloc(nlev*sizeof(int));
	}
    }

  if ( linvert == false ) cdoWarning("No variables with invertable levels found!");

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( vardata[varID] )
	    {    
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      offset   = gridsize*levelID;

	      pstreamReadRecord(streamID1, vardata[varID]+offset, &nmiss);
	      varnmiss[varID][levelID] = nmiss;
	    }
	  else
	    {
	      pstreamDefRecord(streamID2, varID, levelID);
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

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vardata[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  pstreamDefRecord(streamID2, varID, levelID);

		  offset = gridsize*(nlevel-levelID-1);
		  nmiss = varnmiss[varID][nlevel-levelID-1];

		  pstreamWriteRecord(streamID2, vardata[varID]+offset, nmiss);
		}   
	    }
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array ) Free(array);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vardata[varID] )
	{
	  Free(varnmiss[varID]);
	  Free(vardata[varID]);
	}
    }

  Free(varnmiss);
  Free(vardata);

  cdoFinish();

  return 0;
}
