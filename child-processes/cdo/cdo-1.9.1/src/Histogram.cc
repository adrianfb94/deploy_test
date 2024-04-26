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
#include "listarray.h"


void *Histogram(void *argument)
{
  int nrecs, varID, levelID;
  int nmiss;
  int offset;
  int nlevel, zaxisID;
  double missval;

  cdoInitialize(argument);

  // clang-format off
  int HISTCOUNT = cdoOperatorAdd("histcount", 0, 0, NULL);
  int HISTSUM   = cdoOperatorAdd("histsum",   0, 0, NULL);
  int HISTMEAN  = cdoOperatorAdd("histmean",  0, 0, NULL);
  int HISTFREQ  = cdoOperatorAdd("histfreq",  0, 0, NULL);
  // clang-format on

  UNUSED(HISTSUM);

  int operatorID = cdoOperatorID();

  operatorInputArg("bins");

  lista_t *flista = lista_new(FLT_LISTA);
  int nbins = args2flt_lista(operatorArgc(), operatorArgv(), flista) - 1;
  if ( nbins < 1 ) cdoAbort("Too few arguments!");
  double *fltarr = (double *) lista_dataptr(flista);

  if ( cdoVerbose )
    {
      printf("nbins = %d\n", nbins);
      for ( int i = 0; i < nbins; i++ )
	printf("flt %d = %g\n", i+1, fltarr[i]);
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);

  int vlistID2 = vlistDuplicate(vlistID1);

  /* create zaxis for output bins */
  int zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nbins);
  double *bins = (double*) Malloc(nbins*sizeof(double));
  /* for ( int i = 0; i < nbins; i++ ) bins[i] = (fltarr[i]+fltarr[i+1])/2; */
  for ( int i = 0; i < nbins; i++ ) bins[i] = fltarr[i];
  zaxisDefLevels(zaxisID2, bins);
  Free(bins);
  zaxisDefLbounds(zaxisID2, fltarr);
  zaxisDefUbounds(zaxisID2, fltarr+1);
  zaxisDefName(zaxisID2, "bin");
  zaxisDefLongname(zaxisID2, "histogram bins");
  zaxisDefUnits(zaxisID2, "level");

  /* check zaxis: only 2D fields allowed */
  int nzaxis = vlistNzaxis(vlistID1);
  for ( int index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID1, index);
      nlevel = zaxisInqSize(zaxisID);
      if ( nlevel > 1 )
	cdoAbort("Found 3D field with %d levels. Only 2D fields allowed!", nlevel);
      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID2);
  double **vardata   = (double **) Malloc(nvars*sizeof(double *));
  double **varcount  = (double **) Malloc(nvars*sizeof(double *));
  double **vartcount = (double **) Malloc(nvars*sizeof(double *));
  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
      vardata[varID]  = (double*) Malloc(nbins*gridsize*sizeof(double));
      varcount[varID] = (double*) Malloc(nbins*gridsize*sizeof(double));
      vartcount[varID] = (double*) Malloc(gridsize*sizeof(double));
      memset(vardata[varID], 0, nbins*gridsize*sizeof(double));
      memset(varcount[varID], 0, nbins*gridsize*sizeof(double));
      memset(vartcount[varID], 0, gridsize*sizeof(double));
    }

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID1 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array, &nmiss);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  nmiss=0;
	  for ( int i = 0; i < gridsize; i++ )
	    {
	      if ( !DBL_IS_EQUAL(array[i], missval) )
		{
		  vartcount[varID][i] += 1;
		  int index = 0;
		  while( index < nbins )
		    {
		      offset = gridsize*index;
		      if ( !DBL_IS_EQUAL(vardata[varID][offset+i], missval) &&
			   array[i] >= fltarr[index] && array[i] < fltarr[index+1] )
			{
			  vardata[varID][offset+i]  += array[i];
			  varcount[varID][offset+i] += 1;
			  break;
			}
		      index++;
		    }
		}
	      else { /* missing value */
		nmiss++;
	      }
	    }
	}
      tsID1++;
    }


  pstreamDefTimestep(streamID2, 0);

  for ( varID = 0; varID < nvars; varID++ )
    {
      missval = vlistInqVarMissval(vlistID2, varID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));

      /* fix mising values */
      
      for ( int index = 0; index < nbins; index++ )
	{
	  nmiss = 0;
	  offset = gridsize*index;

	  for ( int i = 0; i < gridsize; i++ )
	    {
	      if ( vartcount[varID][i] > 0 )
		{
		  if ( operatorID == HISTMEAN || operatorID == HISTFREQ )
		    {
		      if ( varcount[varID][offset+i] > 0 ) 
			{
			  if ( operatorID == HISTMEAN )
			    vardata[varID][offset+i] /= varcount[varID][offset+i];	    
			  else 
			    vardata[varID][offset+i] = varcount[varID][offset+i] / vartcount[varID][i];
			} 
		    }
		}
	      else
		{
		  nmiss++;
		  varcount[varID][offset+i] = missval;
		  vardata[varID][offset+i] = missval;
		}
	    }

	  pstreamDefRecord(streamID2,  varID,  index);

	  if ( operatorID == HISTCOUNT )
	    pstreamWriteRecord(streamID2, varcount[varID]+offset, nmiss);
	  else
	    pstreamWriteRecord(streamID2, vardata[varID]+offset, nmiss);
	}
    }
  
  pstreamClose(streamID1);
  pstreamClose(streamID2);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  Free(vardata[varID]);
	  Free(varcount[varID]);
	  Free(vartcount[varID]);
	}

      Free(vardata);
      Free(varcount);
      Free(vartcount);
    }

  if ( array ) Free(array);

  lista_destroy(flista);

  cdoFinish();

  return 0;
}
