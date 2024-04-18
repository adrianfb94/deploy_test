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

      Setvals     setvals       Set list of old values to new values
      Setrtoc     setrtoc       Set range to new value
      Setrtoc2    setrtoc2      Set range to new value others to value2
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "listarray.h"


void *Replacevalues(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int nvals = 0;
  lista_t *flista = lista_new(FLT_LISTA);
  double *fltarr = NULL;
  double rmin = 0, rmax = 0;
  double newval = 0, newval2 = 0;

  cdoInitialize(argument);

  // clang-format off
  int SETVALS  = cdoOperatorAdd("setvals",  0, 0, "I1,O1,...,In,On");
  int SETRTOC  = cdoOperatorAdd("setrtoc",  0, 0, "range (min, max), value");
  int SETRTOC2 = cdoOperatorAdd("setrtoc2", 0, 0, "range (min, max), value1, value2");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SETVALS )
    {
      nvals = args2flt_lista(operatorArgc(), operatorArgv(), flista);
      if ( nvals < 2 ) cdoAbort("Too few arguments!");
      if ( nvals % 2 != 0 )  cdoAbort("Need pairs of arguments!");
      fltarr = (double *) lista_dataptr(flista);
      nvals = nvals / 2;
    }
  else if ( operatorID == SETRTOC )
    {
      operatorCheckArgc(3);
      rmin   = parameter2double(operatorArgv()[0]);
      rmax   = parameter2double(operatorArgv()[1]);
      newval = parameter2double(operatorArgv()[2]);
    }
  else if ( operatorID == SETRTOC2 )
    {
      operatorCheckArgc(4);
      rmin    = parameter2double(operatorArgv()[0]);
      rmax    = parameter2double(operatorArgv()[1]);
      newval  = parameter2double(operatorArgv()[2]);
      newval2 = parameter2double(operatorArgv()[3]);
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  double missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operatorID == SETVALS )
	    {
	      for ( int i = 0; i < gridsize; i++ )
		if ( !DBL_IS_EQUAL(array[i], missval) )
		  {
		    /* printf("\nelem %d val %f ",i,array[i]); */
		    for ( int j = 0; j < nvals; j++ )
		      {
			if ( DBL_IS_EQUAL(array[i], fltarr[j*2] ) )
			  {
			    array[i] = fltarr[j*2+1];
			    /* printf("j=%d %f %f ",j,fltarr[j*2],fltarr[j*2+1]); */
			    break;
			  }
		      }
		  }
	    }
	  else if ( operatorID == SETRTOC )
	    {
	      for ( int i = 0; i < gridsize; i++ )
		if ( !DBL_IS_EQUAL(array[i], missval) )
		  {
		    if ( array[i] >= rmin && array[i] <= rmax)
                      array[i] = newval;
		  }
	    }
	  else if ( operatorID == SETRTOC2 )
	    {
	      for ( int i = 0; i < gridsize; i++ )
		if ( !DBL_IS_EQUAL(array[i], missval) )
		  {
		    if ( array[i] >= rmin && array[i] <= rmax )
                      array[i] = newval;
		    else
                      array[i] = newval2;
		  }
	    }

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array ) Free(array);

  lista_destroy(flista);

  cdoFinish();

  return 0;
}
