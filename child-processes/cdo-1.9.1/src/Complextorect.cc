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


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Complextorect(void *argument)
{
  int nrecs;
  int varID, levelID;
  int i;
  int datatype;
  int nmiss;

  cdoInitialize(argument);

  // int COMPLEXTORECT = cdoOperatorAdd("complextorect", 0, 0, NULL);

  // int operatorID = cdoOperatorID();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; ++varID )
    {
      datatype = vlistInqVarDatatype(vlistID2, varID);
      if ( datatype == CDI_DATATYPE_CPX64 )
	datatype = CDI_DATATYPE_FLT64;
      else
	datatype = CDI_DATATYPE_FLT32;

      vlistDefVarDatatype(vlistID2, varID, datatype);
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  int streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID3);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(2*gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));
  double *array3 = (double*) Malloc(gridsize*sizeof(double));
      
  int tsID  = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamDefRecord(streamID3, varID, levelID);
	      
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  pstreamReadRecord(streamID1, array1, &nmiss);

	  for ( i = 0; i < gridsize; ++i )
	    {
	      array2[i] = array1[2*i];
	      array3[i] = array1[2*i+1];
	    }

	  pstreamWriteRecord(streamID2, array2, nmiss);
	  pstreamWriteRecord(streamID3, array3, nmiss);
	}
       
      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);
  if ( array3 ) Free(array3);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  cdoFinish();

  return 0;
}
