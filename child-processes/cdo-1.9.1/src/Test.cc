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


void *Test(void *argument)
{
  /*
  int streamID1, streamID2;
  */

  cdoInitialize(argument);

  /*
  streamID1 = pstreamOpenRead(cdoStreamName(0));
  streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamClose(streamID2);
  pstreamClose(streamID1);
  */
  cdoFinish();

  return 0;
}


void *Test2(void *argument)
{
  /*
  int streamID1, streamID2, streamID3;
  */

  cdoInitialize(argument);

  /*
  streamID1 = pstreamOpenRead(cdoStreamName(0));
  streamID2 = pstreamOpenRead(cdoStreamName(1));
  streamID3 = pstreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
  */
  cdoFinish();

  return 0;
}


void *Testdata(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  int tsID2 = 0;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  int vlistID2 = vlistDuplicate(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));
  float *fval = (float*) Malloc(gridsize*sizeof(float));
  int *ival = (int*) Malloc(gridsize*sizeof(int));
  unsigned char *cval = (unsigned char*) Malloc(gridsize*sizeof(unsigned char)*4);
  unsigned char *cval2 = (unsigned char*) Malloc(gridsize*sizeof(unsigned char)*4);

  FILE *fp = fopen("testdata", "w");

  int tsID1 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID2);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);
	  
	  pstreamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  for ( int i = 0; i < gridsize; ++i )
	    {
	      fval[i] = (float) array[i];

	      memcpy(&ival[i], &fval[i], 4);
	      memcpy(&cval[i*4], &fval[i], 4);

	      cval2[i+gridsize*0] = cval[i*4+0];
	      cval2[i+gridsize*1] = cval[i*4+1];
	      cval2[i+gridsize*2] = cval[i*4+2];
	      cval2[i+gridsize*3] = cval[i*4+3];

	      if ( tsID1 == 0 && recID == 0 )
	      printf("%4d %3d %3d %3d %3d %d %g\n",
		     i, (unsigned int)cval[4*i+0], (unsigned int)cval[4*i+1], (unsigned int)cval[4*i+2], (unsigned int)cval[4*i+3], ival[i], fval[i]);
	    }

	  pstreamWriteRecord(streamID2, array, nmiss);

	  fwrite(cval, 4, gridsize, fp);
	}

      tsID1++;
      tsID2++;
    }
  
  fclose(fp);
  pstreamClose(streamID1);
  pstreamClose(streamID2);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}
