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

      Copy       copy            Copy datasets
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "par_io.h"
#include "pstream.h"


#ifdef __cplusplus
extern "C" {
#endif
int streamGrbInqDataScanningMode(void);
#if defined (__cplusplus)
}
#endif

void *Copy(void *argument)
{
  bool lconstvars = true;
  int streamID2 = CDI_UNDEFID;
  int vlistID2 = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;
  int nrecs;
  int varID, levelID;
  int nmiss;
  int ntsteps, nvars;
  double *array = NULL;
  par_io_t parIO;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
                cdoOperatorAdd("copy",   0, 0, NULL);
  int SELALL  = cdoOperatorAdd("selall", 0, 0, NULL);
  int SZIP    = cdoOperatorAdd("szip",   0, 0, NULL);
  // clang-format on

#ifdef HIRLAM_EXTENSIONS
  // KEEP in mind the difference between copy and selall with respect to unpacking and repacking the GRIB information!
  // Especially when setting the DataScanningMode.
  printf("cdo copy/selall : UNCHANGED_RECORD=%d\n",UNCHANGED_RECORD);
  //if (cdiGribDataScanningMode != -1) lcopy = false;
  printf("cdo copy/selall : cdiGribDataScanningMode=%d; lcopy=%d\n", streamGrbInqDataScanningMode(), lcopy);
#endif //#ifdef HIRLAM_EXTENSIONS

  int operatorID = cdoOperatorID();

  if ( operatorID == SZIP )
    {
      cdoCompType  = CDI_COMPRESS_SZIP;
      cdoCompLevel = 0;
    }

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  int tsID2 = 0;
  for ( int indf = 0; indf < nfiles; indf++ )
    {
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      int streamID1 = pstreamOpenRead(cdoStreamName(indf));

      int vlistID1 = pstreamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

	  vlistID2 = vlistDuplicate(vlistID1);
	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  ntsteps = vlistNtsteps(vlistID1);
	  nvars   = vlistNvars(vlistID1);

	  if ( ntsteps == 1 )
	    {
	      for ( varID = 0; varID < nvars; ++varID )
		if ( vlistInqVarTimetype(vlistID1, varID) != TIME_CONSTANT ) break;
	      
	      if ( varID == nvars ) ntsteps = 0;
	    }

	  if ( ntsteps == 0 && nfiles > 1 )
	    {	      
              lconstvars = false;
	      for ( varID = 0; varID < nvars; ++varID )
		vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
	    }

	  pstreamDefVlist(streamID2, vlistID2);

	  int gridsize = vlistGridsizeMax(vlistID1);
	  array = (double*) Malloc(gridsize*sizeof(double));
	  if ( cdoParIO )
	    {
	      fprintf(stderr, "Parallel reading enabled!\n");
	      parIO.array = (double*) Malloc(gridsize*sizeof(double));
	      parIO.array_size = gridsize;
	    }
	}
      else
	{
	  vlistCompare(vlistID1, vlistID2, CMP_ALL);
	}

      int tsID1 = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  pstreamDefTimestep(streamID2, tsID2);
	       
	  for ( int recID = 0; recID < nrecs; recID++ )
	    { 
	      if ( lcopy && (operatorID == SELALL || operatorID == SZIP) )
		{
		  pstreamInqRecord(streamID1, &varID, &levelID);

                  if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                    if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
                      continue;
                  
		  pstreamDefRecord(streamID2,  varID,  levelID);
		  pstreamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  if ( cdoParIO )
		    {
		      parIO.recID = recID; parIO.nrecs = nrecs;
		      /* fprintf(stderr, "in1 streamID %d varID %d levelID %d\n", streamID1, varID, levelID);*/
		      parReadRecord(streamID1, &varID, &levelID, array, &nmiss, &parIO);
		      /* fprintf(stderr, "in2 streamID %d varID %d levelID %d\n", streamID1, varID, levelID);*/
		    }
		  else
		    {
		      pstreamInqRecord(streamID1, &varID, &levelID);

                      if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                        if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
                          continue;

		      pstreamReadRecord(streamID1, array, &nmiss);
		    }
		  /*
		  if ( cdoParIO )
		    fprintf(stderr, "out1 %d %d %d\n", streamID2,  varID,  levelID);
		  */
		  pstreamDefRecord(streamID2,  varID,  levelID);
		  pstreamWriteRecord(streamID2, array, nmiss);
		  /*
		  if ( cdoParIO )
		    fprintf(stderr, "out2 %d %d %d\n", streamID2,  varID,  levelID);
		  */
		}
	    }

	  tsID1++;
	  tsID2++;
	}

      pstreamClose(streamID1);
    }

  pstreamClose(streamID2);

  if ( array ) Free(array);
  if ( vlistID2 != CDI_UNDEFID ) vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
