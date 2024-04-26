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

#define  MAX_LINE_LEN  4096

void *Setrcaname(void *argument)
{
  int nrecs;
  int varID, levelID;
  char **rcsnames;
  char line[MAX_LINE_LEN];
  char sname[CDI_MAX_NAME], sdescription[CDI_MAX_NAME], sunits[CDI_MAX_NAME];
  int scode, sltype, slevel;
  int zaxisID, ltype, code, nlev;
  int level;
  int gridsize, nmiss;
  double *array = NULL;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  operatorInputArg("file name with RCA names");
  rcsnames = operatorArgv();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID2);

  FILE *fp = fopen(rcsnames[0], "r");
  if ( fp != NULL )
    {
      while ( readline(fp, line, MAX_LINE_LEN) )
	{
	  sscanf(line, "%d\t%d\t%d\t%s\t%s\t%s", &scode, &sltype, &slevel, sname, sdescription, sunits);
	  /*
	  printf("%s\n", line);
	  printf("%d:%d:%d:%s:%s:%s\n", scode, sltype, slevel, sname, sdescription, sunits);
	  */
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      code = vlistInqVarCode(vlistID2, varID);
	      zaxisID = vlistInqVarZaxis(vlistID2, varID);
	      nlev = zaxisInqSize(zaxisID);

	      ltype = zaxis2ltype(zaxisID);

	      if ( code == scode )
		{
		  if ( ltype == 105 )
		    {
		      if ( nlev != 1 )
			{
			  cdoWarning("Number of levels should be 1 for level type 105!");
			  cdoWarning("Maybe environment variable SPLIT_LTYPE_105 is not set.");
			  continue;
			}
		      level = (int) cdoZaxisInqLevel(zaxisID, 0);
		      if ( sltype == 105 && slevel == level )
			{
			  vlistDefVarName(vlistID2, varID, sname);
			  vlistDefVarLongname(vlistID2, varID, sdescription);
			  vlistDefVarUnits(vlistID2, varID, sunits);
			  break;
			}
		    }
		  else if ( sltype != 105 )
		    {
		      vlistDefVarName(vlistID2, varID, sname);
		      vlistDefVarLongname(vlistID2, varID, sdescription);
		      vlistDefVarUnits(vlistID2, varID, sunits);
		      break;
		    }
		}
	    }
	}

      fclose(fp);
    }
  else
    {
      perror(rcsnames[0]);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);

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

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) Free(array);

  cdoFinish();

  return 0;
}
