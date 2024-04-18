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

      Ninfo      npar            Number of parameters
      Ninfo      nlevel          Number of levels
      Ninfo      nyear           Number of years
      Ninfo      nmon            Number of months
      Ninfo      ndate           Number of dates
      Ninfo      ntime           Number of timesteps
      Ninfo      ngridpoints     Number of gridpoints
      Ninfo      ngrids          Number of grids
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Ninfo(void *argument)
{
  enum {NYEAR, NMON, NDATE, NTIME, NPAR, NLEVEL, NGRIDPOINTS, NGRIDS};
  int varID;
  int date0 = 0;
  int day, mon0 = 0, mon, year0 = 0, year;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("nyear"       , NYEAR       , 0 , NULL);
  cdoOperatorAdd("nmon"        , NMON        , 0 , NULL);
  cdoOperatorAdd("ndate"       , NDATE       , 0 , NULL);
  cdoOperatorAdd("ntime"       , NTIME       , 0 , NULL);
  cdoOperatorAdd("npar"        , NPAR        , 0 , NULL);
  cdoOperatorAdd("nlevel"      , NLEVEL      , 0 , NULL);
  cdoOperatorAdd("ngridpoints" , NGRIDPOINTS , 0 , NULL);
  cdoOperatorAdd("ngrids"      , NGRIDS      , 0 , NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  int streamID = pstreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);

  int nvars   = vlistNvars(vlistID);
  int taxisID = vlistInqTaxis(vlistID);
  int ntsteps = vlistNtsteps(vlistID);
  int ngrids  = vlistNgrids(vlistID);

  switch ( operfunc )
    {
    case NYEAR:
      {
      int nyear = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);
	    cdiDecodeDate(vdate, &year, &mon, &day);
	 
	    if ( tsID == 0 || year0 != year )
	      {
		year0 = year;
		nyear++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "%d\n", nyear);
      break;
      }
    case NMON:
      {
      int nmon = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);
	    cdiDecodeDate(vdate, &year, &mon, &day);
	 
	    if ( tsID == 0 || mon0 != mon )
	      {
		mon0 = mon;
		nmon++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "%d\n", nmon);
      break;
      }
    case NDATE:
      {
      int ndate = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);
	    
	    if ( tsID == 0 || date0 != vdate )
	      {
		date0 = vdate;
		ndate++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "%d\n", ndate);
      break;
      }
    case NTIME:
      {
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) ) tsID++;
      fprintf(stdout, "%d\n", tsID);
      break;
      }
    case NPAR:
      fprintf(stdout, "%d\n", nvars);
      break;
    case NLEVEL:
      for ( varID = 0; varID < nvars; varID++ )
	{
	  int zaxisID = vlistInqVarZaxis(vlistID, varID);
	  int levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, "%d\n", levelsize);
	}
      break;
    case NGRIDPOINTS:
      for ( varID = 0; varID < nvars; varID++ )
	{
	  int gridID = vlistInqVarGrid(vlistID, varID);
	  int gridsize = gridInqSize(gridID);
	  fprintf(stdout, "%d\n", gridsize);
	}
      break;
    case NGRIDS:
      fprintf(stdout, "%d\n", ngrids);
      break;
    default:
      cdoAbort("operator not implemented!");
      break;
    }

  pstreamClose(streamID);

  cdoFinish();

  return 0;
}
