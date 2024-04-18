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

      Showinfo   showparam       Show parameters
      Showinfo   showcode        Show code numbers
      Showinfo   showname        Show variable names
      Showinfo   showstdname     Show variable standard names
      Showinfo   showlevel       Show levels
      Showinfo   showyear        Show years
      Showinfo   showmon         Show months
      Showinfo   showdate        Show dates
      Showinfo   showtime        Show timesteps
      Showinfo   showltype       Show level types
      Showinfo   showformat      Show file format
*/


#include <stdio.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "Showattribute.h"

void *Showinfo(void *argument)
{
  int date0 = 0;
  int year, month, day;
  int month0 = 0, year0 = 0;

  cdoInitialize(argument);

  // clang-format off
  int SHOWYEAR      = cdoOperatorAdd("showyear",      0, 0, NULL);
  int SHOWMON       = cdoOperatorAdd("showmon",       0, 0, NULL);
  int SHOWDATE      = cdoOperatorAdd("showdate",      0, 0, NULL);
  int SHOWTIME      = cdoOperatorAdd("showtime",      0, 0, NULL);
  int SHOWTIMESTAMP = cdoOperatorAdd("showtimestamp", 0, 0, NULL);
  int SHOWCODE      = cdoOperatorAdd("showcode",      0, 0, NULL);
  int SHOWUNIT      = cdoOperatorAdd("showunit",      0, 0, NULL);
  int SHOWPARAM     = cdoOperatorAdd("showparam",     0, 0, NULL);
  int SHOWNAME      = cdoOperatorAdd("showname",      0, 0, NULL);
  int SHOWSTDNAME   = cdoOperatorAdd("showstdname",   0, 0, NULL);
  int SHOWLEVEL     = cdoOperatorAdd("showlevel",     0, 0, NULL);
  int SHOWLTYPE     = cdoOperatorAdd("showltype",     0, 0, NULL);
  int SHOWFORMAT    = cdoOperatorAdd("showformat",    0, 0, NULL);
  int SHOWGRID      = cdoOperatorAdd("showgrid",      0, 0, NULL); 
  int SHOWATTS      = cdoOperatorAdd("showatts",      0, 0, NULL);
  int SHOWATTSGLOB  = cdoOperatorAdd("showattsglob",  0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID = pstreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);

  int nvars   = vlistNvars(vlistID);
  int taxisID = vlistInqTaxis(vlistID);
  int ntsteps = vlistNtsteps(vlistID);

  if ( operatorID == SHOWYEAR )
    {
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || year0 != year )
	      {
		year0 = year;
		fprintf(stdout, " %4d", year0);
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWMON )
    {
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || month0 != month )
	      {
		month0 = month;
		fprintf(stdout, " %2d", month0);
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWDATE )
    {
      char vdatestr[32];
      int tsID  = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);
	 
	    date2str(vdate, vdatestr, sizeof(vdatestr));

	    if ( tsID == 0 || date0 != vdate )
	      {
		date0 = vdate;
		fprintf(stdout, " %s", vdatestr);
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIME )
    {
      char vtimestr[32];
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vtime = taxisInqVtime(taxisID);

	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    fprintf(stdout, " %s", vtimestr);

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIMESTAMP )
    {
      char vdatetimestr[64];
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( pstreamInqTimestep(streamID, tsID) )
	  {
	    int vdate = taxisInqVdate(taxisID);
	    int vtime = taxisInqVtime(taxisID);

	    datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
	    fprintf(stdout, " %s", vdatetimestr);

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWCODE )
    {
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  fprintf(stdout, " %d", vlistInqVarCode(vlistID, varID));
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWGRID )  
    {
      fprintf(stdout, "# param nr | grid nr | z-axis nr:   /* Use in combination with operatores: griddes and zaxisdes */ \n");	  	  
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int gridID  = vlistInqVarGrid(vlistID, varID);
	  int zaxisID = vlistInqVarZaxis(vlistID, varID);

	  fprintf(stdout, "      %3d     %3d      %3d\n",
                  vlistInqVarCode(vlistID, varID), vlistGridIndex(vlistID, gridID) + 1, vlistZaxisIndex(vlistID, zaxisID) + 1);
	}
    }    
  else if ( operatorID == SHOWUNIT )
    {
      char varunits[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  varunits[0] = 0;
	  vlistInqVarUnits(vlistID, varID, varunits);
          if ( strlen(varunits) ) fprintf(stdout, " %s", varunits);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWPARAM )
    {
      char paramstr[32];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int param = vlistInqVarParam(vlistID, varID);
	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  fprintf(stdout, " %s", paramstr);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWNAME )
    {
      char varname[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID, varID, varname);
	  fprintf(stdout, " %s", varname);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWSTDNAME )
    {
      char stdname[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarStdname(vlistID, varID, stdname);
          fprintf(stdout, " %s", stdname[0] != 0 ? stdname : "unknown");
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWLEVEL )
    {
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int zaxisID = vlistInqVarZaxis(vlistID, varID);
	  int nlevs = zaxisInqSize(zaxisID);
	  for ( int levelID = 0; levelID < nlevs; levelID++ )
	    fprintf(stdout, " %.9g", cdoZaxisInqLevel(zaxisID, levelID));
	  fprintf(stdout, "\n");
	}
    }
  else if ( operatorID == SHOWLTYPE )
    {
      int nzaxis = vlistNzaxis(vlistID);
      for ( int index = 0; index < nzaxis; index++ )
	{
	  int zaxisID = vlistZaxis(vlistID, index);
	  int ltype = zaxis2ltype(zaxisID);

	  if ( ltype != -1 ) fprintf(stdout, " %d", ltype);
	}
      fprintf(stdout, "\n"); 
    }
  else if ( operatorID == SHOWFORMAT )
    {
      printFiletype(streamID, vlistID);
    }
  else if ( operatorID == SHOWSTDNAME )
    {
      char stdname[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarStdname(vlistID, varID, stdname);
          fprintf(stdout, " %s", stdname[0] != 0 ? stdname : "unknown");
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWATTS || operatorID == SHOWATTSGLOB )
    {
      if ( operatorID == SHOWATTS )
        {
          int vlistID = pstreamInqVlist(streamID);
          int nvars = vlistNvars(vlistID);
          for ( int varID = 0; varID < nvars; varID++ )
            {
              char varname[CDI_MAX_NAME];
              vlistInqVarName(vlistID, varID, varname);
              fprintf(stdout, "%s:\n", varname);

    	      int nattsvar;
    	      cdiInqNatts(vlistID, varID, &nattsvar);
    	      printAtts(vlistID, varID, nattsvar, NULL);
    	    }
        }
      fprintf(stdout, "Global:\n");
      int natts;
      cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
      printAtts(vlistID, CDI_GLOBAL, natts, NULL);
    }
  pstreamClose(streamID);

  cdoFinish();

  return 0;
}
