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

      Diff       diff            Compare two datasets
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Diff(void *argument)
{
  bool lhead = true;
  int nrecs, nrecs2;
  int varID1, varID2;
  int levelID;
  int nmiss1, nmiss2;
  int ndrec = 0, nd2rec = 0, ngrec = 0;
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  double abslim = 0., abslim2 = 1.e-3, rellim = 1.0;

  cdoInitialize(argument);

  // clang-format off
  int DIFF  = cdoOperatorAdd("diff",  0, 0, NULL);
  int DIFFP = cdoOperatorAdd("diffp", 0, 0, NULL);
  int DIFFN = cdoOperatorAdd("diffn", 0, 0, NULL);
  int DIFFC = cdoOperatorAdd("diffc", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if ( operatorArgc() >= 1 ) abslim = parameter2double(operatorArgv()[0]);
  if ( abslim < -1.e33 || abslim > 1.e+33 ) cdoAbort("Abs. limit out of range!");
  if ( operatorArgc() == 2 ) rellim = parameter2double(operatorArgv()[1]);
  if ( rellim < -1.e33 || rellim > 1.e+33 ) cdoAbort("Rel. limit out of range!");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

  int indg = 0;
  int tsID = 0;
  int taxisID = vlistInqTaxis(vlistID1);
  while ( TRUE )
    {
      nrecs = pstreamInqTimestep(streamID1, tsID);
      if ( nrecs > 0 )
	{
	  date2str(taxisInqVdate(taxisID), vdatestr, sizeof(vdatestr));
	  time2str(taxisInqVtime(taxisID), vtimestr, sizeof(vtimestr));
	}

      nrecs2 = pstreamInqTimestep(streamID2, tsID);

      if ( nrecs == 0 || nrecs2 == 0 ) break;

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID1, &levelID);
	  pstreamInqRecord(streamID2, &varID2, &levelID);

	  indg += 1;

	  int param    = vlistInqVarParam(vlistID1, varID1);
	  int code     = vlistInqVarCode(vlistID1, varID1);
	  int gridID   = vlistInqVarGrid(vlistID1, varID1);
	  int zaxisID  = vlistInqVarZaxis(vlistID1, varID1);
	  int gridsize = gridInqSize(gridID);
	  double missval1 = vlistInqVarMissval(vlistID1, varID1);
	  double missval2 = vlistInqVarMissval(vlistID2, varID2);

	  //checkrel = gridInqType(gridID) != GRID_SPECTRAL;
          bool checkrel = true;

	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  pstreamReadRecord(streamID1, array1, &nmiss1);
	  pstreamReadRecord(streamID2, array2, &nmiss2);

	  int ndiff = 0;
	  bool dsgn = false;
          bool zero = false;
          double absm = 0.0;
	  double relm = 0.0;
          double absdiff;

	  for ( int i = 0; i < gridsize; i++ )
	    {
	      if ( (DBL_IS_NAN(array1[i]) && !DBL_IS_NAN(array2[i])) ||
		  (!DBL_IS_NAN(array1[i]) &&  DBL_IS_NAN(array2[i])) )
		{
		  ndiff++;
		  relm = 1.0;
		}
	      else if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) )
		{
		  absdiff = fabs(array1[i] - array2[i]);
		  if ( absdiff > 0. ) ndiff++;

		  absm = MAX(absm, absdiff);

		  if ( array1[i]*array2[i] < 0. )
		    dsgn = true;
		  else if ( IS_EQUAL(array1[i]*array2[i], 0.) )
		    zero = true;
		  else
		    relm = MAX(relm, absdiff / MAX(fabs(array1[i]), fabs(array2[i])));
		}
	      else if ( (DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2)) || 
		       (!DBL_IS_EQUAL(array1[i], missval1) &&  DBL_IS_EQUAL(array2[i], missval2)) )
		{
		  ndiff++;
		  relm = 1.0;
		}
	    }

	  if ( ! cdoSilentMode || cdoVerbose )
	    {
	      if ( absm > abslim || (checkrel && relm >= rellim) || cdoVerbose )
		{
		  if ( lhead )
		    {
		      lhead = false;

		      set_text_color(stdout, BRIGHT, BLACK);
		      fprintf(stdout, "               Date     Time   Level Gridsize    Miss ");
		      fprintf(stdout, "   Diff ");
		      fprintf(stdout, ": S Z  Max_Absdiff Max_Reldiff");

		      if ( operatorID == DIFFN )
			fprintf(stdout, " : Parameter name");
		      else if ( operatorID == DIFF || operatorID == DIFFP )
			fprintf(stdout, " : Parameter ID");
		      else if ( operatorID == DIFFC )
			fprintf(stdout, " : Code number");
		      reset_text_color(stdout);

		      fprintf(stdout, "\n");
		    }

		  if ( operatorID == DIFFN ) vlistInqVarName(vlistID1, varID1, varname);
		  
		  set_text_color(stdout, BRIGHT, BLACK);
		  fprintf(stdout, "%6d ", indg);
		  reset_text_color(stdout);
		  set_text_color(stdout, RESET, BLACK);
		  fprintf(stdout, ":");
		  reset_text_color(stdout);
		
		  set_text_color(stdout, RESET, MAGENTA);
		  fprintf(stdout, "%s %s ", vdatestr, vtimestr);
		  reset_text_color(stdout);
		  set_text_color(stdout, RESET, GREEN);
                  double level = cdoZaxisInqLevel(zaxisID, levelID);
		  fprintf(stdout, "%7g ", level);
		  fprintf(stdout, "%8d %7d ", gridsize, MAX(nmiss1, nmiss2));
		  fprintf(stdout, "%7d ", ndiff);
		  reset_text_color(stdout);
		
		  set_text_color(stdout, RESET, BLACK);
		  fprintf(stdout, ":");
		  reset_text_color(stdout);
		  fprintf(stdout, " %c %c ", dsgn ? 'T' : 'F', zero ? 'T' : 'F');
		  set_text_color(stdout, RESET, BLUE);
		  fprintf(stdout, "%#12.5g%#12.5g", absm, relm);
		  set_text_color(stdout, RESET, BLACK);
		  fprintf(stdout, " : ");
		  reset_text_color(stdout);

		  set_text_color(stdout, BRIGHT, GREEN);
		  if ( operatorID == DIFFN )
		    fprintf(stdout, "%-11s", varname);
		  else if ( operatorID == DIFF || operatorID == DIFFP )
		    fprintf(stdout, "%-11s", paramstr);
		  else if ( operatorID == DIFFC )
		    fprintf(stdout, "%4d", code);
		  reset_text_color(stdout);
		      
		  fprintf(stdout, "\n");
		}
	    }

	  ngrec++;
	  if ( absm > abslim  || (checkrel && relm >= rellim) ) ndrec++;
	  if ( absm > abslim2 || (checkrel && relm >= rellim) ) nd2rec++;
	}
      tsID++;
    }

  if ( ndrec > 0 )
    {
      set_text_color(stdout, BRIGHT, RED);
      fprintf(stdout, "  %d of %d records differ", ndrec, ngrec);
      reset_text_color(stdout);
      fprintf(stdout, "\n");

      if ( ndrec != nd2rec && abslim < abslim2 )
	fprintf(stdout, "  %d of %d records differ more than 0.001\n", nd2rec, ngrec);
      /*  fprintf(stdout, "  %d of %d records differ more then one thousandth\n", nprec, ngrec); */
    }

  if ( nrecs == 0 && nrecs2 > 0 )
    cdoWarning("stream2 has more time steps than stream1!");
  if ( nrecs > 0 && nrecs2 == 0 )
    cdoWarning("stream1 has more time steps than stream2!");

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);

  cdoFinish();

  return 0;
}
