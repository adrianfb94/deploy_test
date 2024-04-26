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

      Fldstat    fldrange        Field range (max-min)
      Fldstat    fldmin          Field minimum
      Fldstat    fldmax          Field maximum
      Fldstat    fldsum          Field sum
      Fldstat    fldmean         Field mean
      Fldstat    fldavg          Field average
      Fldstat    fldstd          Field standard deviation
      Fldstat    fldstd1         Field standard deviation [Normalize by (n-1)]
      Fldstat    fldvar          Field variance
      Fldstat    fldvar1         Field variance [Normalize by (n-1)]
      Fldstat    fldpctl         Field percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"

static
void print_location_LL(int operfunc, int vlistID, int varID, int levelID, int gridID, double sglval, double *fieldptr,
		       int vdate, int vtime)
{
  static bool showHeader = true;
  int code = vlistInqVarCode(vlistID, varID);

  int year, month, day, hour, minute, second;
  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if ( gridInqType(gridID) == GRID_GAUSSIAN ||
       gridInqType(gridID) == GRID_LONLAT )
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      double level = cdoZaxisInqLevel(zaxisID, levelID);
      int nlon  = gridInqXsize(gridID);
      int nlat  = gridInqYsize(gridID);
      for ( int j = 0; j < nlat; ++j )
        for ( int i = 0; i < nlon; ++i )
          {
            if ( DBL_IS_EQUAL(fieldptr[j*nlon+i], sglval) )
              {
                double xval = gridInqXval(gridID,i);
                double yval = gridInqYval(gridID,j);
                if ( showHeader )
                  {
                    if ( operfunc == func_min )
                      fprintf(stdout, "  Date     Time     Code  Level   Lon      Lat          Minval\n");
                    else
                      fprintf(stdout, "  Date     Time     Code  Level   Lon      Lat          Maxval\n");
                    
                    showHeader = false;
                  }
		  
                fprintf(stdout, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d %3d %7g %9.7g %9.7g %12.5g\n",
                        year, month, day, hour, minute, second, code, level, xval, yval, sglval);
              }
          }
    }
}


void *Fldstat(void *argument)
{
  int gridID2, lastgrid = -1;
  int index;
  int nrecs;
  int varID, levelID;
  int nmiss;
  double sglval;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("fldrange", func_range,  0, NULL);
  cdoOperatorAdd("fldmin",   func_min,    0, NULL);
  cdoOperatorAdd("fldmax",   func_max,    0, NULL);
  cdoOperatorAdd("fldsum",   func_sum,    0, NULL);
  cdoOperatorAdd("fldmean",  func_meanw,  1, NULL);
  cdoOperatorAdd("fldavg",   func_avgw,   1, NULL);
  cdoOperatorAdd("fldstd",   func_stdw,   1, NULL);
  cdoOperatorAdd("fldstd1",  func_std1w,  1, NULL);
  cdoOperatorAdd("fldvar",   func_varw,   1, NULL);
  cdoOperatorAdd("fldvar1",  func_var1w,  1, NULL);
  cdoOperatorAdd("fldpctl",  func_pctl,   0, NULL);
  // clang-format on

  int operatorID  = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  bool needWeights = cdoOperatorF2(operatorID) != 0;

  double pn = 0;
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2double(operatorArgv()[0]);
      percentile_check_number(pn);
    }

  bool useweights = true;

  if ( needWeights )
    {
      unsigned npar = operatorArgc();
      if ( npar > 0 )
	{
	  char **parnames = operatorArgv();

	  if ( cdoVerbose )
	    for ( unsigned i = 0; i < npar; i++ )
	      cdoPrint("key %u = %s", i+1, parnames[i]);

	  if ( strcmp(parnames[0], "noweights") == 0 ) useweights = false;
	  else cdoAbort("Parameter >%s< unsupported! Supported parameter are: noweights", parnames[0]);
	}
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  double slon = 0;
  double slat = 0;
  gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &slon);
  gridDefYvals(gridID2, &slat);

  int ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  field_type field;
  field_init(&field);

  int lim = vlistGridsizeMax(vlistID1);
  field.ptr    = (double*) Malloc(lim*sizeof(double));
  field.weight = NULL;
  if ( needWeights )
    {
      field.weight = (double*) Malloc(lim*sizeof(double));
      if ( !useweights )
	{
	  cdoPrint("Using constant grid cell area weights!");
	  for ( int i = 0; i < lim; ++i ) field.weight[i] = 1;
	}
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      /* Precompute date + time for later representation in verbose mode */
      int vdate = 0, vtime = 0;
      if ( cdoVerbose )
        {
          if ( operfunc == func_min || operfunc == func_max )
            {
              vdate = taxisInqVdate(taxisID1);
              vtime = taxisInqVtime(taxisID1);
            }
        }

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field.ptr, &nmiss);

          field.nmiss   = (size_t)nmiss;
          field.grid = vlistInqVarGrid(vlistID1, varID);
	  field.size = gridInqSize(field.grid);

	  if ( needWeights && field.grid != lastgrid )
	    {
	      lastgrid = field.grid;
	      field.weight[0] = 1;
	      if ( useweights && field.size > 1 )
		{
		  bool wstatus = gridWeights(field.grid, field.weight) != 0;
		  if ( wstatus && tsID == 0 && levelID == 0 )
		    {
		      char varname[CDI_MAX_NAME];
		      vlistInqVarName(vlistID1, varID, varname);
		      cdoWarning("Grid cell bounds not available, using constant grid cell area weights for variable %s!", varname);
		    }
		}
	    }

	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operfunc == func_pctl )
	    sglval = fldpctl(field, pn);
	  else  
	    sglval = fldfun(field, operfunc);

	  if ( cdoVerbose && (operfunc == func_min || operfunc == func_max) )
	    print_location_LL(operfunc, vlistID1, varID, levelID, field.grid, sglval, field.ptr, vdate, vtime);

	  if ( DBL_IS_EQUAL(sglval, field.missval) )
	    nmiss = 1;
	  else
	    nmiss = 0;

	  pstreamDefRecord(streamID2, varID,  levelID);
	  pstreamWriteRecord(streamID2, &sglval, nmiss);
	}
      tsID++;
    }


  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  if ( field.ptr )    Free(field.ptr);
  if ( field.weight ) Free(field.weight);

  cdoFinish();

  return 0;
}
