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
#include "grid.h"


static
const char *filetypestr(int filetype)
{
  switch ( filetype )
    {
    case CDI_FILETYPE_GRB:  return ("GRIB");            break;
    case CDI_FILETYPE_GRB2: return ("GRIB2");           break;
    case CDI_FILETYPE_NC:   return ("NetCDF");          break;
    case CDI_FILETYPE_NC2:  return ("NetCDF2");         break;
    case CDI_FILETYPE_NC4:  return ("NetCDF4");         break;
    case CDI_FILETYPE_NC4C: return ("NetCDF4 classic"); break;
    case CDI_FILETYPE_NC5:  return ("NetCDF5");         break;
    case CDI_FILETYPE_SRV:  return ("SERVICE");         break;
    case CDI_FILETYPE_EXT:  return ("EXTRA");           break;
    case CDI_FILETYPE_IEG:  return ("IEG");             break;
    default:            return ("");
    }
}

static
const char *datatypestr(int datatype)
{
  static char str[20];

  str[0] = 0;
  snprintf(str, sizeof(str), "%d bit packed", datatype);

  if      ( datatype == CDI_DATATYPE_PACK   ) return ("P0");
  else if ( datatype > 0 && datatype <= 32  ) return (str);
  else if ( datatype == CDI_DATATYPE_CPX32  ) return ("C32");
  else if ( datatype == CDI_DATATYPE_CPX64  ) return ("C64");
  else if ( datatype == CDI_DATATYPE_FLT32  ) return ("32 bit floats");
  else if ( datatype == CDI_DATATYPE_FLT64  ) return ("64 bit floats");
  else if ( datatype == CDI_DATATYPE_INT8   ) return ("I8");
  else if ( datatype == CDI_DATATYPE_INT16  ) return ("I16");
  else if ( datatype == CDI_DATATYPE_INT32  ) return ("I32");
  else if ( datatype == CDI_DATATYPE_UINT8  ) return ("U8");
  else if ( datatype == CDI_DATATYPE_UINT16 ) return ("U16");
  else if ( datatype == CDI_DATATYPE_UINT32 ) return ("U32");
  else                                        return ("");
}

static
void print_stat(const char *sinfo, int memtype, int datatype, int filetype, off_t nvalues, double data_size, double file_size, double tw)
{
  nvalues /= 1000000;
  data_size /= 1024.*1024.*1024.;

  double rout = 0;
  if ( tw > 0 ) rout = nvalues/tw;

  if ( memtype == MEMTYPE_FLOAT )
    cdoPrint("%s Wrote %.1f GB of 32 bit floats to %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype), filetypestr(filetype), rout);
  else
    cdoPrint("%s Wrote %.1f GB of 64 bit floats to %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype), filetypestr(filetype), rout);

  file_size /= 1024.*1024.*1024.;

  rout = 0;
  if ( tw > 0 ) rout = 1024*file_size/tw;

  cdoPrint("%s Wrote %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, file_size, tw, rout);
}


void *CDIwrite(void *argument)
{
  int memtype = CDO_Memtype;
  int nvars = 10, nlevs = 0, ntimesteps = 30;
  const char *defaultgrid = "global_.2";
  int tsID, varID, levelID;
  int i;
  int vlistID;
  int zaxisID, taxisID;
  int vdate, vtime;
  int filetype = -1, datatype = -1;
  int irun, nruns = 1;
  unsigned int seed = 1;
  char sinfo[64];
  off_t nvalues = 0;
  double file_size = 0, data_size = 0;
  double tw, tw0, t0, twsum = 0;
  double ***vars = NULL;
  float *farray = NULL;

  srand(seed);
  sinfo[0] = 0;

  cdoInitialize(argument);

  if ( cdoVerbose ) cdoPrint("parameter: <nruns, <grid, <nlevs, <ntimesteps, <nvars>>>>>");

  if ( operatorArgc() > 5 ) cdoAbort("Too many arguments!");

  const char *gridfile = defaultgrid;
  if ( operatorArgc() >= 1 ) nruns = parameter2int(operatorArgv()[0]);
  if ( operatorArgc() >= 2 ) gridfile = operatorArgv()[1];
  if ( operatorArgc() >= 3 ) nlevs = parameter2int(operatorArgv()[2]);
  if ( operatorArgc() >= 4 ) ntimesteps = parameter2int(operatorArgv()[3]);
  if ( operatorArgc() >= 5 ) nvars = parameter2int(operatorArgv()[4]);

  if ( nruns <    0 ) nruns = 0;
  if ( nruns > 9999 ) nruns = 9999;

  if ( nlevs <= 0  ) nlevs = 1;
  if ( nlevs > 255 ) nlevs = 255;
  if ( ntimesteps <= 0 ) ntimesteps = 1;
  if ( nvars <= 0 ) nvars = 1;

  int gridID = cdoDefineGrid(gridfile);
  int gridsize = gridInqSize(gridID);

  if ( nlevs == 1 )
    zaxisID  = zaxisCreate(ZAXIS_SURFACE, 1);
  else
    {
      double *levels = (double *) Malloc(nlevs*sizeof(double));
      for ( i = 0; i < nlevs; ++i ) levels[i] = 100*i; 
      zaxisID  = zaxisCreate(ZAXIS_HEIGHT, nlevs);
      zaxisDefLevels(zaxisID, levels);
      Free(levels);
    }

  if ( cdoVerbose )
    {
      cdoPrint("nruns      : %d", nruns);
      cdoPrint("gridsize   : %d", gridsize);
      cdoPrint("nlevs      : %d", nlevs);
      cdoPrint("ntimesteps : %d", ntimesteps);
      cdoPrint("nvars      : %d", nvars);
    }

  double *array = (double*) Malloc(gridsize*sizeof(double));
  double *xvals = (double*) Malloc(gridsize*sizeof(double));
  double *yvals = (double*) Malloc(gridsize*sizeof(double));

  int gridID2 = gridID;
  if ( gridInqType(gridID) == GRID_GME ) gridID2 = gridToUnstructured(gridID, 0);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    gridID2 = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID2, xvals);
  gridInqYvals(gridID2, yvals);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID2, units);
  grid_to_radian(units, gridsize, xvals, "grid center lon");
  gridInqYunits(gridID2, units);
  grid_to_radian(units, gridsize, yvals, "grid center lat");

  for ( i = 0; i < gridsize; i++ )
    array[i] = 2 - cos(acos(cos(xvals[i]) * cos(yvals[i]))/1.2);

  Free(xvals);
  Free(yvals);

  vars = (double ***) Malloc(nvars*sizeof(double **));
  for ( varID = 0; varID < nvars; varID++ )
    {
      vars[varID] = (double **) Malloc(nlevs*sizeof(double *));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  vars[varID][levelID] = (double*) Malloc(gridsize*sizeof(double));
	  for ( i = 0; i < gridsize; ++i )
	    vars[varID][levelID][i] = varID + array[i]*(levelID+1);
	}
    }

  if ( memtype == MEMTYPE_FLOAT ) farray = (float*) Malloc(gridsize*sizeof(float));

  vlistID = vlistCreate();

  for ( i = 0; i < nvars; ++i )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(varID+1, 255, 255));
      //    vlistDefVarName(vlistID, varID, );
    }

  taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  // vlistDefNtsteps(vlistID, 1);

  for ( irun = 0; irun < nruns; ++irun )
    {
      tw0 = timer_val(timer_write);
      data_size = 0;
      nvalues = 0;

      int streamID = pstreamOpenWrite(cdoStreamName(0), cdoFiletype());

      pstreamDefVlist(streamID, vlistID);

      filetype = pstreamInqFiletype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);
      if ( datatype == CDI_UNDEFID ) datatype = CDI_DATATYPE_FLT32;
	  
      int julday = date_to_julday(CALENDAR_PROLEPTIC, 19870101);

      t0 = timer_val(timer_write);

      for ( tsID = 0; tsID < ntimesteps; tsID++ )
	{
	  vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
	  vtime = 0;
	  taxisDefVdate(taxisID, vdate);
	  taxisDefVtime(taxisID, vtime);
	  pstreamDefTimestep(streamID, tsID);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      for ( levelID = 0; levelID < nlevs; levelID++ )
		{
		  nvalues += gridsize;
		  pstreamDefRecord(streamID, varID, levelID);
		  if ( memtype == MEMTYPE_FLOAT )
		    {
		      double *darray = vars[varID][levelID];
		      for ( i = 0; i < gridsize; ++i ) farray[i] = darray[i];
		      pstreamWriteRecordF(streamID, farray, 0);
		      data_size += gridsize*4;
		    }
		  else
		    {
		      pstreamWriteRecord(streamID, vars[varID][levelID], 0);
		      data_size += gridsize*8;
		    }
		}
	    }

	  if ( cdoVerbose )
	    {
	      tw = timer_val(timer_write) - t0;
	      t0 = timer_val(timer_write);
	      cdoPrint("Timestep %d: %.2f seconds", tsID+1, tw);
	    }
	}

      pstreamClose(streamID);

      tw = timer_val(timer_write) - tw0;
      twsum += tw;

      file_size = (double) fileSize(cdoStreamName(0)->args);

      if ( nruns > 1 ) snprintf(sinfo, sizeof(sinfo), "(run %d)", irun+1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, data_size, file_size, tw);
    }

  if ( nruns > 1 )
    print_stat("(mean)", memtype, datatype, filetype, nvalues, data_size, file_size, twsum/nruns);

  vlistDestroy(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {
      for ( levelID = 0; levelID < nlevs; levelID++ ) Free(vars[varID][levelID]);
      Free(vars[varID]);
    }
  Free(vars);

  Free(array);

  if ( farray ) Free(farray);

  cdoFinish();

  return 0;
}
