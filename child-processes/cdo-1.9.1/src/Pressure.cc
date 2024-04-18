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

      Pressure    pressure_fl          Pressure on full hybrid levels
      Pressure    pressure_hl          Pressure on half hybrid levels
      Pressure    deltap               Difference of two half hybrid levels
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"
#include "listarray.h"
#include "stdnametable.h"


void *Pressure(void *argument)
{
  int mode;
  gribcode_t gribcodes = {};
  int nrecs;
  int i, k, offset;
  int varID, levelID;
  int zaxisIDp, zaxisIDh = -1;
  int nhlevf = 0, nhlevh = 0, nlevel = 0;
  int nvct = 0;
  int nmiss;
  int psID = -1, lnpsID = -1;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  double minval, maxval;
  double *pout = NULL;

  cdoInitialize(argument);

  // clang-format off
  int PRESSURE_FL = cdoOperatorAdd("pressure_fl", 0, 0, NULL);
  int PRESSURE_HL = cdoOperatorAdd("pressure_hl", 0, 0, NULL);
  int DELTAP      = cdoOperatorAdd("deltap",      0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  int gridsize = vlist_check_gridsize(vlistID1);

  int nhlev;
  double *vct = vlist_read_vct(vlistID1, &zaxisIDh, &nvct, &nhlev, &nhlevf, &nhlevh);

  bool l3Dvars = (zaxisIDh != -1 && gridsize > 0);
  if ( !l3Dvars ) cdoAbort("No 3D variable with hybrid sigma pressure coordinate found!");
    
  double *ps_prog    = (double*) Malloc(gridsize*sizeof(double));
  double *deltap     = (double*) Malloc(gridsize*nhlevf*sizeof(double));
  double *full_press = (double*) Malloc(gridsize*nhlevf*sizeof(double));
  double *half_press = (double*) Malloc(gridsize*nhlevh*sizeof(double));

  if ( operatorID == PRESSURE_FL || operatorID == DELTAP )
    {
      if ( cdoVerbose ) cdoPrint("Creating ZAXIS_HYBRID .. (nhlevf=%d)", nhlevf);
      zaxisIDp = zaxisCreate(ZAXIS_HYBRID, nhlevf);
    }
  else
    {
      if ( cdoVerbose ) cdoPrint("Creating ZAXIS_HYBRID_HALF .. (nhlevh=%d)", nhlevh);
      zaxisIDp = zaxisCreate(ZAXIS_HYBRID_HALF, nhlevh);
    }

  double *level = (double*) Malloc(nhlevh*sizeof(double));
  for ( int l = 0; l < nhlevh; l++ ) level[l] = l+1;
  zaxisDefLevels(zaxisIDp, level);
  Free(level);

  zaxisDefVct(zaxisIDp, 2*nhlevh, vct);

  int nvars = vlistNvars(vlistID1);

  bool useTable = false;
  for ( varID = 0; varID < nvars; varID++ )
    {
      int tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      if ( tableNum > 0 && tableNum != 255 )
	{
	  useTable = true;
	  break;
	}
    }

  if ( cdoVerbose && useTable ) cdoPrint("Use code tables!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      int zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      int nlevel   = zaxisInqSize(zaxisID);
      int instNum  = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      int tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      int code     = vlistInqVarCode(vlistID1, varID);
      int param    = vlistInqVarParam(vlistID1, varID);

      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( useTable )
	{
	  if ( tableNum == 2 )
	    {
	      mode = WMO_MODE;
	      wmo_gribcodes(&gribcodes);
	    }
	  else if ( tableNum == 128 )
	    {
	      mode = ECHAM_MODE;
	      echam_gribcodes(&gribcodes);
	    }
          //  KNMI: HIRLAM model version 7.2 uses tableNum=1    (LAMH_D11*)
          //  KNMI: HARMONIE model version 36 uses tableNum=1   (grib*)   (opreational NWP version)
          //  KNMI: HARMONIE model version 38 uses tableNum=253 (grib,grib_md) and tableNum=1 (grib_sfx) (research version)
	  else if ( tableNum == 1 || tableNum == 253 )
	    {
	      mode = HIRLAM_MODE;
	      hirlam_harmonie_gribcodes(&gribcodes);
	    }
	  else
	    mode = -1;
	}
      else
	{
	  mode = ECHAM_MODE;
	  echam_gribcodes(&gribcodes);
	}

      if ( cdoVerbose )
        {
	  vlistInqVarName(vlistID1, varID, varname);
	  cdoPrint("Mode = %d  Center = %d TableNum =%d Code = %d Param = %s Varname = %s varID = %d",
                   mode, instNum,tableNum,  code, paramstr, varname, varID);
        }

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  strtolower(varname);

	  /*                        ECHAM                            ECMWF       */
	  if      ( strcmp(varname, "geosp") == 0 || strcmp(varname, "z")    == 0 ) code = 129;
	  else if ( strcmp(varname, "st")    == 0 || strcmp(varname, "t")    == 0 ) code = 130;
	  else if ( strcmp(varname, "aps")   == 0 || strcmp(varname, "sp"  ) == 0 ) code = 134;
	  else if ( strcmp(varname, "ps")    == 0 )                                 code = 134;
	  else if ( strcmp(varname, "lsp")   == 0 || strcmp(varname, "lnsp") == 0 ) code = 152;
	  /* else if ( strcmp(varname, "geopoth") == 0 ) code = 156; */
	}

      if ( mode == ECHAM_MODE )
	{
	  if      ( code == gribcodes.ps   && nlevel == 1 ) psID    = varID;
	  else if ( code == gribcodes.lsp  && nlevel == 1 ) lnpsID  = varID;
	}
      else if ( mode == WMO_MODE )
	{
	  if ( code == gribcodes.ps        && nlevel == 1 ) psID    = varID;
	}
      else if ( mode == HIRLAM_MODE )
	{
	  if ( code == gribcodes.ps        && nlevel == 1 ) psID    = varID;
        }
    }

  int pvarID = lnpsID;
  if ( zaxisIDh != -1 && lnpsID != -1 )
    {
      int gridID = vlistInqVarGrid(vlistID1, lnpsID);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  lnpsID = -1;
	  cdoWarning("Spectral LOG(%s) not supported - using %s!", var_stdname(surface_air_pressure), var_stdname(surface_air_pressure));
	}
    }

  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      pvarID = psID;
      if ( psID == -1 )
	cdoAbort("%s not found!", var_stdname(surface_air_pressure));
    }

  int gridID = vlistInqVarGrid(vlistID1, pvarID);
  if ( gridInqType(gridID) == GRID_SPECTRAL )
    cdoAbort("%s on spectral representation not supported!", var_stdname(surface_air_pressure));

  double *pdata = (double*) Malloc(gridsize*sizeof(double));

  int vlistID2 = vlistCreate();
  varID = vlistDefVar(vlistID2, gridID, zaxisIDp, TIME_VARYING);
  vlistDefVarParam(vlistID2, varID, cdiEncodeParam(1, 255, 255));
  vlistDefVarName(vlistID2, varID, "pressure");
  vlistDefVarStdname(vlistID2, varID, "air_pressure");
  vlistDefVarUnits(vlistID2, varID, "Pa");

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( varID == pvarID )
	    {	  
	      pstreamReadRecord(streamID1, pdata, &nmiss);
	      if ( nmiss > 0 ) cdoAbort("Missing valus unsupported!");
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  if ( lnpsID != -1 )
	    for ( i = 0; i < gridsize; i++ ) ps_prog[i] = exp(pdata[i]);
	  else if ( psID != -1 )
	    memcpy(ps_prog, pdata, gridsize*sizeof(double));

	  /* check range of ps_prog */
	  minmaxval(gridsize, ps_prog, NULL, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);
	    
	  presh(full_press, half_press, vct, ps_prog, nhlevf, gridsize);
	}

      if ( operatorID == PRESSURE_FL )
	{
	  nlevel = nhlevf;
	  pout = full_press;
	}
      else if ( operatorID == DELTAP )
	{
	  nlevel = nhlevf;
	  for ( k = 0; k < nhlevf; ++k )
	    for ( i = 0; i < gridsize; ++i )
	      {
		deltap[k*gridsize+i] = half_press[(k+1)*gridsize+i] - half_press[k*gridsize+i];
	      }

	  pout = deltap;
	}
      else if ( operatorID == PRESSURE_HL )
	{
	  nlevel = nhlevh;
	  pout = half_press;
	}
	  
      varID = 0;
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  pstreamDefRecord(streamID2, varID, levelID);
	  offset = levelID*gridsize;
	  pstreamWriteRecord(streamID2, pout+offset, 0);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( pdata      ) Free(pdata);
  if ( ps_prog    ) Free(ps_prog);
  if ( deltap     ) Free(deltap);
  if ( full_press ) Free(full_press);
  if ( half_press ) Free(half_press);
  if ( vct        ) Free(vct);

  cdoFinish();

  return 0;
}
