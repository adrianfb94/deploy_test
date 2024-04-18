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

      Vertwind    vertwind      Convert the vertical velocity to [m/s]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"


#define R  287.07  /* spezielle Gaskonstante fuer Luft */
#define G  9.80665 /* Erdbeschleunigung */

void *Vertwind(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nvct = 0;
  int nmiss;
  int tempID = -1, sqID = -1, psID = -1, omegaID = -1;
  char varname[CDI_MAX_NAME];
  double *vct = NULL;
  double *hpress = NULL, *ps_prog = NULL;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  vlist_check_gridsize(vlistID1);

  int temp_code  = 130;
  int sq_code    = 133;
  int ps_code    = 134;
  int omega_code = 135;

  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; ++varID )
    {
      int code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  strtolower(varname);

	  if      ( strcmp(varname, "st")    == 0 ) code = temp_code;
	  else if ( strcmp(varname, "sq")    == 0 ) code = sq_code;
	  else if ( strcmp(varname, "aps")   == 0 ) code = ps_code;
	  else if ( strcmp(varname, "omega") == 0 ) code = omega_code;
	}

      if      ( code == temp_code  ) tempID  = varID;
      else if ( code == sq_code    ) sqID    = varID;
      else if ( code == ps_code    ) psID    = varID;
      else if ( code == omega_code ) omegaID = varID;
    }

  if ( tempID == -1 || sqID == -1 || omegaID == -1 )
    {
      if ( tempID  == -1 ) cdoWarning("Temperature (code 130) not found!");
      if ( sqID    == -1 ) cdoWarning("Specific humidity (code 133) not found!");
      if ( omegaID == -1 ) cdoWarning("Vertical velocity (code 135) not found!");
      cdoAbort("Parameter not found!");
    }

  /* Get missing values */
  double missval_t   = vlistInqVarMissval(vlistID1, tempID);
  double missval_sq  = vlistInqVarMissval(vlistID1, sqID);
  double missval_wap = vlistInqVarMissval(vlistID1, omegaID);
  double missval_out = missval_wap;

  int gridID  = vlistInqVarGrid(vlistID1, omegaID);
  int zaxisID = vlistInqVarZaxis(vlistID1, omegaID);

  if ( psID == -1 && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
    cdoAbort("Surface pressure (code 134) not found!");

  int gridsize = gridInqSize(gridID);
  int nlevel = zaxisInqSize(zaxisID);
  double *level = (double*) Malloc(nlevel*sizeof(double));
  cdoZaxisInqLevels(zaxisID, level);

  double *temp    = (double*) Malloc(gridsize*nlevel*sizeof(double));
  double *sq      = (double*) Malloc(gridsize*nlevel*sizeof(double));
  double *omega   = (double*) Malloc(gridsize*nlevel*sizeof(double));
  double *wms     = (double*) Malloc(gridsize*nlevel*sizeof(double));
  double *fpress  = (double*) Malloc(gridsize*nlevel*sizeof(double));


  if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
    {
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  size_t offset = (size_t)levelID*gridsize;
	  for ( int i = 0; i < gridsize; ++i )
	    fpress[offset+i] = level[levelID];
	}
    }
  else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
    {
      ps_prog = (double*) Malloc(gridsize*sizeof(double));
      hpress  = (double*) Malloc(gridsize*(nlevel+1)*sizeof(double));
  
      nvct = zaxisInqVctSize(zaxisID);
      if ( nlevel == (nvct/2 - 1) )
	{
	  vct = (double*) Malloc(nvct*sizeof(double));
	  zaxisInqVct(zaxisID, vct);
	}
      else
	cdoAbort("Unsupported vertical coordinate table format!");
    }
  else
    cdoAbort("Unsupported Z-Axis type!");


  vlistClearFlag(vlistID1);
  for ( levelID = 0; levelID < nlevel; ++levelID )
    vlistDefFlag(vlistID1, omegaID, levelID, TRUE);

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);
  vlistDefVarCode(vlistID2, 0, 40);
  vlistDefVarName(vlistID2, 0, "W");
  vlistDefVarLongname(vlistID2, 0, "Vertical velocity");
  vlistDefVarUnits(vlistID2, 0, "m/s");
  vlistDefVarMissval(vlistID2, 0, missval_out);

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

	  size_t offset = (size_t)levelID*gridsize;

	  if      ( varID == tempID )
	    pstreamReadRecord(streamID1, temp+offset, &nmiss);
	  else if ( varID == sqID )
	    pstreamReadRecord(streamID1, sq+offset, &nmiss);
	  else if ( varID == omegaID )
	    pstreamReadRecord(streamID1, omega+offset, &nmiss);
	  else if ( varID == psID && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	    pstreamReadRecord(streamID1, ps_prog, &nmiss);
	}

      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	presh(fpress, hpress, vct, ps_prog, nlevel, gridsize);

      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  size_t offset = (size_t)levelID*gridsize;

	  for ( int i = 0; i < gridsize; ++i )
	    {
	      if ( DBL_IS_EQUAL(temp[offset+i],missval_t)    || 
		   DBL_IS_EQUAL(omega[offset+i],missval_wap) ||
		   DBL_IS_EQUAL(sq[offset+i],missval_sq) )
		{
		  wms[offset+i] = missval_out;
		}
	      else
		{
	          // Virtuelle Temperatur bringt die Feuchteabhaengigkeit hinein
	          double tv = temp[offset+i] * (1. + 0.608*sq[offset+i]);

	          // Die Dichte erhaelt man nun mit der Gasgleichung rho=p/(R*tv) Level in Pa!
	          double rho = fpress[offset+i] / (R*tv);
	          /*
		    Nun daraus die Vertikalgeschwindigkeit im m/s, indem man die Vertikalgeschwindigkeit
                    in Pa/s durch die Erdbeschleunigung und die Dichte teilt
	          */
	          wms[offset+i] = omega[offset+i]/(G*rho);
	        }
            }
	}

      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  size_t offset = (size_t)levelID*gridsize;

	  int nmiss_out = 0;
	  for ( int i = 0; i < gridsize; i++ )
            if ( DBL_IS_EQUAL(wms[offset+i],missval_out) )
	      nmiss_out++;

	  pstreamDefRecord(streamID2, 0, levelID);
	  pstreamWriteRecord(streamID2, wms+offset, nmiss_out);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);
 
  vlistDestroy(vlistID2);

  Free(temp);
  Free(sq);
  Free(omega);
  Free(wms);
  Free(fpress);

  if ( ps_prog ) Free(ps_prog);
  if ( hpress )  Free(hpress);
  if ( vct ) Free(vct);

  Free(level);

  cdoFinish();

  return 0;
}
