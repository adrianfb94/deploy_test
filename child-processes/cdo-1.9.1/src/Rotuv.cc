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

      Rotuv      rotuvb          Backward rotation
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


static
void rot_uv_back(int gridID, double *us, double *vs)
{
  double xpole, ypole, angle;
  if ( gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL )
    gridInqParamRLL(gridID, &xpole, &ypole, &angle);

  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  double *xvals = (double*) Malloc(nlon*sizeof(double));
  double *yvals = (double*) Malloc(nlat*sizeof(double));
  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_degree(units, 1, &xpole, "xpole");
  grid_to_degree(units, nlon, xvals, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, 1, &ypole, "ypole");
  grid_to_degree(units, nlat, yvals, "grid center lat");

  double u, v;
  for ( long ilat = 0; ilat < nlat; ilat++ )
    for ( long ilon = 0; ilon < nlon; ilon++ )
      {
	long i = ilat*nlon + ilon;

        double xval = lamrot_to_lam(yvals[ilat], xvals[ilon], ypole, xpole, angle);
        double yval = phirot_to_phi(yvals[ilat], xvals[ilon], ypole, angle);

	usvs_to_uv(us[i], vs[i], yval, xval, ypole, xpole, &u, &v);

	us[i] = u;
	vs[i] = v;
      }

  Free(xvals);
  Free(yvals);
}

#define  MAXARG     16384

void *Rotuv(void *argument)
{
  int varID, levelID;
  int varID1, varID2, nlevel1, nlevel2;
  int gridsize;
  int code, gridID;
  int offset;
  int nlevel;
  int chcodes[MAXARG];
  char *chvars[MAXARG];
  char varname[CDI_MAX_NAME];
  char varname2[CDI_MAX_NAME];
  double *single, *usvar = NULL, *vsvar = NULL;

  cdoInitialize(argument);

  operatorInputArg("pairs of u and v in the rotated system");

  int nch = operatorArgc();
  if ( nch%2 ) cdoAbort("Odd number of input arguments!");

  bool lvar = false; // We have a list of codes
  int len = (int)strlen(operatorArgv()[0]);
  int ix = (operatorArgv()[0][0] == '-') ? 1 : 0;
  for ( int i = ix; i < len; ++i )
    if ( !isdigit(operatorArgv()[0][i]) )
      {
        lvar = true; // We have a list of variables
        break;
      }

  if ( lvar )
    {
      for ( int i = 0; i < nch; i++ )
	chvars[i] = operatorArgv()[i];
    }
  else
    {
      for ( int i = 0; i < nch; i++ )
	chcodes[i] = parameter2int(operatorArgv()[i]);
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID1);
  int nrecs = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecs*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecs*sizeof(int));

  int **varnmiss    = (int **) Malloc(nvars*sizeof(int *));
  double **vardata  = (double **) Malloc(nvars*sizeof(double *));

  bool lfound[MAXARG];
  for ( int i = 0; i < nch; i++ ) lfound[i] = false;

  if ( lvar )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID2, varID, varname);
	  for ( int i = 0; i < nch; i++ )
	    if ( strcmp(varname, chvars[i]) == 0 ) lfound[i] = true;
	}
      for ( int i = 0; i < nch; i++ )
	if ( ! lfound[i] ) cdoAbort("Variable %s not found!", chvars[i]);
    }
  else
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  code = vlistInqVarCode(vlistID2, varID);
	  for ( int i = 0; i < nch; i++ )
	    if ( code == chcodes[i] ) lfound[i] = true;
	}
      for ( int i = 0; i < nch; i++ )
	if ( ! lfound[i] ) cdoAbort("Code %d not found!", chcodes[i]);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      if ( ! (gridInqType(gridID) == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL) )
	cdoAbort("Only rotated lon/lat grids supported!");

      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      varnmiss[varID] = (int*) Malloc(nlevel*sizeof(int));
      vardata[varID]  = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

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

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));	  
	  offset  = gridsize*levelID;
	  single  = vardata[varID] + offset;
	  pstreamReadRecord(streamID1, single, &varnmiss[varID][levelID]);
	  if ( varnmiss[varID][levelID] )
	    cdoAbort("Missing values unsupported for this operator!");
	}

      for ( int i = 0; i < nch; i += 2 )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lvar )
		{
		  vlistInqVarName(vlistID2, varID, varname);
		  if ( strcmp(varname, chvars[i]) == 0 ) break;
		}
	      else
		{
		  code = vlistInqVarCode(vlistID2, varID);
		  if ( code == chcodes[i] ) break;
		}
	    }

	  if ( varID == nvars )
	    cdoAbort("u-wind not found!");
	  else
	    usvar = vardata[varID];

	  varID1 = varID;
	  
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lvar )
		{
		  vlistInqVarName(vlistID2, varID, varname);
		  if ( strcmp(varname, chvars[i+1]) == 0 ) break;
		}
	      else
		{
		  code = vlistInqVarCode(vlistID2, varID);
		  if ( code == chcodes[i+1] ) break;
		}
	    }

	  if ( varID == nvars )
	    cdoAbort("v-wind not found!");
	  else
	    vsvar = vardata[varID];

	  varID2 = varID;

	  if ( cdoVerbose )
            {
              if ( lvar )
                {
                  vlistInqVarName(vlistID2, varID1, varname);
                  vlistInqVarName(vlistID2, varID2, varname2);
                  cdoPrint("Using var %s [%s](u) and var %s [%s](v)",
                           varname, chvars[i],
                           varname2, chvars[i+1]);
                }
              else
                cdoPrint("Using code %d [%d](u) and code %d [%d](v)",
                         vlistInqVarCode(vlistID1, varID1), chcodes[i],
                         vlistInqVarCode(vlistID1, varID2), chcodes[i+1]);
            }

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  nlevel1  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));
	  nlevel2  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID2));

	  if ( nlevel1 != nlevel2 )
	    cdoAbort("u-wind and v-wind have different number of levels!");

	  for ( levelID = 0; levelID < nlevel1; levelID++ )
	    {
	      offset = gridsize*levelID;
	      rot_uv_back(gridID, usvar + offset, vsvar + offset);
	    }
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata[varID] + offset;

	  pstreamDefRecord(streamID2, varID,  levelID);
	  pstreamWriteRecord(streamID2, single, varnmiss[varID][levelID]);     
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(varnmiss[varID]);
      Free(vardata[varID]);
    }

  cdoFinish();

  return 0;
}
