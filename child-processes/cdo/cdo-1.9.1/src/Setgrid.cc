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

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Setgrid(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridID2 = -1;
  int gridtype = -1;
  int nmiss;
  int areasize = 0;
  int  masksize = 0;
  bool lregular = false;
  bool lregularnn = false;
  bool ldereference = false;
  int number = 0, position = 0;
  int grid2_nvgp;
  int lbounds = TRUE;
  int *grid2_vgpm = NULL;
  char *gridname = NULL;
  char *griduri = NULL;
  double *gridmask = NULL;
  double *areaweight = NULL;

  cdoInitialize(argument);

  // clang-format off
  int SETGRID       = cdoOperatorAdd("setgrid",       0, 0, "grid description file or name");
  int SETGRIDTYPE   = cdoOperatorAdd("setgridtype",   0, 0, "grid type");
  int SETGRIDAREA   = cdoOperatorAdd("setgridarea",   0, 0, "filename with area weights");
  int SETGRIDMASK   = cdoOperatorAdd("setgridmask",   0, 0, "filename with grid mask");
  int UNSETGRIDMASK = cdoOperatorAdd("unsetgridmask", 0, 0, NULL);
  int SETGRIDNUMBER = cdoOperatorAdd("setgridnumber", 0, 0, "grid number and optionally grid position");
  int SETGRIDURI    = cdoOperatorAdd("setgriduri",    0, 0, "reference URI of the horizontal grid");
  int USEGRIDNUMBER = cdoOperatorAdd("usegridnumber", 0, 0, "use existing grid identified by grid number");
  // clang-format on

  int operatorID = cdoOperatorID();

  // open stream before calling cdoDefineGrid!!!
  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  if ( operatorID != UNSETGRIDMASK )
    operatorInputArg(cdoOperatorEnter(operatorID));  

  if ( operatorID == SETGRID )
    {
      operatorCheckArgc(1);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      operatorCheckArgc(1);
      gridname = operatorArgv()[0];

      if      ( strcmp(gridname, "curvilinear0") == 0 )  {gridtype = GRID_CURVILINEAR; lbounds = 0;}
      else if ( strcmp(gridname, "curvilinear") == 0 )   {gridtype = GRID_CURVILINEAR; lbounds = 1;}
      else if ( strcmp(gridname, "cell") == 0 )           gridtype = GRID_UNSTRUCTURED;
      else if ( strcmp(gridname, "unstructured0") == 0 ) {gridtype = GRID_UNSTRUCTURED; lbounds = 0;}
      else if ( strcmp(gridname, "unstructured") == 0 )  {gridtype = GRID_UNSTRUCTURED; lbounds = 1;}
      else if ( strcmp(gridname, "generic") == 0 )        gridtype = GRID_GENERIC;
      else if ( strcmp(gridname, "dereference") == 0 )    ldereference = true;
      else if ( strcmp(gridname, "lonlat") == 0 )         gridtype = GRID_LONLAT;
      else if ( strcmp(gridname, "gaussian") == 0 )       gridtype = GRID_GAUSSIAN;
      else if ( strcmp(gridname, "regularnn") == 0 )     {gridtype = GRID_GAUSSIAN; lregularnn = true;}
      else if ( strcmp(gridname, "regular") == 0 )       {gridtype = GRID_GAUSSIAN; lregular = true;}
      else cdoAbort("Unsupported grid name: %s", gridname);
    }
  else if ( operatorID == SETGRIDAREA )
    {
      operatorCheckArgc(1);
      char *areafile = operatorArgv()[0];

      argument_t *fileargument = file_argument_new(areafile);
      int streamID = pstreamOpenRead(fileargument);
      file_argument_free(fileargument);

      int vlistID = pstreamInqVlist(streamID);

      nrecs = pstreamInqTimestep(streamID, 0);
      pstreamInqRecord(streamID, &varID, &levelID);

      int gridID = vlistInqVarGrid(vlistID, varID);
      areasize = gridInqSize(gridID);
      areaweight = (double*) Malloc(areasize*sizeof(double));
  
      pstreamReadRecord(streamID, areaweight, &nmiss);
      pstreamClose(streamID);

      if ( cdoVerbose )
	{
	  double arrmean = areaweight[0];
	  double arrmin  = areaweight[0];
	  double arrmax  = areaweight[0];
	  for ( int i = 1; i < areasize; i++ )
	    {
	      if ( areaweight[i] < arrmin ) arrmin = areaweight[i];
	      if ( areaweight[i] > arrmax ) arrmax = areaweight[i];
	      arrmean += areaweight[i];
	    }
	  arrmean = arrmean/areasize;

	  cdoPrint("areaweights: %d %#12.5g%#12.5g%#12.5g", areasize, arrmin, arrmean, arrmax);
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      operatorCheckArgc(1);
      char *maskfile = operatorArgv()[0];
      argument_t *fileargument = file_argument_new(maskfile);
      int streamID = pstreamOpenRead(fileargument);
      file_argument_free(fileargument);

      int vlistID = pstreamInqVlist(streamID);

      nrecs = pstreamInqTimestep(streamID, 0);
      pstreamInqRecord(streamID, &varID, &levelID);

      double missval  = vlistInqVarMissval(vlistID, varID);
      int gridID   = vlistInqVarGrid(vlistID, varID);
      masksize = gridInqSize(gridID);
      gridmask = (double*) Malloc(masksize*sizeof(double));
  
      pstreamReadRecord(streamID, gridmask, &nmiss);
      pstreamClose(streamID);

      for ( int i = 0; i < masksize; i++ )
	if ( DBL_IS_EQUAL(gridmask[i], missval) ) gridmask[i] = 0;
    }
  else if ( operatorID == USEGRIDNUMBER )
    {
      operatorCheckArgc(1);
      number = parameter2int(operatorArgv()[0]);
    }
  else if ( operatorID == SETGRIDNUMBER )
    {
      if ( operatorArgc() >= 1 && operatorArgc() <= 2 )
	{
	  number = parameter2int(operatorArgv()[0]);
	  if ( operatorArgc() == 2 ) position = parameter2int(operatorArgv()[1]);
	}
      else
	{
	  operatorCheckArgc(1);
	}
    }
  else if ( operatorID == SETGRIDURI )
    {
      operatorCheckArgc(1);
      griduri = operatorArgv()[0];
    }

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETGRID )
    {
      int found = 0;
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No grid with %d points found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDNUMBER || operatorID == SETGRIDURI || operatorID == USEGRIDNUMBER )
    {
      if ( operatorID == SETGRIDNUMBER )
	{
          int gridID1 = vlistGrid(vlistID1, 0);
	  gridID2 = gridCreate(GRID_UNSTRUCTURED, gridInqSize(gridID1));
	  gridDefNumber(gridID2, number);
	  gridDefPosition(gridID2, position);
	}
      else if ( operatorID == USEGRIDNUMBER ) 
        {
	  if ( number < 1 || number > vlistNgrids(vlistID1) )
	    cdoAbort("Invalid grid number: %d (max = %d)!", number, vlistNgrids(vlistID1));
          
	  gridID2 = vlistGrid(vlistID1, number-1);
        }
      else
	{
          int gridID1 = vlistGrid(vlistID1, 0);
	  gridID2 = gridDuplicate(gridID1);
	  gridDefReference(gridID2, griduri);
	}

      int found = 0;
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No horizontal grid with %d cells found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      bool lrgrid = false;
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1 = vlistGrid(vlistID1, index);
	  gridID2 = -1;

	  if ( gridInqType(gridID1) == GRID_GENERIC && gridInqSize(gridID1) == 1 ) continue;
	  
	  if ( lregular || lregularnn )
	    {
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
                gridID2 = gridToRegular(gridID1);
	    }
	  else if ( ldereference )
	    {
	      gridID2 = referenceToGrid(gridID1);
	      if ( gridID2 == -1 ) cdoAbort("Reference to horizontal grid not found!");
	    }
	  else
	    {
	      if      ( gridtype == GRID_CURVILINEAR  )
		{
		  gridID2 = gridToCurvilinear(gridID1, lbounds);
		}
	      else if ( gridtype == GRID_UNSTRUCTURED )
		{
                  bool ligme = false;
		  if ( gridInqType(gridID1) == GRID_GME ) ligme = true;
		  gridID2 = gridToUnstructured(gridID1, 1);

		  if ( ligme )
		    {
		      grid2_nvgp = gridInqSize(gridID2);
		      grid2_vgpm = (int*) Malloc(grid2_nvgp*sizeof(int));
		      gridInqMaskGME(gridID2, grid2_vgpm);
		      gridCompress(gridID2);
		    }
		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_CURVILINEAR )
		{
		  gridID2 = gridCurvilinearToRegular(gridID1);
		  if ( gridID2 == -1 ) cdoWarning("Conversion of curvilinear grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_UNSTRUCTURED )
		{
		  gridID2 = -1;
		  if ( gridID2 == -1 ) cdoWarning("Conversion of unstructured grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_GENERIC )
		{
		  gridID2 = -1;
		  if ( gridID2 == -1 ) cdoWarning("Conversion of generic grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_LONLAT )
		{
		  gridID2 = gridID1;
		}
	      else cdoAbort("Unsupported grid name: %s", gridname);
	    }

	  if ( gridID2 == -1 )
            {
              if ( !(lregular || lregularnn) )
                cdoAbort("Unsupported grid type!");
            }

	  if ( gridID2 != -1 )
            {
              if ( lregular || lregularnn ) lrgrid = true;
              vlistChangeGridIndex(vlistID2, index, gridID2);
            }
        }

      if ( (lregular || lregularnn) && !lrgrid )
        cdoWarning("No reduced Gaussian grid found!");
    }
  else if ( operatorID == SETGRIDAREA )
    {
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1  = vlistGrid(vlistID1, index);
	  int gridsize = gridInqSize(gridID1);
	  if ( gridsize == areasize )
	    {
	      gridID2 = gridDuplicate(gridID1);
	      gridDefArea(gridID2, areaweight);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	    }
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1  = vlistGrid(vlistID1, index);
	  int gridsize = gridInqSize(gridID1);
	  if ( gridsize == masksize )
	    {
	      int *mask = (int*) Malloc(masksize*sizeof(int));
	      for ( int i = 0; i < masksize; i++ )
		{
		  if ( gridmask[i] < 0 || gridmask[i] > 255 )
		    mask[i] = 0;
		  else
		    mask[i] = (int)lround(gridmask[i]);
		}
	      gridID2 = gridDuplicate(gridID1);
	      gridDefMask(gridID2, mask);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      Free(mask);
	    }
	}
    }
  else if ( operatorID == UNSETGRIDMASK )
    {
      int ngrids = vlistNgrids(vlistID1);
      for ( int index = 0; index < ngrids; index++ )
	{
	  int gridID1  = vlistGrid(vlistID1, index);
	  gridID2 = gridDuplicate(gridID1);
	  gridDefMask(gridID2, NULL);
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);
  //vlistPrint(vlistID2);

  int gridsize = (lregular || lregularnn) ? vlistGridsizeMax(vlistID2) : vlistGridsizeMax(vlistID1);

  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);
	  
	  pstreamReadRecord(streamID1, array, &nmiss);

	  int gridID1 = vlistInqVarGrid(vlistID1, varID);
	  if ( lregular || lregularnn )
	    {
	      gridID2 = vlistInqVarGrid(vlistID2, varID);
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  double missval = vlistInqVarMissval(vlistID1, varID);
                  int lnearst = lregularnn ? 1 : 0;
                  field2regular(gridID1, gridID2, missval, array, nmiss, lnearst);
		}
	    }
	  else if ( gridInqType(gridID1) == GRID_GME )
	    {
	      int gridsize = gridInqSize(gridID1);
	      int j = 0;
	      for ( int i = 0; i < gridsize; i++ )
		if ( grid2_vgpm[i] ) array[j++] = array[i];
	    }

	  pstreamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( gridmask ) Free(gridmask);
  if ( areaweight ) Free(areaweight);
  if ( array ) Free(array);
  if ( grid2_vgpm ) Free(grid2_vgpm);

  cdoFinish();

  return 0;
}
