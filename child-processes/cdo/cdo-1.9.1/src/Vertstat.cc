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

      Vertstat   vertrange       Vertical range
      Vertstat   vertmin         Vertical minimum
      Vertstat   vertmax         Vertical maximum
      Vertstat   vertsum         Vertical sum
      Vertstat   vertint         Vertical integral
      Vertstat   vertmean        Vertical mean
      Vertstat   vertavg         Vertical average
      Vertstat   vertvar         Vertical variance
      Vertstat   vertvar1        Vertical variance [Normalize by (n-1)]
      Vertstat   vertstd         Vertical standard deviation
      Vertstat   vertstd1        Vertical standard deviation [Normalize by (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define IS_SURFACE_LEVEL(zaxisID)  (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)


int getSurfaceID(int vlistID)
{
  int surfID = -1;
  int nzaxis = vlistNzaxis(vlistID);

  for ( int index = 0; index < nzaxis; ++index )
    {
      int zaxisID = vlistZaxis(vlistID, index);
      if ( IS_SURFACE_LEVEL(zaxisID) )
	{
	  surfID = vlistZaxis(vlistID, index);
	  break;
	}
    }

  if ( surfID == -1 ) surfID = zaxisCreate(ZAXIS_SURFACE, 1);

  return surfID;
}

static
void setSurfaceID(int vlistID, int surfID)
{
  int nzaxis = vlistNzaxis(vlistID);

  for ( int index = 0; index < nzaxis; ++index )
    {
      int zaxisID = vlistZaxis(vlistID, index);
      if ( zaxisID != surfID || !IS_SURFACE_LEVEL(zaxisID) )
	vlistChangeZaxisIndex(vlistID, index, surfID);
    }
}


void genLayerBounds(int nlev, double *levels, double *lbounds, double *ubounds)
{
  if ( nlev == 1 )
    {
      lbounds[0] = 0.;
      ubounds[0] = 1.;
    }
  else
    {
      lbounds[0]      = levels[0];
      ubounds[nlev-1] = levels[nlev-1];
      for ( int i = 0; i < nlev-1; ++i )
        {
          double bound = 0.5*(levels[i] + levels[i+1]);
          lbounds[i+1] = bound;
          ubounds[i]   = bound;
        }
    }
}


int getLayerThickness(bool genbounds, int index, int zaxisID, int nlev, double *thickness, double *weights)
{
  int status = 0;
  int i;
  double *levels  = (double *) Malloc(nlev*sizeof(double));
  double *lbounds = (double *) Malloc(nlev*sizeof(double));
  double *ubounds = (double *) Malloc(nlev*sizeof(double));

  cdoZaxisInqLevels(zaxisID, levels);
  if ( genbounds )
    {
      status = 2;
      genLayerBounds(nlev, levels, lbounds, ubounds);
    }
  else if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      status = 1;
      zaxisInqLbounds(zaxisID, lbounds);
      zaxisInqUbounds(zaxisID, ubounds);
    }
  else
    {
      for ( i = 0; i < nlev; ++i )
	{
	  lbounds[i] = 0.;
	  ubounds[i] = 1.;
	}
    }
  
  for ( i = 0; i < nlev; ++i ) thickness[i] = fabs(ubounds[i]-lbounds[i]);

  double lsum = 0;
  for ( i = 0; i < nlev; ++i ) lsum += thickness[i];

  for ( i = 0; i < nlev; ++i ) weights[i] = thickness[i];
  
  for ( i = 0; i < nlev; ++i ) weights[i] /= (lsum/nlev);

  double wsum = 0;
  for ( i = 0; i < nlev; ++i ) wsum += weights[i];

  if ( cdoVerbose )
    {
      cdoPrint("zaxisID=%d  nlev=%d  layersum=%g  weightsum=%g", index, nlev, lsum, wsum);
      printf("         level     bounds   thickness  weight\n");
      for ( i = 0; i < nlev; ++i )
	printf("   %3d  %6g  %6g/%-6g  %6g  %6g\n", i+1, levels[i], lbounds[i], ubounds[i], thickness[i], weights[i]);
    }

  Free(levels);
  Free(lbounds);
  Free(ubounds);

  return status;
}


void *Vertstat(void *argument)
{
  int nrecs;
  int gridID;
  int varID, levelID;
  int nmiss;
  typedef struct {
    int zaxisID;
    int status;
    int numlevel;
    double *thickness;
    double *weights;
  }
  vert_t;

  cdoInitialize(argument);

  // clang-format off
                 cdoOperatorAdd("vertrange", func_range, 0, NULL);
                 cdoOperatorAdd("vertmin",   func_min,   0, NULL);
                 cdoOperatorAdd("vertmax",   func_max,   0, NULL);
                 cdoOperatorAdd("vertsum",   func_sum,   0, NULL);
  int VERTINT  = cdoOperatorAdd("vertint",   func_sum,   1, NULL);
                 cdoOperatorAdd("vertmean",  func_mean,  1, NULL);
                 cdoOperatorAdd("vertavg",   func_avg,   1, NULL);
                 cdoOperatorAdd("vertvar",   func_var,   1, NULL);
                 cdoOperatorAdd("vertvar1",  func_var1,  1, NULL);
                 cdoOperatorAdd("vertstd",   func_std,   1, NULL);
                 cdoOperatorAdd("vertstd1",  func_std1,  1, NULL);

  int operatorID   = cdoOperatorID();
  int operfunc     = cdoOperatorF1(operatorID);
  bool needWeights = cdoOperatorF2(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  //int applyWeights = lmean;

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    vlistDefFlag(vlistID1, varID, 0, TRUE);

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int surfID = getSurfaceID(vlistID1);
  setSurfaceID(vlistID2, surfID);

  int nzaxis = vlistNzaxis(vlistID1);
  int nlev, zaxisID;
  vert_t *vert = (vert_t*) Malloc(nzaxis*sizeof(vert_t));
  if ( needWeights )
    {
      bool genbounds = false;
      unsigned npar = operatorArgc();
      if ( npar > 0 )
	{
	  char **parnames = operatorArgv();

	  if ( cdoVerbose )
	    for ( unsigned i = 0; i < npar; i++ )
	      cdoPrint("key %d = %s", i+1, parnames[i]);

	  if ( strcmp(parnames[0], "genbounds") == 0 ) genbounds = true;
	  else cdoAbort("Parameter >%s< unsupported! Supported parameter are: genbounds", parnames[0]);
	}
      
      for ( int index = 0; index < nzaxis; ++index )
	{
	  zaxisID = vlistZaxis(vlistID1, index);
	  nlev = zaxisInqSize(zaxisID);
	  vert[index].numlevel = 0;
	  vert[index].status   = 0;
	  vert[index].zaxisID  = zaxisID;
	  // if ( nlev > 1 )
	    {
	      vert[index].numlevel = nlev;
	      vert[index].thickness = (double *) Malloc(nlev*sizeof(double));
	      vert[index].weights = (double *) Malloc(nlev*sizeof(double));
	      vert[index].status = getLayerThickness(genbounds, index, zaxisID, nlev, vert[index].thickness, vert[index].weights); 
	    }
	}
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type *vars1 = (field_type*) Malloc(nvars*sizeof(field_type));
  field_type *samp1 = (field_type*) Malloc(nvars*sizeof(field_type));
  field_type *vars2 = NULL;
  if ( lvarstd || lrange )
    vars2 = (field_type*) Malloc(nvars*sizeof(field_type));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      double missval = vlistInqVarMissval(vlistID1, varID);

      field_init(&vars1[varID]);
      field_init(&samp1[varID]);
      vars1[varID].grid    = gridID;
      vars1[varID].zaxis   = zaxisID;
      vars1[varID].nsamp   = 0;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double*) Malloc(gridsize*sizeof(double));
      samp1[varID].grid    = gridID;
      samp1[varID].nmiss   = 0;
      samp1[varID].missval = missval;
      samp1[varID].ptr     = NULL;
      if ( lvarstd || lrange )
	{
	  field_init(&vars2[varID]);
	  vars2[varID].grid    = gridID;
	  vars2[varID].nmiss   = 0;
	  vars2[varID].missval = missval;
	  vars2[varID].ptr     = (double*) Malloc(gridsize*sizeof(double));
	}
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

          vars1[varID].nsamp++;
          if ( lrange ) vars2[varID].nsamp++;
	  gridsize = gridInqSize(vars1[varID].grid);
	  zaxisID  = vars1[varID].zaxis;
	  nlev = zaxisInqSize(zaxisID);

	  double layer_weight = 1.0;
	  double layer_thickness = 1.0;
	  if ( needWeights )
	    {
	      for ( int index = 0; index < nzaxis; ++index )
		if ( vert[index].zaxisID == zaxisID )
		  {		    
		    if ( vert[index].status == 0 && tsID == 0 && levelID == 0 && nlev > 1 )
		      {
			char varname[CDI_MAX_NAME];
			vlistInqVarName(vlistID1, varID, varname);
			cdoWarning("Layer bounds not available, using constant vertical weights for variable %s!", varname);
		      }
		    else
		      {
			layer_weight    = vert[index].weights[levelID];
			layer_thickness = vert[index].thickness[levelID];
		      }

		    break;
		  }
	    }

	  if ( levelID == 0 )
	    {
	      pstreamReadRecord(streamID1, vars1[varID].ptr, &nmiss);
	      vars1[varID].nmiss = (size_t)nmiss;
              if ( lrange )
                {
                  vars2[varID].nmiss = (size_t)nmiss;
                  for ( int i = 0; i < gridsize; i++ )
                    vars2[varID].ptr[i] = vars1[varID].ptr[i];
                }

	      if ( operatorID == VERTINT && IS_NOT_EQUAL(layer_thickness, 1.0) ) farcmul(&vars1[varID], layer_thickness);
	      if ( lmean && IS_NOT_EQUAL(layer_weight, 1.0) ) farcmul(&vars1[varID], layer_weight);

	      if ( lvarstd )
                {
                  if ( IS_NOT_EQUAL(layer_weight, 1.0) )
                    {
                      farmoqw(&vars2[varID], vars1[varID], layer_weight);
                      farcmul(&vars1[varID], layer_weight);
                    }
                  else
                    {
                      farmoq(&vars2[varID], vars1[varID]);
                    }
                }

	      if ( nmiss > 0 || samp1[varID].ptr || needWeights )
		{
		  if ( samp1[varID].ptr == NULL )
		    samp1[varID].ptr = (double *) Malloc(gridsize*sizeof(double));

		  for ( int i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[varID].ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i] = 0.;
		    else
		      samp1[varID].ptr[i] = layer_weight;
		}
	    }
	  else
	    {
	      pstreamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t)nmiss;
	      field.grid    = vars1[varID].grid;
	      field.missval = vars1[varID].missval;

	      if ( operatorID == VERTINT && IS_NOT_EQUAL(layer_thickness, 1.0) ) farcmul(&field, layer_thickness);
	      if ( lmean && IS_NOT_EQUAL(layer_weight, 1.0) ) farcmul(&field, layer_weight);

	      if ( field.nmiss > 0 || samp1[varID].ptr )
		{
		  if ( samp1[varID].ptr == NULL )
		    {
		      samp1[varID].ptr = (double*) Malloc(gridsize*sizeof(double));
		      for ( int i = 0; i < gridsize; i++ )
			samp1[varID].ptr[i] = vars1[varID].nsamp;
		    }

		  for ( int i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i] += layer_weight;
		}

	      if ( lvarstd )
		{
                  if ( IS_NOT_EQUAL(layer_weight, 1.0) )
                    {
                      farsumqw(&vars2[varID], field, layer_weight);
                      farsumw(&vars1[varID], field, layer_weight);
                    }
                  else
                    {
                      farsumq(&vars2[varID], field);
                      farsum(&vars1[varID], field);
                    }
		}
              else if ( lrange )
                {
                  farmin(&vars2[varID], field);
                  farmax(&vars1[varID], field);
                }
	      else
		{
		  farfun(&vars1[varID], field, operfunc);
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars1[varID].nsamp )
	    {
	      if ( lmean )
		{
		  if ( samp1[varID].ptr == NULL )
		    farcdiv(&vars1[varID], (double)vars1[varID].nsamp);
		  else
		    fardiv(&vars1[varID], samp1[varID]);
		}
	      else if ( lvarstd )
		{
		  if ( samp1[varID].ptr == NULL )
		    {
		      if ( lstd )
			farcstd(&vars1[varID], vars2[varID], vars1[varID].nsamp, divisor);
		      else
			farcvar(&vars1[varID], vars2[varID], vars1[varID].nsamp, divisor);
		    }
		  else
		    {
		      if ( lstd )
			farstd(&vars1[varID], vars2[varID], samp1[varID], divisor);
		      else
			farvar(&vars1[varID], vars2[varID], samp1[varID], divisor);
		    }
		}
              else if ( lrange )
                {
                  farsub(&vars1[varID], vars2[varID]);
                }

	      pstreamDefRecord(streamID2, varID, 0);
	      pstreamWriteRecord(streamID2, vars1[varID].ptr, (int)vars1[varID].nmiss);
	      vars1[varID].nsamp = 0;
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(vars1[varID].ptr);
      if ( samp1[varID].ptr ) Free(samp1[varID].ptr);
      if ( lvarstd ) Free(vars2[varID].ptr);
    }

  Free(vars1);
  Free(samp1);
  if ( lvarstd ) Free(vars2);

  if ( field.ptr ) Free(field.ptr);

  if ( needWeights )
    for ( int index = 0; index < nzaxis; ++index )
      if ( vert[index].numlevel > 1 )  Free(vert[index].weights);

  Free(vert);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
