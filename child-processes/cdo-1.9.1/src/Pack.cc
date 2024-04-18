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

      Pack    pack         Pack
*/

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024


static
int get_type_values(int datatype, double *tmin, double *tmax)
{
  int status = 0;

  switch (datatype)
    {
    case CDI_DATATYPE_INT8:    *tmin = -SCHAR_MAX+1; *tmax = SCHAR_MAX;   break;
    case CDI_DATATYPE_UINT8:   *tmin = 0;            *tmax = UCHAR_MAX-1; break;
    case CDI_DATATYPE_INT16:   *tmin = -SHRT_MAX+1;  *tmax = SHRT_MAX;    break;
    case CDI_DATATYPE_UINT16:  *tmin = 0;            *tmax = USHRT_MAX-1; break;
    case CDI_DATATYPE_INT32:   *tmin = -INT_MAX+1;   *tmax = INT_MAX;     break;
    case CDI_DATATYPE_UINT32:  *tmin = 0;            *tmax = UINT_MAX-1;  break;
    default: status = 1; break;
    }

  return status;
}

static
int compute_scale(int datatype, double fmin, double fmax, double *scale_factor, double *add_offset)
{
  double tmin, tmax;
  double ao = 0.0, sf = 1.0;

  *scale_factor = sf;
  *add_offset = ao;

  if ( get_type_values(datatype, &tmin, &tmax) ) return 1;

  if ( IS_NOT_EQUAL(fmin, fmax) )
    {
      sf = (fmax - fmin) / (tmax - tmin);
      ao = ((fmax + fmin) - sf * (tmin + tmax)) / 2;
    }

  *scale_factor = sf;
  *add_offset   = ao;

  return 0;
}


void *Pack(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID;
  int nalloc = 0;
  int nmiss;
  int nlevel;
  int datatype = CDI_DATATYPE_INT16;
  dtlist_type *dtlist = dtlist_new();
  double missval1, missval2;
  field_type ***vars = NULL;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vars   = (field_type ***) Realloc(vars, nalloc*sizeof(field_type **));
	}

      dtlist_taxisInqTimestep(dtlist, taxisID1, tsID);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
	  pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  int nts = tsID;

  if ( cdoDefaultDataType != CDI_UNDEFID )
    {
      if ( cdoDefaultDataType == CDI_DATATYPE_FLT64 || cdoDefaultDataType == CDI_DATATYPE_FLT32 )
	{
	  cdoWarning("Changed default output datatype to int16");
	  cdoDefaultDataType = datatype;
	}
      else
	{
	  datatype = cdoDefaultDataType;
	}
    }

  cdoDefaultDataType = datatype;


  for ( varID = 0; varID < nvars; varID++ )
    {
      double fmin =  1.e300;
      double fmax = -1.e300;
      double sf, ao;
      long ivals = 0;
      long nmisspv = 0;

      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval1 = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  for ( tsID = 0; tsID < nts; tsID++ )
	    {
	      double *array =  vars[tsID][varID][levelID].ptr;
	      nmiss = vars[tsID][varID][levelID].nmiss;
	      if ( nmiss > 0 )
		{
		  nmisspv += nmiss;
		  for ( int i = 0; i < gridsize; ++i )
		    {
		      if ( !DBL_IS_EQUAL(array[i], missval1) )
			{
			  if ( array[i] < fmin ) fmin = array[i];
			  if ( array[i] > fmax ) fmax = array[i];
			  ivals++;
			}
		    }
		}
	      else
		{
		  for ( int i = 0; i < gridsize; ++i )
		    {
		      if ( array[i] < fmin ) fmin = array[i];
		      if ( array[i] > fmax ) fmax = array[i];
		    }
		  ivals = gridsize;
		}
	    }
	}

      vlistDefVarDatatype(vlistID2, varID, datatype);
      missval2 = vlistInqVarMissval(vlistID2, varID);

      if ( nmisspv > 0 )
	{
	  double tmin, tmax;
	  if ( !get_type_values(datatype, &tmin, &tmax) )
	    {
	      if ( !(missval2 < tmin || missval2 > tmax) )
		cdoWarning("new missing value %g is inside data range (%g - %g)!", missval2, tmin, tmax); 

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  for ( tsID = 0; tsID < nts; tsID++ )
		    {
		      double *array =  vars[tsID][varID][levelID].ptr;
		      nmiss = vars[tsID][varID][levelID].nmiss;
		      if ( nmiss > 0 )
			for ( int i = 0; i < gridsize; ++i )
			  if ( DBL_IS_EQUAL(array[i], missval1) ) array[i] = missval2;
		    }
		}	  
	    }
	}

      // printf("fmin %g fmax %g missval %g\n", fmin, fmax, vlistInqVarMissval(vlistID2, varID));
      if ( !compute_scale(datatype, fmin, fmax, &sf, &ao) )
	{
	  // printf("sf = %g ao = %g \n", sf, ao);
	  // printf("smin %g smax %g\n", (fmin - ao)/sf, (fmax -ao)/sf);

	  vlistDefVarScalefactor(vlistID2, varID, sf);
	  vlistDefVarAddoffset(vlistID2, varID, ao);
	}
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      dtlist_taxisDefTimestep(dtlist, taxisID2, tsID);
      pstreamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( vars[tsID][varID][levelID].ptr )
		{
		  nmiss = vars[tsID][varID][levelID].nmiss;
		  pstreamDefRecord(streamID2, varID, levelID);
		  pstreamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		  Free(vars[tsID][varID][levelID].ptr);
		  vars[tsID][varID][levelID].ptr = NULL;
		}
	    }
	}

      field_free(vars[tsID], vlistID1);
    }

  if ( vars  ) Free(vars);

  dtlist_delete(dtlist);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
