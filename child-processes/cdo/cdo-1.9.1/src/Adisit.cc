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

      Adisit      adisit          compute insitu from potential temperature
      Adisit      adipot          compute potential from insitu temperature
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


/*
!>
!! transformation from potential to in situ temperature
!! according to Bryden, 1973, "New polynomials for thermal expansion,
!! adiabatic temperature gradient and potential temperature of sea
!! water". Deep Sea Research and Oceanographic Abstracts. 20, 401-408
!! (GILL P.602), which gives the inverse transformation for an
!! approximate value, all terms linear in t are taken after that one
!! newton step.  for the check value 8.4678516 the accuracy is 0.2
!! mikrokelvin.
!!
*/

/* compute insitu temperature from potential temperature */
static
double adisit_1(double tpot, double sal, double p)
{
  double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9,
         a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
         a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10,
         a_d = 4.1057E-9,
         a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double qc = p * (a_a1 + p * (a_c1 - a_e1 * p));
  double qv = p * (a_b1 - a_d * p);
  double dc = 1. + p * (-a_a2 + p * (a_c2 - a_e2 * p));
  double dv = a_b2 * p;
  double qnq  = -p * (-a_a3 + p * a_c3);
  double qn3  = -p * a_a4;

  double tpo = tpot;
  double qvs = qv*(sal - 35.) + qc;
  double dvs = dv*(sal - 35.) + dc;
  double t   = (tpo + qvs)/dvs;
  double fne = - qvs + t*(dvs + t*(qnq + t*qn3)) - tpo;
  double fst = dvs + t*(2.*qnq + 3.*qn3*t);
  t = t - fne/fst;

  return t;
}

/* compute potential temperature from insitu temperature */
/* Ref: Gill, p. 602, Section A3.5:Potential Temperature */
static
double adipot(double t, double s, double p)
{
  double a_a1 = 3.6504E-4, a_a2 = 8.3198E-5, a_a3 = 5.4065E-7, a_a4 = 4.0274E-9,
         a_b1 = 1.7439E-5, a_b2 = 2.9778E-7,
         a_c1 = 8.9309E-7, a_c2 = 3.1628E-8, a_c3 = 2.1987E-10,
         a_d = 4.1057E-9,
         a_e1 = 1.6056E-10, a_e2 = 5.0484E-12;

  double s_rel = s - 35.0;

  double aa = (a_a1+ t*(a_a2 - t*(a_a3 - a_a4*t)));
  double bb = s_rel*(a_b1 -a_b2*t)     ;
  double cc = (a_c1 + t*(-a_c2 + a_c3*t));
  double cc1 = a_d*s_rel;
  double dd = (-a_e1 + a_e2*t);

  double tpot = t-p*(aa + bb + p*(cc - cc1 + p*dd));

  return tpot;
}

static
void calc_adisit(long gridsize, long nlevel, double *pressure, field_type tho, field_type sao, field_type tis)
{
  /* pressure units: hPa     */
  /* tho units:      Celsius */
  /* sao units:      psu     */

  for ( long levelID = 0; levelID < nlevel; ++levelID )
    {
      long offset = gridsize*levelID;
      double *thoptr = tho.ptr + offset;
      double *saoptr = sao.ptr + offset;
      double *tisptr = tis.ptr + offset;

      for ( long i = 0; i < gridsize; ++i )
	{
	  if ( DBL_IS_EQUAL(thoptr[i], tho.missval) ||
	       DBL_IS_EQUAL(saoptr[i], sao.missval) )
	    {
	      tisptr[i] = tis.missval;
	    }
	  else
	    {
	      tisptr[i] = adisit_1(thoptr[i], saoptr[i], pressure[levelID]); 
	    }
	}
    }
}

static
void calc_adipot(long gridsize, long nlevel, double *pressure, field_type t, field_type s, field_type tpot)
{
  /* pressure units: hPa     */
  /* t units:      Celsius */
  /* s units:      psu     */

  for ( long levelID = 0; levelID < nlevel; ++levelID )
    {
      long offset = gridsize*levelID;
      double *tptr = t.ptr + offset;
      double *sptr = s.ptr + offset;
      double *tpotptr = tpot.ptr + offset;

      for ( long i = 0; i < gridsize; ++i )
	{
	  if ( DBL_IS_EQUAL(tptr[i], t.missval) ||
	       DBL_IS_EQUAL(sptr[i], s.missval) )
	    {
	      tpotptr[i] = tpot.missval;
	    }
	  else
	    {
	      tpotptr[i] = adipot(tptr[i], sptr[i], pressure[levelID]); 
	    }
	}
    }
}


void *Adisit(void *argument)
{
  int nrecs;
  int varID, levelID;
  int offset;
  int i;
  int nmiss;
  int thoID = -1, saoID = -1;
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  double pin = -1;
  double *single;

  cdoInitialize(argument);
  int ADISIT = cdoOperatorAdd("adisit", 1, 1, "");
  int ADIPOT = cdoOperatorAdd("adipot", 1, 1, "");

  UNUSED(ADIPOT);

  int operatorID = cdoOperatorID();

  if ( operatorArgc() == 1 ) pin = parameter2double(operatorArgv()[0]);
  
  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  int nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      int code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  vlistInqVarStdname(vlistID1,varID, stdname);
	  strtolower(varname);

               if ( strcmp(varname, "s")     == 0 ) code = 5;
          else if ( strcmp(varname, "t")     == 0 ) code = 2;
	  else if ( strcmp(stdname, "sea_water_salinity") == 0 ) code = 5;
          
          if ( operatorID == ADISIT )
          {
	   if ( strcmp(stdname, "sea_water_potential_temperature") == 0 ) code = 2;
          }
          else {
	   if ( strcmp(stdname, "sea_water_temperature") == 0 ) code = 2;
          }
	}

      if      ( code == 2 ) thoID = varID;
      else if ( code == 5 ) saoID = varID;
    }

  if ( saoID == -1 ) cdoAbort("Sea water salinity not found!");
  if ( thoID == -1 ) cdoAbort("Potential or Insitu temperature not found!");

  int gridID = vlistGrid(vlistID1, 0);
  int gridsize = vlist_check_gridsize(vlistID1);

  int zaxisID = vlistInqVarZaxis(vlistID1, saoID);
  int nlevel1 = zaxisInqSize(zaxisID);
      zaxisID = vlistInqVarZaxis(vlistID1, thoID);
  int nlevel2 = zaxisInqSize(zaxisID);

  if ( nlevel1 != nlevel2 ) cdoAbort("temperature and salinity have different number of levels!");
  int nlevel = nlevel1;

  double *pressure = (double*) Malloc(nlevel*sizeof(double));
  cdoZaxisInqLevels(zaxisID, pressure);

  if ( pin >= 0 ) 
    for ( i = 0; i < nlevel; ++i ) pressure[i] = pin;
  else
    for ( i = 0; i < nlevel; ++i ) pressure[i] /= 10;

  if ( cdoVerbose )
    {
      cdoPrint("Level Pressure");
      for ( i = 0; i < nlevel; ++i )
	cdoPrint("%5d  %g", i+1, pressure[i]);
    }

  field_type tho, sao, tis;
  field_init(&tho);
  field_init(&sao);
  field_init(&tis);
  tho.ptr = (double*) Malloc(gridsize*nlevel*sizeof(double));
  sao.ptr = (double*) Malloc(gridsize*nlevel*sizeof(double));
  tis.ptr = (double*) Malloc(gridsize*nlevel*sizeof(double));

  tho.nmiss = 0;
  sao.nmiss = 0;
  tis.nmiss = 0;
  
  tho.missval = vlistInqVarMissval(vlistID1, thoID);
  sao.missval = vlistInqVarMissval(vlistID1, saoID);
  tis.missval = tho.missval;

  int datatype = CDI_DATATYPE_FLT32;
  if ( vlistInqVarDatatype(vlistID1, thoID) == CDI_DATATYPE_FLT64 &&
       vlistInqVarDatatype(vlistID1, saoID) == CDI_DATATYPE_FLT64 )
    datatype = CDI_DATATYPE_FLT64;

  int vlistID2 = vlistCreate();

  int tisID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);
  if ( operatorID == ADISIT )
    {
      vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(20, 255, 255));
      vlistDefVarName(vlistID2, tisID2, "to");
      vlistDefVarLongname(vlistID2, tisID2, "Sea water temperature");
      vlistDefVarStdname(vlistID2, tisID2, "sea_water_temperature");
    }
  else
    {
      vlistDefVarParam(vlistID2, tisID2, cdiEncodeParam(2, 255, 255));
      vlistDefVarName(vlistID2, tisID2, "tho");
      vlistDefVarLongname(vlistID2, tisID2, "Sea water potential temperature");
      vlistDefVarStdname(vlistID2, tisID2, "sea_water_potential_temperature");
    }
  vlistDefVarUnits(vlistID2, tisID2, "K");
  vlistDefVarMissval(vlistID2, tisID2, tis.missval);
  vlistDefVarDatatype(vlistID2, tisID2, datatype);

  int saoID2 = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);
  vlistDefVarParam(vlistID2, saoID2, cdiEncodeParam(5, 255, 255));
  vlistDefVarName(vlistID2, saoID2, "s");
  vlistDefVarLongname(vlistID2, saoID2, "Sea water salinity");
  vlistDefVarStdname(vlistID2, saoID2, "sea_water_salinity");
  vlistDefVarUnits(vlistID2, saoID2, "psu");
  vlistDefVarMissval(vlistID2, saoID2, sao.missval);
  vlistDefVarDatatype(vlistID2, saoID2, datatype);


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
	       
      for ( int recID = 0; recID < nrecs; ++recID )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  offset = gridsize*levelID;

	  if ( varID == thoID )
            {
              pstreamReadRecord(streamID1, tho.ptr+offset, &nmiss);
              tho.nmiss = (size_t) nmiss;
            }
	  if ( varID == saoID )
            {
              pstreamReadRecord(streamID1, sao.ptr+offset, &nmiss);
              sao.nmiss = (size_t) nmiss;
            }
        }

      if ( operatorID == ADISIT )
        calc_adisit(gridsize, nlevel, pressure, tho, sao, tis); 
      else
        calc_adipot(gridsize, nlevel, pressure, tho, sao, tis); 


      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  offset = gridsize*levelID;

	  single = tis.ptr+offset;

	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(single[i], tis.missval) ) nmiss++;
 
	  pstreamDefRecord(streamID2, tisID2, levelID);
	  pstreamWriteRecord(streamID2, single, nmiss);     

	  single = sao.ptr+offset;

	  nmiss = 0;
	  for ( i = 0; i < gridsize; ++i )
	    if ( DBL_IS_EQUAL(single[i], sao.missval) ) nmiss++;
 
	  pstreamDefRecord(streamID2, saoID2, levelID);
	  pstreamWriteRecord(streamID2, single, nmiss);     
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  Free(pressure);
  Free(tis.ptr);
  Free(tho.ptr);
  Free(sao.ptr);

  cdoFinish();

  return 0;
}
