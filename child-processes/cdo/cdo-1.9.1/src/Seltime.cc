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

      Seltime    seltimestep     Select timesteps
      Seltime    seltime         Select times
      Seltime    selhour         Select hours
      Seltime    selday          Select days
      Seltime    selmonth        Select months
      Seltime    selyear         Select years
      Seltime    selseason       Select seasons
      Seltime    seldate         Select dates
      Seltime    selsmon         Select single month
*/


#include <cdi.h>
#include "cdo_int.h"
#include "pstream.h"
#include "listarray.h"


void season_to_months(const char *season, int *imonths)
{
  const char *smons = "JFMAMJJASONDJFMAMJJASOND";
  const int imons[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  assert(strlen(smons)==(sizeof(imons)/sizeof(int)));

  size_t len = strlen(season);
  if ( len == 3 && strcmp(season, "ANN") == 0 )
    {
      for ( size_t k = 0; k < 12; ++k ) imonths[k+1] = 1;
    }
  else
    {
      if ( len > 12 ) cdoAbort("Too many months %d (limit=12)!", (int)len);
      char *season_u = strdup(season);
      strtoupper(season_u);
      const char *sstr = strstr(smons, season_u);
      free(season_u);
      if ( sstr == NULL ) cdoAbort("Season %s not available!", season);
      size_t ks = (size_t)(sstr-smons);
      size_t ke = ks + len;
      for ( size_t k = ks; k < ke; ++k ) imonths[imons[k]]++;
    }
}

static
int seaslist(lista_t *ilista)
{
  int imon[13]; /* 1-12 ! */
  for ( int i = 0; i < 13; ++i ) imon[i] = 0;

  int nsel = operatorArgc();
  if ( isdigit(*operatorArgv()[0]))
    {
      for ( int i = 0; i < nsel; i++ )
        {
          int ival = parameter2int(operatorArgv()[i]);
          if      ( ival == 1 ) { imon[12]++; imon[ 1]++; imon[ 2]++; }
          else if ( ival == 2 ) { imon[ 3]++; imon[ 4]++; imon[ 5]++; }
          else if ( ival == 3 ) { imon[ 6]++; imon[ 7]++; imon[ 8]++; }
          else if ( ival == 4 ) { imon[ 9]++; imon[10]++; imon[11]++; }
          else cdoAbort("Season %d not available!", ival);
        }
    }
  else
    {
      for ( int i = 0; i < nsel; i++ )
        season_to_months(operatorArgv()[i], imon);
    }

  nsel = 0;
  for ( int i = 1; i < 13; ++i ) if ( imon[i] ) lista_set_int(ilista, nsel++, i);

  return nsel;
}


double datestr_to_double(const char *datestr, int opt)
{
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  double fval = 0;

  if ( opt ) {hour = 23; minute = 59; second = 59;}

  if ( strchr(datestr, '-') == NULL )
    {
      fval = parameter2double(datestr);
    }
  else if ( strchr(datestr, 'T') )
    {
      sscanf(datestr, "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
      fval = cdiEncodeTime(hour, minute, second);
      if ( fabs(fval) > 0 ) fval /= 1000000;
      fval += cdiEncodeDate(year, month, day);
    }
  else
    {
      sscanf(datestr, "%d-%d-%d", &year, &month, &day);
      fval = cdiEncodeTime(hour, minute, second);
      if ( fabs(fval) > 0 ) fval /= 1000000;
      fval += cdiEncodeDate(year, month, day);
    }

  return fval;
}

static
int datelist(lista_t *flista)
{
  bool set2 = true;
  double fval = 0;

  int nsel = operatorArgc();
  if ( nsel < 1 ) cdoAbort("Too few arguments!");
  if ( nsel > 2 ) cdoAbort("Too many arguments!");

  for ( int i = 0; i < nsel; i++ )
    {
      if ( operatorArgv()[i][0] == '-' && operatorArgv()[i][1] == 0 )
	{
	  if ( i == 0 )
	    fval = -99999999999.;
	  else
	    fval =  99999999999.;

	  lista_set_flt(flista, i,  fval);
	}
      else
	{
	  fval = datestr_to_double(operatorArgv()[i], 0);
	  if ( strchr(operatorArgv()[i], 'T') )
	    set2 = false;
	  else if ( nsel > 1 && i > 0 )
	    fval += 0.999;
	}

      lista_set_flt(flista, i, fval);
    }

  if ( nsel == 1 && set2 )
    {
      fval += 0.999;
      lista_set_flt(flista, nsel, fval);
      nsel = 2;
    }

  return nsel;
}


void *Seltime(void *argument)
{
  int streamID2 = -1;
  int nrecs;
  int varID, levelID;
  int nsel = 0;
  int i;
  int nmiss;
  int ncts = 0, nts, it;
  int hour = 0, minute = 0, second = 0;
  int nts1 = 0, nts2 = 0;
  int its1 = 0, its2 = 0;
  double selfval = 0;
  double *array = NULL;
  lista_t *ilista = lista_new(INT_LISTA);
  lista_t *flista = lista_new(FLT_LISTA);
  bool copy_nts2 = false;
  bool lconstout = false;
  bool process_nts1 = false, process_nts2 = false;
  bool *selfound = NULL;
  int *vdate_list = NULL, *vtime_list = NULL;
  double *single;
  field_type ***vars = NULL;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int SELTIMESTEP = cdoOperatorAdd("seltimestep", func_step,     0, "timesteps");
  int SELDATE     = cdoOperatorAdd("seldate",     func_datetime, 0, "start date and end date (format YYYY-MM-DDThh:mm:ss)");
  int SELTIME     = cdoOperatorAdd("seltime",     func_time,     0, "times (format hh:mm:ss)");
  int SELHOUR     = cdoOperatorAdd("selhour",     func_time,     0, "hours");
  int SELDAY      = cdoOperatorAdd("selday",      func_date,     0, "days");
  int SELMONTH    = cdoOperatorAdd("selmonth",    func_date,     0, "months");
  int SELYEAR     = cdoOperatorAdd("selyear",     func_date,     0, "years");
  int SELSEASON   = cdoOperatorAdd("selseason",   func_date,     0, "seasons");
  int SELSMON     = cdoOperatorAdd("selsmon",     func_date,     0, "month[,nts1[,nts2]]");
  // clang-format on

  int operatorID = cdoOperatorID();

  int operfunc = cdoOperatorF1(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SELSEASON )
    {
      nsel = seaslist(ilista);
    }
  else if ( operatorID == SELDATE )
    {
      nsel = datelist(flista);
    }
  else if ( operatorID == SELTIME )
    {
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Too few arguments!");
      for ( i = 0; i < nsel; i++ )
	{
	  if ( strchr(operatorArgv()[i], ':') )
	    {
	      sscanf(operatorArgv()[i], "%d:%d:%d", &hour, &minute, &second);
	      lista_set_int(ilista, i, cdiEncodeTime(hour, minute, second));
	    }
	  else
	    {
	      lista_set_int(ilista, i, parameter2int(operatorArgv()[i]));
	    }
	}
    }
  else
    {
      nsel = args2int_lista(operatorArgc(), operatorArgv(), ilista);
    }

  if ( nsel < 1 ) cdoAbort("No timestep selected!");

  int *intarr = (int *) lista_dataptr(ilista);
  double *fltarr = (double *) lista_dataptr(flista);

  if ( operatorID == SELSMON )
    {
      if ( nsel > 1 ) nts1 = intarr[1];
      if ( nsel > 2 ) nts2 = intarr[2];
      else            nts2 = nts1;

      if ( nsel > 3 ) cdoAbort("Too many parameters");

      if ( cdoVerbose )
	cdoPrint("mon=%d  nts1=%d  nts2=%d", intarr[0], nts1, nts2);

      nsel = 1;
    }

  if ( nsel )
    {
      selfound = (bool*) Malloc(nsel*sizeof(bool));
      for ( i = 0; i < nsel; i++ ) selfound[i] = false;
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if ( nsel == 1 && operfunc == func_step )  vlistDefNtsteps(vlistID2,  1);
  else                                       vlistDefNtsteps(vlistID2, -1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( ! lcopy )
    {
      int gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int ntsteps = vlistNtsteps(vlistID1);

  /* add support for negative timestep values */
  if ( operatorID == SELTIMESTEP && ntsteps > 0 )
    {
      for ( i = 0; i < nsel; i++ )
	{
	  if ( intarr[i] < 0 )
	    {
	      if ( cdoVerbose )
		cdoPrint("timestep %d changed to %d", intarr[i], ntsteps + 1 + intarr[i]);
	      intarr[i] = ntsteps + 1 + intarr[i];
	    }
	}
    }

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	if ( operatorID == SELDATE )
	  cdoPrint("fltarr entry: %d %14.4f", i+1, fltarr[i]);
	else
	  cdoPrint("intarr entry: %d %d", i+1, intarr[i]);
    }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) nconst++;
      
  bool lnts1 = (operatorID == SELSMON) && (nts1 > 0);

  if ( lnts1 || nconst )
    {
      if ( lnts1 )
	{
	  vdate_list = (int*) Malloc(nts1*sizeof(int));
	  vtime_list = (int*) Malloc(nts1*sizeof(int));
	}
      else
	{
	  nts1 = 1;
	}

      vars  = (field_type ***) Malloc(nts1*sizeof(field_type **));

      for ( int tsID = 0; tsID < nts1; tsID++ )
	{
	  vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lnts1 || (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) )
		{
		  int gridID  = vlistInqVarGrid(vlistID1, varID);
		  int nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  int gridsize = gridInqSize(gridID);
		  
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      vars[tsID][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}
    }

  int tsID  = 0;
  int tsID2 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      bool copytimestep = false;
      int selival = -1;

      if ( operfunc == func_step )
	{
	  selival = tsID + 1;
	  if ( selival > intarr[nsel-1] ) break;
	}
      else if ( operfunc == func_date )
	{
          int year, month, day;
          cdiDecodeDate(vdate, &year, &month, &day);
          if      ( operatorID == SELYEAR ) selival = year;
          else if ( operatorID == SELDAY  ) selival = day;
          else                              selival = month;
	}
      else if ( operfunc == func_time )
	{
          int hour, minute, second;
          cdiDecodeTime(vtime, &hour, &minute, &second);
          if      ( operatorID == SELHOUR ) selival = hour;
          else if ( operatorID == SELTIME ) selival = vtime;
	}
      else if ( operfunc == func_datetime )
	{
	  selfval = vdate + vtime/1000000.;
	}

      if ( operatorID == SELDATE )
	{
	  if ( selfval >= fltarr[0] && selfval <= fltarr[nsel-1] )
	    {
	      copytimestep = true;
	      selfound[0]      = true;
	      selfound[nsel-1] = true;
	    }
          else if ( selfval > fltarr[nsel-1] )
            {
              break;
            }
	}
      else
	{
	  for ( i = 0; i < nsel; i++ )
	    if ( selival == intarr[i] )
	      {
		copytimestep = true;
		selfound[i] = true;
		break;
	      }
	}

      if ( operatorID == SELSMON && copytimestep == false )
	{
	  copy_nts2 = false;

	  if ( process_nts1 )
	    {
	      process_nts2 = true;
	      its2 = 0;
	      process_nts1 = false;
	    }

	  if ( process_nts2 )
	    {
	      if ( its2++ < nts2 )
		copy_nts2 = true;
	      else
		process_nts2 = false;
	    }
	}

      if ( copytimestep || copy_nts2 )
	{
	  if ( tsID2 == 0 )
	    {
	      streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
	      pstreamDefVlist(streamID2, vlistID2);
	    }

	  if ( lnts1 && ncts == 0 )
	    {
	      nts = nts1;
	      if ( its1 < nts1 )
		{
		  nts = its1;
		  cdoWarning("%d timesteps missing before month %d!", nts1-its1, intarr[0]);
		}

	      for ( it = 0; it < nts; it++ )
		{
		  taxisDefVdate(taxisID2, vdate_list[it]);
		  taxisDefVtime(taxisID2, vtime_list[it]);
		  pstreamDefTimestep(streamID2, tsID2++);
		  
		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT && tsID2 > 1 ) continue;
		      int nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  pstreamDefRecord(streamID2, varID, levelID);
			  single = vars[it][varID][levelID].ptr;
			  nmiss  = vars[it][varID][levelID].nmiss;
			  pstreamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}

	      its1 = 0;
	    }

	  ncts++;
	  if ( !process_nts2 )
	    {
	      its2 = 0;
	      process_nts1 = true;
	    }

	  taxisCopyTimestep(taxisID2, taxisID1);
	  pstreamDefTimestep(streamID2, tsID2++);

	  if ( tsID > 0 && lconstout )
	    {
	      lconstout = false;
	      nts = nts1 - 1;
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
		    {
		      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  pstreamDefRecord(streamID2, varID, levelID);
			  single = vars[nts][varID][levelID].ptr;
			  nmiss  = vars[nts][varID][levelID].nmiss;
			  pstreamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}
	    }

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);
	      pstreamDefRecord(streamID2, varID, levelID);
	      if ( lcopy )
		{
		  pstreamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  pstreamReadRecord(streamID1, array, &nmiss);
		  pstreamWriteRecord(streamID2, array, nmiss);
		}
	    }
	}
      else
	{
	  ncts = 0;

	  if ( lnts1 || tsID == 0 )
	    {
	      if ( tsID == 0 && nconst && (!lnts1) ) lconstout = true;

	      nts = nts1-1;
	      if ( lnts1 )
		{
		  if ( its1 <= nts )
		    nts = its1;
		  else
		    for ( it = 0; it < nts; it++ )
		      {
			vdate_list[it] = vdate_list[it+1];
			vtime_list[it] = vtime_list[it+1];
			for ( varID = 0; varID < nvars; varID++ )
			  {
			    if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT ) continue;
			    int gridID   = vlistInqVarGrid(vlistID1, varID);
			    int gridsize = gridInqSize(gridID);
			    int nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
			    for ( levelID = 0; levelID < nlevel; levelID++ )
			      {
				memcpy(vars[it][varID][levelID].ptr,
				       vars[it+1][varID][levelID].ptr,
				       gridsize*sizeof(double));
				vars[it][varID][levelID].nmiss = vars[it+1][varID][levelID].nmiss;
			      }
			  }
		      }

		  vdate_list[nts] = taxisInqVdate(taxisID1);
		  vtime_list[nts] = taxisInqVtime(taxisID1);

		  its1++;
		}

	      for ( int recID = 0; recID < nrecs; recID++ )
		{
		  pstreamInqRecord(streamID1, &varID, &levelID);
		  if ( lnts1 || (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) )
		    {
		      single = vars[nts][varID][levelID].ptr;
		      pstreamReadRecord(streamID1, single, &nmiss);
		      vars[nts][varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}
       
      tsID++;
    }

  if ( streamID2 != -1 ) pstreamClose(streamID2);
  pstreamClose(streamID1);
 
  if ( operatorID == SELSMON )
    if ( its2 < nts2 )
      cdoWarning("%d timesteps missing after the last month!", nts2-its2);

  if ( array ) Free(array);

  for ( int isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == false )
	{
	 
	  int isel2;
	  int lcont = FALSE;
	  for ( isel2 = isel+1; isel2 < nsel; isel2++ )
	    if ( selfound[isel2] ) break;
	  if ( isel2 == nsel && (nsel-isel) > 1 ) lcont = TRUE;

	  if ( operatorID == SELTIMESTEP )
	    {
	      int lcont2 = FALSE;
	      if ( lcont )
		{    
		  for ( isel2 = isel+1; isel2 < nsel; isel2++ )
		    if ( intarr[isel2-1] != intarr[isel2]-1 ) break;
		  if ( isel2 == nsel ) lcont2 = TRUE;
		}

	      if ( lcont2 )
		{
		  cdoWarning("Time steps %d-%d not found!", intarr[isel], intarr[nsel-1]);
		  break;
		}
	      else
		cdoWarning("Time step %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDATE )
	    {
	      if ( isel == 0 )
		cdoWarning("Date between %14.4f and %14.4f not found!", fltarr[0], fltarr[nsel-1]);
	    }
	  else if ( operatorID == SELTIME )
	    {
	      cdoWarning("Time %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELHOUR )
	    {
	      cdoWarning("Hour %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDAY )
	    {
	      cdoWarning("Day %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELMONTH )
	    {
	      cdoWarning("Month %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELYEAR )
	    {
	      cdoWarning("Year %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELSEASON )
	    {
	      if ( isel < 3 )
		cdoWarning("Month %d not found!", intarr[isel]);
	    }
	}
    }

  if ( selfound ) Free(selfound);

  lista_destroy(ilista);
  lista_destroy(flista);

  if ( lnts1 || nconst )
    {
      for ( tsID = 0; tsID < nts1; tsID++ ) field_free(vars[tsID], vlistID2);

      if ( vars  ) Free(vars);
      if ( vdate_list ) Free(vdate_list);
      if ( vtime_list ) Free(vtime_list);
    }

  vlistDestroy(vlistID2);

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return 0;
}
