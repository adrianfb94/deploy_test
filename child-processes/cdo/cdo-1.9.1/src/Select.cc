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

      Select      select         Select fields
*/

#include <cdi.h>
#include "cdo_int.h"
#include "pstream.h"
#include "sellist.h"


double datestr_to_double(const char *datestr, int opt);

bool *cdo_read_timestepmask(const char *maskfile, int *n);

static
void write_const_vars(int streamID2, int vlistID2, int nvars, double **vardata2)
{
  for ( int varID2c = 0; varID2c < nvars; ++varID2c )
    {
      if ( vardata2[varID2c] )
        {
          double missval = vlistInqVarMissval(vlistID2, varID2c);
          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2c));
          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2c));
          for ( int levelID2c = 0; levelID2c < nlevel; ++levelID2c )
            {
              double *pdata = vardata2[varID2c]+gridsize*levelID2c;
              int nmiss = 0;
              for ( int i = 0; i < gridsize; ++i )
                if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

              // if ( levelID2c == 0 ) printf("Write varID %d\n", varID2c);
              pstreamDefRecord(streamID2, varID2c, levelID2c);
              pstreamWriteRecord(streamID2, pdata, nmiss);
            }
          Free(vardata2[varID2c]);
          vardata2[varID2c] = NULL;
        }
    }
}

static
void eval_timestepmask(const char *maskfile, list_t *kvlist)
{
  int n = 0;
  bool *imask = cdo_read_timestepmask(maskfile, &n);

  int nvals = 0;
  for ( int i = 0; i < n; ++i ) if ( imask[i] ) nvals++;
  if ( nvals == 0 ) cdoPrint("timestepmask has no values!");
  else
    {
      char **values = (char**) Malloc(nvals*sizeof(char*));
      int j = 0;
      for ( int i = 0; i < n; ++i )
        {
          if ( imask[i] )
            {
              size_t length = (size_t)log10(j+1)+2;
              values[j] = (char*) Malloc(length*sizeof(char));
              sprintf(values[j++], "%d", i+1);
            }
        }

      kvlist_append(kvlist, "timestep", (const char **)values, nvals);

      for ( int i = 0; i < nvals; ++i ) Free(values[i]);
      Free(values);
    }
      
  Free(imask);
}


void *Select(void *argument)
{
  bool lconstvars = true;
  int streamID2 = CDI_UNDEFID;
  int nrecs;
  int nvars, nvars2;
  int varID, levelID;
  int last_year = -999999999;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gname[CDI_MAX_NAME];
  char zname[CDI_MAX_NAME];
  int vlistID0 = -1, vlistID2 = -1;
  int taxisID2 = CDI_UNDEFID;
  int ntsteps2 = 0;
  bool ltimsel = false;
  bool *vars = NULL;
  double **vardata2 = NULL;
  double *array = NULL;
  double fstartdate = -99999999999.;
  double fenddate   = -99999999999.;

  cdoInitialize(argument);

  int SELECT = cdoOperatorAdd("select", 0, 0, "parameter list");
  int DELETE = cdoOperatorAdd("delete", 0, 0, "parameter list");

  bool lcopy = UNCHANGED_RECORD;

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nsel = operatorArgc();
  char **argnames = operatorArgv();

  if ( nsel == 0 ) cdoAbort("Parameter missing!");

  list_t *kvlist = kvlist_new("SELECT");
  if ( kvlist_parse_cmdline(kvlist, nsel, argnames) != 0 ) cdoAbort("Parse error!");
  if ( cdoVerbose ) kvlist_print(kvlist);

  keyValues_t *kv = kvlist_search(kvlist, "timestepmask");
  if ( kv && kv->nvalues > 0 )
    {
      if ( kvlist_search(kvlist, "timestep") ) cdoAbort("Parameter timestep and timestepmask can't be combined!");
      eval_timestepmask(kv->values[0], kvlist);
    }

  sellist_t *sellist = sellist_create(kvlist);

  // clang-format off
  SELLIST_ADD_INT(timestep_of_year, "Timestep of year");
  SELLIST_ADD_INT(timestep,         "Timestep");
  SELLIST_ADD_INT(year,             "Year");
  SELLIST_ADD_INT(month,            "Month");
  SELLIST_ADD_INT(day,              "Day");
  SELLIST_ADD_INT(hour,             "Hour");
  SELLIST_ADD_INT(minute,           "Minute");
  SELLIST_ADD_INT(code,             "Code number");
  SELLIST_ADD_INT(levidx,           "Level index");
  SELLIST_ADD_INT(ltype,            "Level type");
  SELLIST_ADD_INT(zaxisnum,         "Zaxis number");
  SELLIST_ADD_INT(gridnum,          "Grid number");
  SELLIST_ADD_FLT(level,            "Level");
  SELLIST_ADD_WORD(name,            "Variable name");
  SELLIST_ADD_WORD(param,           "Parameter");
  SELLIST_ADD_WORD(zaxisname,       "Zaxis name");
  SELLIST_ADD_WORD(gridname,        "Grid name");
  SELLIST_ADD_WORD(steptype,        "Time step type");
  SELLIST_ADD_WORD(startdate,       "Start date");
  SELLIST_ADD_WORD(enddate,         "End date");
  SELLIST_ADD_WORD(season,          "Season");
  SELLIST_ADD_WORD(date,            "Date");
  SELLIST_ADD_WORD(timestepmask,    "Timestep mask");
  // clang-format on

  if ( cdoVerbose ) sellist_print(sellist);

  sellist_verify(sellist);

  if ( SELLIST_NVAL(timestepmask) > 1 ) cdoAbort("Key timestepmask has too many values!");
  UNUSED(timestepmask);

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  dtlist_type *dtlist = dtlist_new();

  if ( !cdoVerbose && nfiles > 1 ) progressInit();

  timestep = 0;
  int tsID2 = 0;
  for ( int indf = 0; indf < nfiles; ++indf )
    {
      if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, (indf+1.)/nfiles);
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      int streamID1 = pstreamOpenRead(cdoStreamName(indf));

      int vlistID1 = pstreamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      bool lcopy_const = false;

      if ( indf == 0 )
	{
          bool xresult = true;

	  // vlistID0 = vlistDuplicate(vlistID1);

	  vlistClearFlag(vlistID1);
	  nvars = vlistNvars(vlistID1);
	  vars  = (bool*) Malloc(nvars*sizeof(bool));

	  if ( operatorID == DELETE )
	    {
	      xresult = false;
	      for ( varID = 0; varID < nvars; ++varID )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  int nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; ++levID )
		    vlistDefFlag(vlistID1, varID, levID, TRUE);
		}
	    }

          bool lvarsel = SELLIST_NVAL(code) || SELLIST_NVAL(ltype) || SELLIST_NVAL(zaxisnum) ||
            SELLIST_NVAL(gridnum) || SELLIST_NVAL(name) || SELLIST_NVAL(param) ||
            SELLIST_NVAL(zaxisname) || SELLIST_NVAL(gridname) || SELLIST_NVAL(steptype);

          bool llevsel = SELLIST_NVAL(level) || SELLIST_NVAL(levidx);

	  ltimsel = SELLIST_NVAL(date) || SELLIST_NVAL(startdate) || SELLIST_NVAL(enddate) || SELLIST_NVAL(season) ||
            SELLIST_NVAL(timestep_of_year) || SELLIST_NVAL(timestep) || SELLIST_NVAL(year) || SELLIST_NVAL(month) ||
            SELLIST_NVAL(day) || SELLIST_NVAL(hour) || SELLIST_NVAL(minute);
          
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int iparam = vlistInqVarParam(vlistID1, varID);
	      code = vlistInqVarCode(vlistID1, varID);
	      vlistInqVarName(vlistID1, varID, varname);
	      vlistInqVarStdname(vlistID1, varID, stdname);

	      cdiParamToString(iparam, paramstr, sizeof(paramstr));

	      name  = varname;
	      param = paramstr;

	      int gridID  = vlistInqVarGrid(vlistID1, varID);
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      ltype = zaxis2ltype(zaxisID);

              zaxisnum = vlistZaxisIndex(vlistID1, zaxisID)+1;
              int zaxistype = zaxisInqType(zaxisID);
              zaxisName(zaxistype, zname);
              zaxisname = zname;

              gridnum = vlistGridIndex(vlistID1, gridID)+1;
              int gridtype = gridInqType(gridID);
              gridName(gridtype, gname);
              gridname = gname;

              int tsteptype = vlistInqVarTsteptype(vlistID1, varID);
              if      ( tsteptype == TSTEP_INSTANT  ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT2 ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT3 ) steptype = "instant";
              else if ( tsteptype == TSTEP_MIN      ) steptype = "min";
              else if ( tsteptype == TSTEP_MAX      ) steptype = "max";
              else if ( tsteptype == TSTEP_AVG      ) steptype = "avg";
              else if ( tsteptype == TSTEP_ACCUM    ) steptype = "accum";
              else if ( tsteptype == TSTEP_RANGE    ) steptype = "range";
              else if ( tsteptype == TSTEP_DIFF     ) steptype = "diff";
              else                                    steptype = "unknown";
              
	      vars[varID] = false;

              bool found_code  = SELLIST_CHECK(code);
              bool found_name  = SELLIST_CHECK(name);
              bool found_param = SELLIST_CHECK(param);
              bool found_grid  = SELLIST_CHECK(gridnum);
              bool found_gname = SELLIST_CHECK(gridname);
              bool found_stype = SELLIST_CHECK(steptype);
              bool found_ltype = SELLIST_CHECK(ltype);
              bool found_zaxis = SELLIST_CHECK(zaxisnum);
              bool found_zname = SELLIST_CHECK(zaxisname);

              bool lvar  = found_code || found_name || found_param;
              bool lstep = SELLIST_NVAL(steptype) ? found_stype : true;
              bool lgrid = (SELLIST_NVAL(gridnum) || SELLIST_NVAL(gridname)) ? (found_grid || found_gname) : true;
              bool lvert = (SELLIST_NVAL(ltype) || SELLIST_NVAL(zaxisnum) || SELLIST_NVAL(zaxisname)) ? (found_ltype || found_zaxis || found_zname) : true;
	     
              if ( !vars[varID] && lgrid && lvar ) vars[varID] = true;
              if ( !vars[varID] && lvert && lvar ) vars[varID] = true;
              if ( !vars[varID] && lstep && lvar ) vars[varID] = true;

              if ( !vars[varID] && !lvar )
                {
                  if      ( found_grid || found_gname ) vars[varID] = true;
                  else if ( found_stype ) vars[varID] = true;
                  else if ( found_ltype || found_zaxis || found_zname ) vars[varID] = true;
                  else if ( !lvarsel && (SELLIST_NVAL(levidx) || SELLIST_NVAL(level)) )
                    {
                      for ( int levID = 0; levID < nlevs; ++levID )
                        {
                          levidx = levID + 1;
                          level = cdoZaxisInqLevel(zaxisID, levID);
                          if ( !vars[varID] && SELLIST_CHECK(levidx) ) vars[varID] = true;
                          if ( !vars[varID] && SELLIST_CHECK(level)  ) vars[varID] = true;
                        }
                    }
                }
	    }

	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      if ( vars[varID] )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                  if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
                    {
                      int psvarid = vlist_get_psvarid(vlistID1, zaxisID);
                      if ( psvarid != -1 && !vars[psvarid] ) vars[psvarid] = true;
                    }
                }
            }

	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      if ( vars[varID] )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  int nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; ++levID )
		    {
		      levidx = levID + 1;
                      level = cdoZaxisInqLevel(zaxisID, levID);
		      
		      if ( nlevs == 1 && IS_EQUAL(level, 0) )
			{
                          SELLIST_CHECK(level);
			  vlistDefFlag(vlistID1, varID, levID, xresult);
			}
		      else
			{
			  if ( SELLIST_NVAL(levidx) )
			    {
			      if ( SELLIST_CHECK(levidx) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			  else if ( SELLIST_NVAL(level) )
			    {
			      if ( SELLIST_CHECK(level) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			  else
			    {
			      vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			}
		    }
		}
	    }

	  SELLIST_CHECK_FLAG(code);
	  SELLIST_CHECK_FLAG(levidx);
	  SELLIST_CHECK_FLAG(ltype);
	  SELLIST_CHECK_FLAG(zaxisnum);
	  SELLIST_CHECK_FLAG(gridnum);
	  SELLIST_CHECK_FLAG(level);
	  SELLIST_CHECK_FLAG(name);
	  SELLIST_CHECK_FLAG(param);
	  SELLIST_CHECK_FLAG(zaxisname);
	  SELLIST_CHECK_FLAG(gridname);
	  SELLIST_CHECK_FLAG(steptype);

	  int npar = 0;
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; ++levID )
		if ( vlistInqFlag(vlistID1, varID, levID) == TRUE )
                  {
                    npar++;
                    break;
                  }
	    }

	  if ( npar == 0 )
	    {
	      if ( (! lvarsel) && (! llevsel) && ltimsel )
		{
                  lcopy_const = true;

		  for ( varID = 0; varID < nvars; ++varID )
		    {
		      vars[varID] = true;
		      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		      int nlevs   = zaxisInqSize(zaxisID);
		      for ( int levID = 0; levID < nlevs; ++levID )
			vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	      else
		{
		  cdoAbort("No variable selected!");
		}
	    }
          else
            {
              if ( lvarsel && ltimsel )
                {
                  for ( varID = 0; varID < nvars; ++varID )
                    {
                      if ( vars[varID] == true && vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
                        {
                          lcopy_const = true;
                          break;
                        }
                    }
                }
            }

	  //if ( cdoVerbose ) vlistPrint(vlistID1);

	  vlistID0 = vlistDuplicate(vlistID1);
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; ++levID )
		vlistDefFlag(vlistID0, varID, levID, vlistInqFlag(vlistID1, varID, levID));
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID0);

	  vlistID2 = vlistCreate();
	  cdoVlistCopyFlag(vlistID2, vlistID0);

	  //if ( cdoVerbose ) vlistPrint(vlistID2);

	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  int ntsteps = vlistNtsteps(vlistID1);

	  nvars2 = vlistNvars(vlistID2);

	  if ( ntsteps == 1 && nfiles == 1 )
	    {
	      for ( varID = 0; varID < nvars2; ++varID )
		if ( vlistInqVarTimetype(vlistID2, varID) != TIME_CONSTANT ) break;

	      if ( varID == nvars2 ) ntsteps = 0;
	    }

	  ntsteps2 = (operatorID == SELECT && SELLIST_NVAL(timestep) == 1) ? 1 : ntsteps;

	  if ( ntsteps2 == 0 && nfiles > 1 )
	    {
              lconstvars = false;
	      for ( varID = 0; varID < nvars2; ++varID )
		vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
	    }

	  // support for negative timestep values
	  if ( SELLIST_NVAL(timestep) > 0 && ntsteps > 0 && nfiles == 1 )
	    {
	      for ( int i = 0; i < SELLIST_NVAL(timestep); ++i )
		{
                  int ptimestep;
                  SELLIST_GET_VAL(timestep, i, &ptimestep);
		  if ( ptimestep < 0 )
		    {
		      if ( cdoVerbose )
			cdoPrint("timestep %d changed to %d", ptimestep, ptimestep + ntsteps + 1);
		      ptimestep += ntsteps + 1;
                      SELLIST_DEF_VAL(timestep, i, &ptimestep);
		    }
		}
	    }

	  if ( ! lcopy )
	    {
	      int gridsize = vlistGridsizeMax(vlistID1);
	      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
	      array = (double*) Malloc(gridsize*sizeof(double));
	    }

	  SELLIST_GET_VAL(startdate, 0, &startdate);
	  SELLIST_GET_VAL(enddate, 0, &enddate);
	  if ( SELLIST_NVAL(startdate) ) fstartdate = datestr_to_double(startdate, 0);
	  if ( SELLIST_NVAL(enddate)   ) fenddate   = datestr_to_double(enddate, 1);
	}
      else
	{
	  vlistCompare(vlistID0, vlistID1, CMP_ALL);
	}


      if ( nvars2 == 0 )
	{
	  cdoWarning("No variable selected!");
	  goto END_LABEL;
	}

      if ( lcopy_const )
        {
          vardata2 = (double**) Malloc(nvars2*sizeof(double));
          for ( varID = 0; varID < nvars2; ++varID ) vardata2[varID] = NULL;
        }

      bool lstop = false;
      int tsID1 = 0;
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
	{
          timestep++;
	  bool copytimestep = true;

	  if ( ltimsel == true )
	    {
	      copytimestep = false;

	      if ( operatorID == SELECT && SELLIST_NVAL(timestep) > 0 )
		{
                  int ptimestep;
                  SELLIST_GET_VAL(timestep, SELLIST_NVAL(timestep)-1, &ptimestep);
                  if ( timestep > ptimestep )
                    {
                      lstop = true;
                      break;
                    }
                }

              dtlist_taxisInqTimestep(dtlist, taxisID1, 0);
              int vdate = dtlist_get_vdate(dtlist, 0);
              int vtime = dtlist_get_vtime(dtlist, 0);
              int second;
	      cdiDecodeDate(vdate, &year, &month, &day);
	      cdiDecodeTime(vtime, &hour, &minute, &second);
              UNUSED(season);

	      if ( year != last_year )
		{
		  timestep_of_year = 0;
		  last_year = year;
		}

	      timestep_of_year++;

	      if ( SELLIST_CHECK(timestep) ) copytimestep = true;
	      if ( SELLIST_CHECK(timestep_of_year) ) copytimestep = true;

	      if ( !copytimestep && SELLIST_NVAL(date) == 0 && SELLIST_NVAL(timestep) == 0 && SELLIST_NVAL(timestep_of_year) == 0 )
		{
		  bool lseason = false, lyear = false, lmonth = false, lday = false, lhour = false, lminute = false;

		  if ( SELLIST_NVAL(season) == 0 || SELLIST_CHECK_SEASON(season, month) ) lseason   = true;
		  if ( SELLIST_NVAL(year)   == 0 || SELLIST_CHECK(year)   ) lyear   = true;
		  if ( SELLIST_NVAL(month)  == 0 || SELLIST_CHECK(month)  ) lmonth  = true;
		  if ( SELLIST_NVAL(day)    == 0 || SELLIST_CHECK(day)    ) lday    = true;
		  if ( SELLIST_NVAL(hour)   == 0 || SELLIST_CHECK(hour)   ) lhour   = true;
		  if ( SELLIST_NVAL(minute) == 0 || SELLIST_CHECK(minute) ) lminute = true;

		  if ( lseason && lyear && lmonth && lday && lhour && lminute ) copytimestep = true;
		}

	      double fdate = ((double)vdate) + ((double)vtime)/1000000.;

	      if ( SELLIST_NVAL(enddate) )
		{
                  copytimestep = (fdate <= fenddate);
		  if ( fdate > fenddate )
		    {
                      SELLIST_DEF_FLAG(enddate, 0, true);
		      if ( operatorID == SELECT )
			{
			  lstop = true;
			  break;
			}
		    }
		}

	      if ( SELLIST_NVAL(startdate) )
		{
                  copytimestep = (fdate >= fstartdate);
		  if ( fdate >= fstartdate ) SELLIST_DEF_FLAG(startdate, 0, true);
		}

              if ( SELLIST_NVAL(date) )
                {
                  char vdatetimestr[64];
                  datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
                  date = vdatetimestr;
                  if ( SELLIST_CHECK_DATE(date) ) copytimestep = true;
                }

	      if ( operatorID == DELETE ) copytimestep = !copytimestep;

              if ( copytimestep && indf == 0 && tsID1 == 0 ) lcopy_const = false;
	    }

	  if ( copytimestep == true )
	    {
	      if ( streamID2 == CDI_UNDEFID )
		{
                  bool lasttimestep = (nfiles == 1) && (ntsteps2 > 1) && (ntsteps2 == (tsID1+1));
                  if ( lasttimestep && tsID2 == 0 ) ntsteps2 = 1;
                  if ( ntsteps2 == 0 || ntsteps2 == 1 ) vlistDefNtsteps(vlistID2, ntsteps2);
		  streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
		  pstreamDefVlist(streamID2, vlistID2);
		}
	      taxisCopyTimestep(taxisID2, taxisID1);
	      pstreamDefTimestep(streamID2, tsID2);
              
	      for ( int recID = 0; recID < nrecs; ++recID )
		{
		  pstreamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
                      if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                        if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
                          continue;

		      int varID2   = vlistFindVar(vlistID2, varID);
		      int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
		      
                      if ( lcopy_const && tsID2 == 0 )
                        write_const_vars(streamID2, vlistID2, varID2, vardata2);

                      pstreamDefRecord(streamID2, varID2, levelID2);
                      // if ( levelID2 == 0 ) printf("Write varID %d\n", varID2);
		      if ( lcopy )
			{
			  pstreamCopyRecord(streamID2, streamID1);
			}
		      else
			{
                          int nmiss;
			  pstreamReadRecord(streamID1, array, &nmiss);
			  pstreamWriteRecord(streamID2, array, nmiss);
			}
		    }
		}

              if ( lcopy_const && tsID2 == 0 )
                write_const_vars(streamID2, vlistID2, nvars2, vardata2);

	      tsID2++;              
	    }
          else if ( lcopy_const && indf == 0 && tsID1 == 0 )
            {
	      for ( int recID = 0; recID < nrecs; ++recID )
		{
		  pstreamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
		      int varID2 = vlistFindVar(vlistID2, varID);
                      if ( vlistInqVarTimetype(vlistID2, varID2) == TIME_CONSTANT )
                        {
                          int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          if ( levelID == 0 )
                            {
                              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                              vardata2[varID2] = (double*) Malloc(gridsize*nlevel*sizeof(double));
                            }
                          int nmiss;
                          pstreamReadRecord(streamID1, vardata2[varID2]+gridsize*levelID2, &nmiss);
                        }
		    }
		}
            }

	  tsID1++;
	}

      pstreamClose(streamID1);

      if ( lstop ) break;
    }

 END_LABEL:

  if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, 1);    

  SELLIST_CHECK_FLAG(timestep_of_year);
  SELLIST_CHECK_FLAG(timestep);
  SELLIST_CHECK_FLAG(year);
  SELLIST_CHECK_FLAG(month);
  SELLIST_CHECK_FLAG(day);
  SELLIST_CHECK_FLAG(hour);
  SELLIST_CHECK_FLAG(minute);
  SELLIST_CHECK_FLAG(startdate);
  //  SELLIST_CHECK_FLAG(enddate);
  SELLIST_CHECK_FLAG(season);
  SELLIST_CHECK_FLAG(date);

  if ( streamID2 != CDI_UNDEFID ) pstreamClose(streamID2);

  dtlist_delete(dtlist);

  vlistDestroy(vlistID0);
  vlistDestroy(vlistID2);

  sellist_destroy(sellist);
  kvlist_destroy(kvlist);

  if ( array ) Free(array);
  if ( vars ) Free(vars);
  if ( vardata2 )
    {
      for ( varID = 0; varID < nvars2; ++varID )
        if ( vardata2[varID] ) Free(vardata2[varID]);

      Free(vardata2);
    }

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return 0;
}
