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

      Setpartab  setpartab       Set parameter table
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"
#include "pmlist.h"
#include "convert_units.h"

int stringToParam(const char *paramstr);

void paramToString(int param, char *paramstr, int maxlen)
{
  int len;

  int dis, cat, num;
  cdiDecodeParam(param, &num, &cat, &dis);

  size_t umaxlen = maxlen >= 0 ? (unsigned)maxlen : 0U;
  if ( dis == 255 && (cat == 255 || cat == 0 ) )
    len = snprintf(paramstr, umaxlen, "%d", num);
  else  if ( dis == 255 )
    len = snprintf(paramstr, umaxlen, "%d.%d", num, cat);
  else
    len = snprintf(paramstr, umaxlen, "%d.%d.%d", num, cat, dis);

  if ( len >= maxlen || len < 0)
    fprintf(stderr, "Internal problem (%s): size of input string is too small!\n", __func__);
}

typedef enum {CODE_NUMBER, PARAMETER_ID, VARIABLE_NAME, STANDARD_NAME} pt_mode_t;

typedef struct
{
  bool convert;
  bool remove;
  // missing value
  bool changemissval;
  double missval_old;
  //
  bool lfactor;
  double factor;
  //
  bool checkvalid;
  double valid_min;
  double valid_max;
  //
  bool check_min_mean_abs;
  double ok_min_mean_abs;
  //
  bool check_max_mean_abs;
  double ok_max_mean_abs;
  // units
  bool changeunits;
  char units_old[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  // varname
  char name[CDI_MAX_NAME];
  // converter
  void *ut_converter;

  double amean;
  long nvals, n_lower_min, n_greater_max;
} var_t;


void cdo_define_var_units(var_t *var, int vlistID2, int varID, const char *units);

void cmor_check_init(int nvars, var_t *vars);
void cmor_check_eval(int vlistID, int nvars, var_t *vars);
void cmor_check_prep(var_t *var, long gridsize, double missval, double *array);


static
void apply_parameterlist(pt_mode_t ptmode, list_t *pmlist, int nvars, int vlistID2, var_t *vars)
{
  const char *hentry[] = {"Header"};
  const char *ventry[] = {"variable_entry", "parameter"};
  int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);
  int nhentry = (int) sizeof(hentry)/sizeof(hentry[0]);
  char valstr[CDI_MAX_NAME];
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  int codenum;

  // search for global missing value
  bool lmissval = false;
  double missval;
  list_t *kvlist = pmlist_get_kvlist_ventry(pmlist, nhentry, hentry);
  if ( kvlist )
    {
      keyValues_t *kv = kvlist_search(kvlist, "missing_value");
      if ( kv && kv->nvalues > 0 )
        {
          lmissval = true;
          missval = parameter2double(kv->values[0]);
        }
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      var_t *var = &vars[varID];
      vlistInqVarName(vlistID2, varID, varname);

      strcpy(var->name, varname);
      if ( lmissval )
        {
          double missval_old = vlistInqVarMissval(vlistID2, varID);
          if ( ! DBL_IS_EQUAL(missval, missval_old) )
            {
              var->changemissval = true;
              var->missval_old = missval_old;
              vlistDefVarMissval(vlistID2, varID, missval);
            }
        }

      list_t *kvlist = NULL;
      keyValues_t *kv;
      if ( ptmode == CODE_NUMBER )
        {
          codenum = vlistInqVarCode(vlistID2, varID);
          snprintf(valstr, sizeof(valstr), "%d", codenum);
          kvlist = pmlist_search_kvlist_ventry(pmlist, "code", valstr, nventry, ventry);
          if ( kvlist )
            {
              int tableID = vlistInqVarTable(vlistID2, varID);
              int tabnum  = tableInqNum(tableID);
              int levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
              kv = kvlist_search(kvlist, "table");
              int table = (kv && kv->nvalues == 1) ? parameter2int(kv->values[0]) : tabnum;
              kv = kvlist_search(kvlist, "ltype");
              int ltype = (kv && kv->nvalues == 1) ? parameter2int(kv->values[0]) : levtype;
              if ( !(tabnum == table && levtype == ltype) ) kvlist = NULL;
            }
        }
      else if ( ptmode == PARAMETER_ID )
        {
          int param = vlistInqVarParam(vlistID2, varID);
          paramToString(param, paramstr, sizeof(paramstr));
          snprintf(valstr, sizeof(valstr), "%s", paramstr);
          kvlist = pmlist_search_kvlist_ventry(pmlist, "param", valstr, nventry, ventry);
          if ( kvlist )
            {
              int levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
              kv = kvlist_search(kvlist, "ltype");
              int ltype = (kv && kv->nvalues == 1) ? parameter2int(kv->values[0]) : levtype;
              if ( !(levtype == ltype) ) kvlist = NULL;
            }  
        }
      else if ( ptmode == VARIABLE_NAME )
        {
          kvlist = pmlist_search_kvlist_ventry(pmlist, "name", varname, nventry, ventry);
        }

      if ( kvlist )
        {
          int pnum, ptab, pdum;
          cdiDecodeParam(vlistInqVarParam(vlistID2, varID), &pnum, &ptab, &pdum);
          bool lvalid_min = false, lvalid_max = false;

          for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
            {
              kv = *(keyValues_t **)kvnode->data;
              const char *key = kv->key;
              const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;
              bool lv1 = (kv->nvalues == 1);
              
              // printf("key=%s  value=%s\n", key, value ? value : "");

              if      ( lv1 && STR_IS_EQ(key, "standard_name") ) vlistDefVarStdname(vlistID2, varID, value);
              else if ( lv1 && STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID2, varID, value);
              else if ( lv1 && STR_IS_EQ(key, "units")         ) cdo_define_var_units(var, vlistID2, varID, value);
              else if ( lv1 && STR_IS_EQ(key, "name") && ptmode != VARIABLE_NAME ) vlistDefVarName(vlistID2, varID, parameter2word(value));
              else if ( lv1 && STR_IS_EQ(key, "out_name")      )
                {
                  const char *outname = parameter2word(value);
                  if ( !STR_IS_EQ(var->name, outname) )
                    {
                      vlistDefVarName(vlistID2, varID, outname);
                      cdiDefAttTxt(vlistID2, varID, "original_name", (int)strlen(var->name), var->name);
                    }
                }
              else if ( lv1 && STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( lv1 && STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( lv1 && STR_IS_EQ(key, "code")          ) vlistDefVarParam(vlistID2, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              else if ( lv1 && STR_IS_EQ(key, "out_code")      ) vlistDefVarParam(vlistID2, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              else if ( lv1 && STR_IS_EQ(key, "comment")       ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( lv1 && STR_IS_EQ(key, "chunktype")     ) ;
              else if ( lv1 && STR_IS_EQ(key, "cell_methods")  ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( lv1 && STR_IS_EQ(key, "cell_measures") ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( lv1 && STR_IS_EQ(key, "delete")        ) var->remove = parameter2bool(value);
              else if ( lv1 && STR_IS_EQ(key, "convert")       ) var->convert = parameter2bool(value);
              else if ( lv1 && STR_IS_EQ(key, "factor")        )
                {
                  var->lfactor = true;
                  var->factor = parameter2double(value);
                  if ( cdoVerbose ) cdoPrint("%s - scale factor %g", varname, var->factor);
                }
              else if ( lv1 && (STR_IS_EQ(key, "missval") || STR_IS_EQ(key, "missing_value")) )
                {
                  double missval = parameter2double(value);
                  double missval_old = vlistInqVarMissval(vlistID2, varID);
                  if ( ! DBL_IS_EQUAL(missval, missval_old) )
                    {
                      if ( cdoVerbose ) cdoPrint("%s - change missval from %g to %g", varname, missval_old, missval);
                      var->changemissval = true;
                      var->missval_old = missval_old;
                      vlistDefVarMissval(vlistID2, varID, missval);
                    }
                }
              else if ( lv1 && STR_IS_EQ(key, "valid_min") )
                {
                  lvalid_min = true;
                  var->valid_min = parameter2double(value);
                }
              else if ( lv1 && STR_IS_EQ(key, "valid_max") )
                {
                  lvalid_max = true;
                  var->valid_max = parameter2double(value);
                }
              else if ( lv1 && STR_IS_EQ(key, "ok_min_mean_abs") )
                {
                  var->check_min_mean_abs = true;
                  var->ok_min_mean_abs = parameter2double(value);
                }
              else if ( lv1 && STR_IS_EQ(key, "ok_max_mean_abs") )
                {
                  var->check_max_mean_abs = true;
                  var->ok_max_mean_abs = parameter2double(value);
                }
              else if ( lv1 && (STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type")) )
                {
                  int datatype = str2datatype(parameter2word(value));
                  if ( datatype != -1 ) vlistDefVarDatatype(vlistID2, varID, datatype);
                }
              else if ( lv1 && STR_IS_EQ(key, "dimensions") )
                {
                }
              else
                {
                  int nvalues = kv->nvalues;
                  if ( nvalues == 1 && !*value ) nvalues = 0;
                  int dtype = literals_find_datatype(nvalues, kv->values);
                  
                  if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
                    {
                      int *ivals = (int*) Malloc(nvalues*sizeof(int));
                      for ( int i = 0; i < nvalues; ++i ) ivals[i] = literal_to_int(kv->values[i]);
                      cdiDefAttInt(vlistID2, varID, key, dtype, nvalues, ivals);
                      Free(ivals);
                    }
                  else if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
                    {
                      double *dvals = (double*) Malloc(nvalues*sizeof(double));
                      for ( int i = 0; i < nvalues; ++i ) dvals[i] = literal_to_double(kv->values[i]);
                      cdiDefAttFlt(vlistID2, varID, key, dtype, nvalues, dvals);
                      Free(dvals);
                    }
                  else
                    {
                      int len = (value && *value) ? (int) strlen(value) : 0;
                      cdiDefAttTxt(vlistID2, varID, key, len, value);
                    }
                }
            }

          if ( lvalid_min && lvalid_max ) var->checkvalid = true;
        }
      else
        {
          if      ( ptmode == CODE_NUMBER )   cdoPrint("Code number %d not found in parameter table!", codenum);
          else if ( ptmode == PARAMETER_ID )  cdoPrint("Parameter ID %s not found in parameter table!", paramstr);
          else if ( ptmode == VARIABLE_NAME ) cdoPrint("Variable %s not found in parameter table!", varname);
        }
    }
}


void *Setpartab(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int tableID = -1;
  int tableformat = 0;
  bool delvars = false;
  double missval;

  cdoInitialize(argument);

  int SETCODETAB = cdoOperatorAdd("setcodetab",  0, 0, "parameter code table name");
  int SETPARTABC = cdoOperatorAdd("setpartabc",  0, 0, "parameter table name");
  int SETPARTABP = cdoOperatorAdd("setpartabp",  0, 0, "parameter table name");
  int SETPARTABN = cdoOperatorAdd("setpartabn",  0, 0, "parameter table name");

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");

  bool convert_data = false;
  if ( operatorArgc() == 2 )
    {
      if ( strcmp("convert", operatorArgv()[1]) == 0 ) convert_data = true;
      else cdoAbort("Unknown parameter: >%s<", operatorArgv()[1]); 
    }

  if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");

  pt_mode_t ptmode = CODE_NUMBER;
  if      ( operatorID == SETCODETAB ) ptmode = CODE_NUMBER;
  else if ( operatorID == SETPARTABC ) ptmode = CODE_NUMBER;
  else if ( operatorID == SETPARTABP ) ptmode = PARAMETER_ID;
  else if ( operatorID == SETPARTABN ) ptmode = VARIABLE_NAME;

  if ( ptmode == CODE_NUMBER )
    {
      char *partab = operatorArgv()[0];
      FILE *fp = NULL;
      if ( fileExists(partab) ) fp = fopen(partab, "r");
      if ( fp != NULL )
	{
	  fseek(fp, 0L, SEEK_END);
	  size_t fsize = (size_t) ftell(fp);
	  char *parbuf = (char *) Malloc(fsize+1);
	  fseek(fp, 0L, SEEK_SET);
	  fread(parbuf, fsize, 1, fp);
	  parbuf[fsize] = 0;
	  fseek(fp, 0L, SEEK_SET);

	  if ( atoi(parbuf) == 0 ) tableformat = 1;

	  fclose(fp);
	  Free(parbuf);
	}

      if ( tableformat == 0 ) tableID = defineTable(partab);
    }
  else if (  ptmode == PARAMETER_ID )
    {
      tableformat = 1;
    }
  else if (  ptmode == VARIABLE_NAME )
    {
      tableformat = 1;
    }

  if ( cdoVerbose ) cdoPrint("Table format version %d", tableformat);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  // vlistPrint(vlistID2);

  int nvars = vlistNvars(vlistID2);
  var_t *vars = (var_t *) Malloc(nvars*sizeof(var_t));
  memset(vars, 0, nvars*sizeof(var_t));

  if ( convert_data )
    for ( varID = 0; varID < nvars; ++varID ) vars[varID].convert = true;

  if ( tableformat == 0 )
    {
      /*
        for ( int varID = 0; varID < nvars; varID++ )
	  vlistDefVarTable(vlistID2, varID, tableID);
      */
      char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
        {
          int param = vlistInqVarParam(vlistID2, varID);
          int pdis, pcat, pnum;
          cdiDecodeParam(param, &pnum, &pcat, &pdis);
          if ( pdis == 255 )
            {
              int code = pnum;
              int ltype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
              name[0] = 0;
              longname[0] = 0;
              units[0] = 0;
              tableInqEntry(tableID, code, ltype, name, longname, units);
              if ( name[0] )
                {
                  vlistDefVarName(vlistID2, varID, name);
                  if ( longname[0] ) vlistDefVarLongname(vlistID2, varID, longname);
                  if ( units[0] )    vlistDefVarUnits(vlistID2, varID, units);
                }
            }
          vlistDefVarTable(vlistID2, varID, tableID);
        }
    }
  else
    {
      const char *filename = operatorArgv()[0];
      FILE *fp = fopen(filename, "r");
      if ( fp == NULL ) cdoAbort("Open failed on: %s\n", filename);
      
      list_t *pmlist = namelist_to_pmlist(fp, filename);
      fclose(fp);

      apply_parameterlist(ptmode, pmlist, nvars, vlistID2, vars);
      list_destroy(pmlist);

      for ( int varID = 0; varID < nvars; ++varID )
	if ( vars[varID].remove )
          {
            delvars = true;
            break;
          }

      if ( delvars )
	{
	  vlistClearFlag(vlistID1);
	  vlistClearFlag(vlistID2);

	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID2, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		{
		  vlistDefFlag(vlistID1, varID, levID, TRUE);
		  vlistDefFlag(vlistID2, varID, levID, TRUE);
		  if ( vars[varID].remove )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      vlistDefFlag(vlistID2, varID, levID, FALSE);
		    }
		}
	    }

	  int vlistIDx = vlistCreate();
	  cdoVlistCopyFlag(vlistIDx, vlistID2);

	  vlistDestroy(vlistID2);

	  vlistID2 = vlistIDx;
          if ( vlistNvars(vlistID2) == 0 ) cdoAbort("No variable selected!");
	}

      for ( int varID = 0; varID < nvars; ++varID )
        {
          var_t *var = &vars[varID];
          if ( var->convert == false ) var->changeunits = false;
          if ( var->changeunits )
            cdoConvertUnits(&var->ut_converter, &var->changeunits, (char*)&var->units, (char*)&var->units_old, var->name);
        }
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // vlistPrint(vlistID2);
  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  long gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double *) Malloc(gridsize*sizeof(double));

  int tsID1 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID1);
	       
      cmor_check_init(nvars, vars);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

          var_t *var = &vars[varID];
	  int varID2 = varID;
	  int levelID2 = levelID;

	  if ( delvars )
	    {
	      if ( var->remove ) continue;

	      if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
		{
		  varID2   = vlistFindVar(vlistID2, varID);
		  levelID2 = vlistFindLevel(vlistID2, varID, levelID);
		}
	    }

	  pstreamDefRecord(streamID2,  varID2,  levelID2);

	  pstreamReadRecord(streamID1, array, &nmiss);

	  missval = vlistInqVarMissval(vlistID2, varID2);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  if ( vlistInqVarNumber(vlistID2, varID2) != CDI_REAL ) gridsize *= 2;

	  if ( nmiss > 0 && var->changemissval )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( DBL_IS_EQUAL(array[i], var->missval_old) ) array[i] = missval;
		}
	    }

	  if ( var->lfactor )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) ) array[i] *= var->factor;
		}
	    }

#if defined(HAVE_UDUNITS2)
	  if ( var->changeunits )
	    {
	      int nerr = 0;
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    {
		      array[i] = cv_convert_double((const cv_converter*)var->ut_converter, array[i]);
		      if ( ut_get_status() != UT_SUCCESS ) nerr++;
		    }
		}
	      if ( nerr )
		{
		  cdoWarning("Udunits: Error converting units from [%s] to [%s], parameter: %s",
			     var->units_old, var->units, var->name);
		  var->changeunits = false;
		}
	    }
#endif
	  
	  pstreamWriteRecord(streamID2, array, nmiss);

          cmor_check_prep(var, gridsize, missval, array);
	}

      cmor_check_eval(vlistID2, nvars, vars);

      tsID1++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

#if defined(HAVE_UDUNITS2)
  for ( int varID = 0; varID < nvars; varID++ )
    if ( vars[varID].ut_converter ) cdoConvertFree(vars[varID].ut_converter);

  cdoConvertDestroy();
#endif

  if ( array ) Free(array);
  if ( vars  ) Free(vars);

  cdoFinish();

  return 0;
}
