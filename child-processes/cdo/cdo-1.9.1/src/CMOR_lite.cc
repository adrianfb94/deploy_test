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

#if  defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBUDUNITS2) && (defined(HAVE_UDUNITS2_H) || defined(HAVE_UDUNITS2_UDUNITS2_H))
#define HAVE_UDUNITS2
#endif

#if defined(HAVE_UDUNITS2)
#if defined(HAVE_UDUNITS2_UDUNITS2_H)
#  include <udunits2/udunits2.h>
#else
#  include <udunits2.h>
#endif
#endif

#include <errno.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"
#include "pmlist.h"
#include "convert_units.h"

int stringToParam(const char *paramstr);

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


void cdo_define_var_units(var_t *var, int vlistID2, int varID, const char *units)
{
  char units_old[CDI_MAX_NAME];

  vlistInqVarUnits(vlistID2, varID, units_old);
  size_t len1 = strlen(units_old);
  size_t len2 = strlen(units);

  if ( strcmp(units, units_old) != 0 )
    {
      if ( len1 > 0 && len2 > 0 )
	{
	  var->changeunits = true;
	  strcpy(var->units_old, units_old);
	  strcpy(var->units, units);
	}

      vlistDefVarUnits(vlistID2, varID, units);
      cdiDefAttTxt(vlistID2, varID, "original_units", (int)strlen(units_old), units_old);
    }
}


void cmor_check_init(int nvars, var_t *vars)
{
  for ( int varID = 0; varID < nvars; ++varID )
    {
      var_t *var = &vars[varID];
      if ( var->checkvalid || var->check_min_mean_abs || var->check_max_mean_abs )
        {
          var->amean = 0;
          var->nvals = 0;
          var->n_lower_min = 0;
          var->n_greater_max = 0;
        }
    }
}


void cmor_check_eval(int vlistID, int nvars, var_t *vars)
{
  char varname[CDI_MAX_NAME];

  for ( int varID = 0; varID < nvars; ++varID )
    {
      var_t *var = &vars[varID];
      if ( var->checkvalid || var->check_min_mean_abs || var->check_max_mean_abs )
        {
          double amean = var->amean;
          long nvals = var->nvals;

          if ( nvals > 0 ) amean /= nvals;

          long n_lower_min = var->n_lower_min;
          long n_greater_max = var->n_greater_max;

          vlistInqVarName(vlistID, varID, varname);

          if ( n_lower_min > 0 )
            cdoWarning("Invalid value(s) detected for variable '%s': %i values were lower than minimum valid value (%.4g).",
                       varname, n_lower_min, var->valid_min);
          if ( n_greater_max > 0 )
            cdoWarning("Invalid value(s) detected for variable '%s': %i values were greater than maximum valid value (%.4g).",
                       varname, n_greater_max, var->valid_max);

          if ( var->check_min_mean_abs )
            {
              if ( amean < .1*var->ok_min_mean_abs )
                cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is lower by more than an order of magnitude than minimum allowed: %.4g",
                           varname, amean, var->ok_min_mean_abs);

              if ( amean < var->ok_min_mean_abs)
                cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is lower than minimum allowed: %.4g",
                           varname, amean, var->ok_min_mean_abs);
            }

          if ( var->check_max_mean_abs )
            {
              if ( amean > 10.*var->ok_max_mean_abs )
                cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is greater by more than an order of magnitude than maximum allowed: %.4g",
                           varname, amean, var->ok_max_mean_abs);
      
              if ( amean > var->ok_max_mean_abs )
                cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is greater than maximum allowed: %.4g",
                           varname, amean, var->ok_max_mean_abs);
            }
        }
    }
}


void cmor_check_prep(var_t *var, long gridsize, double missval, double *array)
{
  if ( var->checkvalid || var->check_min_mean_abs || var->check_max_mean_abs )
    {
      double aval;
      double amean = 0;
      long nvals = 0;
  
      for ( long i = 0; i < gridsize; ++i )
        {
          aval = array[i];
          if ( !DBL_IS_EQUAL(aval, missval) )
            {
              amean += fabs(aval);
              nvals++;
            }
        }

      var->amean += amean;
      var->nvals += nvals;

      long n_lower_min = 0;
      long n_greater_max = 0;

      for ( long i = 0; i < gridsize; ++i )
        {
          aval = array[i];
          if ( !DBL_IS_EQUAL(aval, missval) )
            {
              if ( aval < var->valid_min ) n_lower_min++;
              if ( aval > var->valid_max ) n_greater_max++;
            }
        }

      var->n_lower_min += n_lower_min;
      var->n_greater_max += n_greater_max;
    }
}

static
void apply_cmorlist(list_t *pmlist, int nvars, int vlistID2, var_t *vars)
{
  const char *hentry[] = {"Header"};
  const char *ventry[] = {"variable_entry", "parameter"};
  int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);
  int nhentry = (int) sizeof(hentry)/sizeof(hentry[0]);
  char varname[CDI_MAX_NAME];

  // search for global missing value
  bool lmissval = false;
  double missval;
  list_t *kvlist = pmlist_get_kvlist_ventry(pmlist, nhentry, hentry);
  if ( kvlist )
    {
      for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
        {
          keyValues_t *kv = *(keyValues_t **)kvnode->data;
          const char *key = kv->key;
          const char *value = (kv->nvalues == 1) ? kv->values[0] : NULL;
          if ( !value || (value && !*value) ) continue;

          if ( STR_IS_EQ(key, "missing_value") )
            {
              lmissval = true;
              missval = parameter2double(kv->values[0]);
            }
          else if ( STR_IS_EQ(key, "table_id") ||
                    STR_IS_EQ(key, "modeling_realm") ||
                    STR_IS_EQ(key, "realm") ||
                    STR_IS_EQ(key, "project_id") ||
                    STR_IS_EQ(key, "frequency") )
            {
              cdiDefAttTxt(vlistID2, CDI_GLOBAL, key, (int)strlen(value), value);
            }
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

      list_t *kvlist = pmlist_search_kvlist_ventry(pmlist, "name", varname, nventry, ventry);
      if ( kvlist )
        {
          bool lvalid_min = false, lvalid_max = false;

          for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
            {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *key = kv->key;
              const char *value = (kv->nvalues == 1) ? kv->values[0] : NULL;
              if ( !value || (value && !*value) ) continue;
              
              // printf("key=%s  value=>%s<\n", key, value);

              if      ( STR_IS_EQ(key, "standard_name") ) vlistDefVarStdname(vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "units")         ) cdo_define_var_units(var, vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "name")          ) /*vlistDefVarName(vlistID2, varID, parameter2word(value))*/;
              else if ( STR_IS_EQ(key, "out_name")      )
                {
                  const char *outname = parameter2word(value);
                  if ( !STR_IS_EQ(var->name, outname) )
                    {
                      vlistDefVarName(vlistID2, varID, outname);
                      cdiDefAttTxt(vlistID2, varID, "original_name", (int)strlen(var->name), var->name);
                    }
                }
              else if ( STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( STR_IS_EQ(key, "comment")       ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( STR_IS_EQ(key, "cell_methods")  ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( STR_IS_EQ(key, "cell_measures") ) cdiDefAttTxt(vlistID2, varID, key, (int)strlen(value), value);
              else if ( STR_IS_EQ(key, "delete")        ) var->remove = parameter2bool(value);
              else if ( STR_IS_EQ(key, "convert")       ) var->convert = parameter2bool(value);
              else if ( STR_IS_EQ(key, "factor")        )
                {
                  var->lfactor = true;
                  var->factor = parameter2double(value);
                  if ( cdoVerbose ) cdoPrint("%s - scale factor %g", varname, var->factor);
                }
              else if ( STR_IS_EQ(key, "missval") || STR_IS_EQ(key, "missing_value") )
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
              else if ( STR_IS_EQ(key, "valid_min") )
                {
                  lvalid_min = true;
                  var->valid_min = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "valid_max") )
                {
                  lvalid_max = true;
                  var->valid_max = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "ok_min_mean_abs") )
                {
                  var->check_min_mean_abs = true;
                  var->ok_min_mean_abs = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "ok_max_mean_abs") )
                {
                  var->check_max_mean_abs = true;
                  var->ok_max_mean_abs = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type") )
                {
                  int datatype = str2datatype(parameter2word(value));
                  if ( datatype != -1 ) vlistDefVarDatatype(vlistID2, varID, datatype);
                }
              else
                {
                  if ( cdoVerbose ) cdoPrint("Attribute %s:%s not supported!", varname,  key);
                }
            }

          if ( lvalid_min && lvalid_max ) var->checkvalid = true;
        }
      else
        {
          cdoPrint("Variable %s not found in CMOR table!", varname);
        }
    }
}


void *CMOR_lite(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  bool delvars = false;
  double missval;

  cdoInitialize(argument);

  CDO_CMOR_Mode = 1;
  if ( CDO_CMOR_Mode ) cdiDefGlobal("CMOR_MODE", CDO_CMOR_Mode);

  cdoOperatorAdd("cmorlite",  0, 0, "parameter table name");

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

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  int nvars = vlistNvars(vlistID2);
  var_t *vars = (var_t *) Malloc(nvars*sizeof(var_t));
  memset(vars, 0, nvars*sizeof(var_t));

  if ( convert_data )
    for ( varID = 0; varID < nvars; ++varID ) vars[varID].convert = true;

  const char *filename = operatorArgv()[0];
  FILE *fp = fopen(filename, "r");
  if ( fp == NULL ) cdoAbort("Open failed on: %s\n", filename);
      
  list_t *pmlist = cmortable_to_pmlist(fp, filename);
  fclose(fp);

  apply_cmorlist(pmlist, nvars, vlistID2, vars);
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

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  /* vlistPrint(vlistID2);*/
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
