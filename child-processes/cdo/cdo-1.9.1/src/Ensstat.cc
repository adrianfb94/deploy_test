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

      Ensstat    ensrange        Ensemble range
      Ensstat    ensmin          Ensemble minimum
      Ensstat    ensmax          Ensemble maximum
      Ensstat    enssum          Ensemble sum
      Ensstat    ensmean         Ensemble mean
      Ensstat    ensavg          Ensemble average
      Ensstat    ensstd          Ensemble standard deviation
      Ensstat    ensstd1         Ensemble standard deviation
      Ensstat    ensvar          Ensemble variance
      Ensstat    ensvar1         Ensemble variance
      Ensstat    enspctl         Ensemble percentiles
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "cdo_task.h"
#include "pstream.h"
#include "util.h"


typedef struct
{
  int streamID;
  int vlistID;
  int nmiss[2];
  double missval[2];
  double *array[2];
} ens_file_t;


typedef struct
{
  int t;
  int varID[2];
  int levelID[2];
  int vlistID1;
  int streamID2;
  int nfiles;
  ens_file_t *ef;
  double *array2;
  double *count2;
  field_type *field;
  int operfunc;
  double pn;
  bool lpctl;
  bool count_data;
  int nvars;
} ensstat_arg_t;


static
void *ensstat_func(void *ensarg)
{
  if ( CDO_task ) cdo_omp_set_num_threads(ompNumThreads);

  ensstat_arg_t *arg = (ensstat_arg_t*) ensarg;
  int t = arg->t;
  int nfiles = arg->nfiles;
  ens_file_t *ef = arg->ef;
  field_type *field = arg->field;

  bool lmiss = false;
  for ( int fileID = 0; fileID < nfiles; fileID++ ) if ( ef[fileID].nmiss[t] > 0 ) lmiss = true;

  int gridID = vlistInqVarGrid(arg->vlistID1, arg->varID[t]);
  int gridsize = gridInqSize(gridID);
  double missval = vlistInqVarMissval(arg->vlistID1, arg->varID[t]);

  int nmiss = 0;
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
  for ( int i = 0; i < gridsize; ++i )
    {
      int ompthID = cdo_omp_get_thread_num();

      field[ompthID].missval = missval;
      field[ompthID].nmiss = 0;
      for ( int fileID = 0; fileID < nfiles; fileID++ )
        {
          field[ompthID].ptr[fileID] = ef[fileID].array[t][i];
          if ( lmiss && DBL_IS_EQUAL(field[ompthID].ptr[fileID], ef[fileID].missval[t]) )
            {
              field[ompthID].ptr[fileID] = missval;
              field[ompthID].nmiss++;
            }
        }

      arg->array2[i] = arg->lpctl ? fldpctl(field[ompthID], arg->pn) : fldfun(field[ompthID], arg->operfunc);

      if ( DBL_IS_EQUAL(arg->array2[i], field[ompthID].missval) )
        {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
          nmiss++;
        }

      if ( arg->count_data ) arg->count2[i] = nfiles - field[ompthID].nmiss;
    }

  pstreamDefRecord(arg->streamID2, arg->varID[t], arg->levelID[t]);
  pstreamWriteRecord(arg->streamID2, arg->array2, nmiss);

  if ( arg->count_data )
    {
      pstreamDefRecord(arg->streamID2, arg->varID[t]+arg->nvars, arg->levelID[t]);
      pstreamWriteRecord(arg->streamID2, arg->count2, 0);
    }

  return NULL;
}


void *Ensstat(void *argument)
{
  void *task = CDO_task ? cdo_task_new() : NULL;
  ensstat_arg_t ensstat_arg;
  int nrecs0;

  cdoInitialize(argument);

  // clang-format off
  cdoOperatorAdd("ensrange", func_range, 0, NULL);
  cdoOperatorAdd("ensmin",   func_min,   0, NULL);
  cdoOperatorAdd("ensmax",   func_max,   0, NULL);
  cdoOperatorAdd("enssum",   func_sum,   0, NULL);
  cdoOperatorAdd("ensmean",  func_mean,  0, NULL);
  cdoOperatorAdd("ensavg",   func_avg,   0, NULL);
  cdoOperatorAdd("ensstd",   func_std,   0, NULL);
  cdoOperatorAdd("ensstd1",  func_std1,  0, NULL);
  cdoOperatorAdd("ensvar",   func_var,   0, NULL);
  cdoOperatorAdd("ensvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("enspctl",  func_pctl,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lpctl = operfunc == func_pctl;

  int argc = operatorArgc();
  int nargc = argc;

  double pn = 0;
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2int(operatorArgv()[0]);
      percentile_check_number(pn);
      argc--;
    }

  bool count_data = false;
  if ( argc == 1 )
    {
      if ( strcmp("count", operatorArgv()[nargc-1]) == 0 ) count_data = true;
      else cdoAbort("Unknown parameter: >%s<", operatorArgv()[nargc-1]); 
    }
    
  int nfiles = cdoStreamCnt() - 1;

  if ( cdoVerbose ) cdoPrint("Ensemble over %d files.", nfiles);

  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename) )
    cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t *) Malloc(nfiles*sizeof(ens_file_t));

  field_type *field = (field_type *) Malloc(ompNumThreads*sizeof(field_type));
  for ( int i = 0; i < ompNumThreads; i++ )
    {
      field_init(&field[i]);
      field[i].size = nfiles;
      field[i].ptr  = (double*) Malloc(nfiles*sizeof(double));
    }

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].streamID = pstreamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID  = pstreamInqVlist(ef[fileID].streamID);
    }

  /* check that the contents is always the same */
  for ( int fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  int vlistID1 = ef[0].vlistID;
  int vlistID2 = vlistDuplicate(vlistID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsizemax = vlistGridsizeMax(vlistID1);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].array[0] = (double*) Malloc(gridsizemax*sizeof(double));
      ef[fileID].array[1] = task ? (double*) Malloc(gridsizemax*sizeof(double)) : NULL;
    }

  double *array2 = (double *) Malloc(gridsizemax*sizeof(double));

  int nvars = vlistNvars(vlistID2);
  double *count2 = NULL;
  if ( count_data )
    {
      count2 = (double *) Malloc(gridsizemax*sizeof(double));
      for ( int varID = 0; varID < nvars; ++varID )
	{
	  char name[CDI_MAX_NAME];
	  vlistInqVarName(vlistID2, varID, name);
	  strcat(name, "_count");
	  int gridID = vlistInqVarGrid(vlistID2, varID);
	  int zaxisID = vlistInqVarZaxis(vlistID2, varID);
	  int timetype = vlistInqVarTimetype(vlistID2, varID);
	  int cvarID = vlistDefVar(vlistID2, gridID, zaxisID, timetype);
	  vlistDefVarName(vlistID2, cvarID, name);
	  vlistDefVarDatatype(vlistID2, cvarID, CDI_DATATYPE_INT16);
	  if ( cvarID != (varID+nvars) ) cdoAbort("Internal error, varIDs do not match!");
	}
    }

  int streamID2 = pstreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  ensstat_arg.vlistID1 = vlistID1;
  ensstat_arg.streamID2 = streamID2;
  ensstat_arg.nfiles = nfiles;
  ensstat_arg.array2 = array2;
  ensstat_arg.count2 = count2;
  ensstat_arg.field = field;
  ensstat_arg.operfunc = operfunc;
  ensstat_arg.pn = pn;
  ensstat_arg.lpctl = lpctl;
  ensstat_arg.count_data = count_data;
  ensstat_arg.nvars = nvars;
  ensstat_arg.t = 0;

  bool lwarning = false;
  bool lerror = false;
  int t = 0;
  int tsID = 0;
  do
    {
      nrecs0 = pstreamInqTimestep(ef[0].streamID, tsID);
      for ( int fileID = 1; fileID < nfiles; fileID++ )
	{
	  int streamID = ef[fileID].streamID;
	  int nrecs = pstreamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    {
	      if ( nrecs == 0 )
                {
                  lwarning = true;
                  cdoWarning("Inconsistent ensemble file, too few time steps in %s!", cdoStreamName(fileID)->args);
                }
	      else if ( nrecs0 == 0 )
                {
                  lwarning = true;
                  cdoWarning("Inconsistent ensemble file, too few time steps in %s!", cdoStreamName(0)->args);
                }
	      else
                {
                  lerror = true;
                  cdoWarning("Inconsistent ensemble file, number of records at time step %d of %s and %s differ!",
                             tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
                }
              goto CLEANUP;
	    }
	}

      if ( nrecs0 > 0 )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);
	  pstreamDefTimestep(streamID2, tsID);
	}

      for ( int recID = 0; recID < nrecs0; recID++ )
	{
          int varID = 0, levelID;

	  for ( int fileID = 0; fileID < nfiles; fileID++ )
	    {
	      pstreamInqRecord(ef[fileID].streamID, &varID, &levelID);
              ef[fileID].missval[t] = vlistInqVarMissval(ef[fileID].vlistID, varID);
	    }
          //#pragma omp parallel for default(none) shared(ef, nfiles)
	  for ( int fileID = 0; fileID < nfiles; fileID++ )
	    {
	      pstreamReadRecord(ef[fileID].streamID, ef[fileID].array[t], &ef[fileID].nmiss[t]);
	    }

          ensstat_arg.ef = ef;
          ensstat_arg.varID[t] = varID;
          ensstat_arg.levelID[t] = levelID;
          if ( task )
            {
              cdo_task_start(task, ensstat_func, &ensstat_arg);
              cdo_task_wait(task);
              //  t = !t;
            }
          else
            {
              ensstat_func(&ensstat_arg);
            }
        }

      tsID++;
    }
  while ( nrecs0 > 0 );

 CLEANUP:

  if ( lwarning ) cdoWarning("Inconsistent ensemble, processed only the first %d timesteps!", tsID);
  if ( lerror ) cdoAbort("Inconsistent ensemble, processed only the first %d timesteps!", tsID);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    pstreamClose(ef[fileID].streamID);

  pstreamClose(streamID2);

  for ( int fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( ef[fileID].array[0] ) Free(ef[fileID].array[0]);
      if ( ef[fileID].array[1] ) Free(ef[fileID].array[1]);
    }

  if ( ef ) Free(ef);
  if ( array2 ) Free(array2);
  if ( count2 ) Free(count2);

  for ( int i = 0; i < ompNumThreads; i++ ) if ( field[i].ptr ) Free(field[i].ptr);
  if ( field ) Free(field);

  if ( task ) cdo_task_delete(task);

  cdoFinish();

  return 0;
}
