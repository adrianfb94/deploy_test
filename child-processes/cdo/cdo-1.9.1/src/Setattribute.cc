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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
void set_attributes(list_t *kvlist, int vlistID)
{
  enum {Undefined=-99};
  const int delim = '@';
  int nvars = vlistNvars(vlistID);
  int ngrids = vlistNgrids(vlistID);
  int nzaxis = vlistNzaxis(vlistID);
  int maxvars = nvars+ngrids*2+nzaxis;
  int *varIDs = (int*) Malloc(maxvars*sizeof(int));

  int kvn = list_size(kvlist);
  char **wname = (char**) Malloc(kvn*sizeof(char*));
  for ( int i = 0; i < kvn; ++i ) wname[i] = NULL;

  char name[CDI_MAX_NAME];
  char buffer[CDI_MAX_NAME];
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      char *varname = NULL, *attname = NULL;
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      strcpy(buffer, kv->key);
      char *result = strrchr(buffer, delim);
      if ( result == NULL )
        {
          attname = buffer;
        }
      else
        {
          attname = result+1;
          *result = 0;
          varname = buffer;
        }

      if ( *attname == 0 ) cdoAbort("Attribute name missing in >%s<!", kv->key);

      int nv = 0;
      int cdiID = Undefined;
      if ( varname && *varname )
        {
          for ( int idx = 0; idx < nvars; idx++ )
            {
              vlistInqVarName(vlistID, idx, name);
              if ( wildcardmatch(varname, name) == 0 )
                {
                  cdiID = vlistID;
                  varIDs[nv++] = idx;
                }
            }

          if ( cdiID == Undefined )
            {
              /*
              for ( int idx = 0; idx < ngrids; idx++ )
                {
                  int gridID = vlistGrid(vlistID, idx);
                  gridInqXname(gridID, name);
                  if ( wildcardmatch(varname, name) == 0 )
                    {
                      cdiID = gridID;
                      varIDs[nv++] = CDI_GLOBAL;
                    }
                  gridInqYname(gridID, name);
                  if ( wildcardmatch(varname, name) == 0 )
                    {
                      cdiID = gridID;
                      varIDs[nv++] = CDI_GLOBAL;
                    }
                }
              */
              for ( int idx = 0; idx < nzaxis; idx++ )
                {
                  int zaxisID = vlistZaxis(vlistID, idx);
                  zaxisInqName(zaxisID, name);
                  if ( wildcardmatch(varname, name) == 0 )
                    {
                      cdiID = zaxisID;
                      varIDs[nv++] = CDI_GLOBAL;
                    }
                }
            }

          if ( cdiID == Undefined )
            {
              bool lwarn = true;
              for ( int i = 0; i < kvn; ++i )
                {
                  if ( wname[i] == NULL )
                    {
                      wname[i] = strdup(varname);
                      break;
                    }
                  if ( STR_IS_EQ(wname[i], varname) )
                    {
                      lwarn = false;
                      break;
                    }
                }
              if ( lwarn )
                {
                  cdoWarning("Variable >%s< not found!", varname);
                }
            }
        }
      else
        {
          cdiID = vlistID;
          varIDs[nv++] = CDI_GLOBAL;
        }

      if ( cdiID != Undefined && nv > 0 )
        {
          const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;
          int nvalues = kv->nvalues;
          if ( nvalues == 1 && !*value ) nvalues = 0;
          int dtype = literals_find_datatype(nvalues, kv->values);

          for ( int idx = 0; idx < nv; ++idx )
            {
              int varID = varIDs[idx];
              // if ( cdoVerbose ) printf("varID, cdiID, attname %d %d %s %d\n", varID, cdiID, attname, (int)strlen(attname));
              if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
                {
                  int *ivals = (int*) Malloc(nvalues*sizeof(int));
                  for ( int i = 0; i < nvalues; ++i ) ivals[i] = literal_to_int(kv->values[i]);
                  cdiDefAttInt(cdiID, varID, attname, dtype, nvalues, ivals);
                  Free(ivals);
                }
              else if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
                {
                  double *dvals = (double*) Malloc(nvalues*sizeof(double));
                  for ( int i = 0; i < nvalues; ++i ) dvals[i] = literal_to_double(kv->values[i]);
                  cdiDefAttFlt(cdiID, varID, attname, dtype, nvalues, dvals);
                  Free(dvals);
                }
              else
                {
                  int len = (value && *value) ? (int) strlen(value) : 0;
                  cdiDefAttTxt(cdiID, varID, attname, len, value);
                }
            }
         }
    }

  Free(varIDs);
  for ( int i = 0; i < kvn; ++i ) if ( wname[i] ) free(wname[i]);
}


void *Setattribute(void *argument)
{
  int nrecs;
  int varID, levelID;

  cdoInitialize(argument);

  cdoOperatorAdd("setattribute", 0, 0, "attributes");

  bool lcopy = UNCHANGED_RECORD;

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int natts = operatorArgc();
  if ( natts == 0 ) cdoAbort("Parameter missing!");

  list_t *pmlist = NULL;
  list_t *kvlist = kvlist_new("SETATTRIBUTES");
  if ( kvlist_parse_cmdline(kvlist, natts, operatorArgv()) != 0 ) cdoAbort("Parse error!");
  if ( cdoVerbose ) kvlist_print(kvlist);

  if ( natts == 1 )
    {
      keyValues_t *kv = *(keyValues_t **)kvlist->head->data;
      if ( STR_IS_EQ(kv->key, "FILE") )
        {
          if ( cdoVerbose ) cdoPrint("Reading attributes from: %s", kv->values[0]);
          const char *filename = parameter2word(kv->values[0]);
          FILE *fp = fopen(filename, "r");
          if ( fp == NULL ) cdoAbort("Open failed on: %s\n", filename);

          pmlist = namelist_to_pmlist(fp, filename);
          if ( pmlist == NULL ) cdoAbort("Parse error!");
          list_destroy(kvlist);
          kvlist = *(list_t **)pmlist->head->data;
          if ( kvlist == NULL ) cdoAbort("Parse error!");;
          fclose(fp);
          if ( cdoVerbose ) kvlist_print(kvlist);
        }
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  set_attributes(kvlist, vlistID2);

  if ( pmlist ) list_destroy(pmlist);
  else          list_destroy(kvlist);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  double *array = NULL;
  if ( ! lcopy )
    {
      size_t gridsize = (size_t) vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamDefRecord(streamID2,  varID,  levelID);
	  
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

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}
