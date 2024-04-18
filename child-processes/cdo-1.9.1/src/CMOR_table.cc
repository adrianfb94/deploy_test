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

*/

#include "cdo.h"
#include "cdo_int.h"
#include "pmlist.h"


static
void dump_cmor_table(list_t *pmlist)
{  
  printf("# Number of lists: %d\n", list_size(pmlist));
  int i = 0;
  for ( listNode_t *pmnode = pmlist->head; pmnode; pmnode = pmnode->next )
    {
      list_t *kvlist = *(list_t **)pmnode->data;
      printf("# list ID: %d;   Number of elements: %d\n", i, list_size(kvlist));
      printf("&%s\n", list_name(kvlist));
      for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
        {
          keyValues_t *kv = *(keyValues_t **)kvnode->data;
          if ( kv ) printf("  %s = %s\n", kv->key, kv->values[0]);
        }
      printf("/\n");
      ++i;
    }
}

static
void conv_cmor_table(list_t *pmlist)
{
  const char *hname = "Header";
  const char *vname = "variable";
  //const char *aname = "axis";

  bool hasmissval = false;
  double missval;

  for ( listNode_t *pmnode = pmlist->head; pmnode; pmnode = pmnode->next )
    {
      list_t *kvlist = *(list_t **)pmnode->data;
      const char *listname = list_name(kvlist);

      if ( strncmp(listname, hname, strlen(hname)) == 0  )
	{
          for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
	    {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *ename  = kv->key;
	      const char *evalue = kv->values[0];
              size_t len = strlen(ename);

	      if ( strncmp("missing_value", ename, len) == 0 )
		{
		  missval = atof(evalue);
		  hasmissval = true;
		}
	    }
	}
      else if ( strncmp(listname, vname, strlen(vname)) == 0 )
	{
	  printf("&%s\n", "parameter");
          for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
	    {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *ename  = kv->key;
	      const char *evalue = kv->values[0];
	      int len = strlen(ename);
	      int vlen = strlen(evalue);

	      if ( vlen > 1 && evalue[0] == '"' && evalue[vlen-1] == '"' ) 
		{
		  vlen -= 2;
		  evalue++;
		}

	      char *ovalue = strdup(evalue);
	      for ( int i = 1; i < vlen; ++i )
		{
		  if ( ovalue[i-1] == '"' && ovalue[i] == '"' )
		    {
		      ovalue [i-1] = '\'';
		      for ( int j = i+1; j < vlen; ++j ) ovalue[j-1] = ovalue[j];
		      vlen -= 1;
		    }
		}

              if ( vlen )
                {
                  if ( strncmp("name", ename, len)            == 0 ||
                       strncmp("standard_name", ename, len)   == 0 ||
                       strncmp("out_name", ename, len)        == 0 ||
                       strncmp("type", ename, len)            == 0 ||
                       strncmp("valid_min", ename, len)       == 0 ||
                       strncmp("valid_max", ename, len)       == 0 ||
                       strncmp("ok_min_mean_abs", ename, len) == 0 ||
                       strncmp("ok_max_mean_abs", ename, len) == 0 )
                    printf("  %-15s = %s\n", ename, ovalue);
                  else if ( strncmp("long_name", ename, len)  == 0 ||
                            strncmp("units", ename, len)           == 0 ||
                            strncmp("cell_methods", ename, len)    == 0 ||
                            strncmp("cell_measures", ename, len)   == 0 ||
                            strncmp("comment", ename, len)         == 0 )
                    printf("  %-15s = \"%.*s\"\n", ename, vlen, ovalue);
                }
              
	      Free(ovalue);
	    }
	  if ( hasmissval ) printf("  %-15s = %g\n", "missing_value", missval);
	  printf("/\n");
	}
    }
}


void *CMOR_table(void *argument)
{
  cdoInitialize(argument);

  int DUMP_CMOR_TABLE = cdoOperatorAdd("dump_cmor_table",   0,   0, NULL);
  int CONV_CMOR_TABLE = cdoOperatorAdd("conv_cmor_table",   0,   0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
  const char *filename = operatorArgv()[0];

  if ( cdoVerbose ) cdoPrint("Parse file: %s", filename);
  
  FILE *fp = fopen(filename, "r");
  if ( fp == NULL ) cdoAbort("Open failed on: %s\n", filename);
      
  list_t *pmlist = cmortable_to_pmlist(fp, filename);
  fclose(fp);

  if      ( operatorID == DUMP_CMOR_TABLE ) dump_cmor_table(pmlist);
  else if ( operatorID == CONV_CMOR_TABLE ) conv_cmor_table(pmlist);

  list_destroy(pmlist);

  cdoFinish();

  return 0;
}
