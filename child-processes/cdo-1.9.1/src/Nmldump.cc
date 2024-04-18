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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"


static
void print_values(int nvalues, char **values)
{
  char fltstr[128];
  if ( nvalues && values )
    {
      int dtype = literals_find_datatype(nvalues, values);
      for ( int i = 0; i < nvalues; ++i )
        {
          if ( i ) printf(", ");
          switch (dtype)
            {
            case CDI_DATATYPE_INT8:  printf("%db", literal_to_int(values[i])); break;
            case CDI_DATATYPE_INT16: printf("%ds", literal_to_int(values[i])); break;
            case CDI_DATATYPE_INT32: printf("%d",  literal_to_int(values[i])); break;
            case CDI_DATATYPE_FLT32: printf("%sf", double_to_attstr(CDO_flt_digits, fltstr, sizeof(fltstr), literal_to_double(values[i]))); break;
            case CDI_DATATYPE_FLT64: printf("%s",  double_to_attstr(CDO_dbl_digits, fltstr, sizeof(fltstr), literal_to_double(values[i]))); break;
            default: printf("\"%s\"", values[i]);
            }
        }
    }
}


void kvldump(list_t *pmlist)
{
  if ( pmlist )
    {
      for ( listNode_t *pmnode = pmlist->head; pmnode; pmnode = pmnode->next )
        {
          if ( pmnode->data )
            {
              list_t *kvlist = *(list_t **)pmnode->data;
              if ( kvlist )
                {
                  const char *listname = list_name(kvlist);
                  if ( listname ) printf("&%s\n", list_name(kvlist));
                  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
                    {
                      keyValues_t *kv = *(keyValues_t **)kvnode->data;
                      const char *key = kv->key;
                      if ( listname ) printf("  ");
                      printf("%s = ", key);
                      print_values(kv->nvalues, kv->values);
                      printf("\n");
                    }
                  if ( listname ) printf("/\n");
                }
            }
        }
    }
}


void *Nmldump(void *argument)
{
  cdoInitialize(argument);

  int NMLDUMP = cdoOperatorAdd("nmldump",  0, 0, NULL);
  int KVLDUMP = cdoOperatorAdd("kvldump",  0, 0, NULL);

  int operatorID = cdoOperatorID();

  list_t *pmlist = namelist_to_pmlist(stdin, "STDIN");

  if ( operatorID == NMLDUMP )
    list_for_each(pmlist, pmlist_print_iter);
  else if ( operatorID == KVLDUMP )
    kvldump(pmlist);

  list_destroy(pmlist);

  cdoFinish();

  return 0;
}
