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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cdi.h>
#include "util.h"


int defineTable(const char *tablearg)
{
  const char *tablename = tablearg;

  int tableID = fileExists(tablename) ? tableRead(tablename) : CDI_UNDEFID;

  if ( tableID == CDI_UNDEFID )
    {
      char *tablepath = getenv("CD_TABLEPATH");
      if ( tablepath )
	{
	  int len = sizeof(tablepath) + sizeof(tablename) + 3;
	  char *tablefile = (char*) malloc(len*sizeof(char));
	  strcpy(tablefile, tablepath);
	  strcat(tablefile, "/");
	  strcat(tablefile, tablename);
	  if ( fileExists(tablename) ) tableID = tableRead(tablefile);
          free(tablefile);
	}
    }

  if ( tableID == CDI_UNDEFID ) tableID = tableInq(-1, 0, tablename);

  if ( tableID == CDI_UNDEFID ) cdoAbort("table <%s> not found", tablename);

  return tableID;
}
