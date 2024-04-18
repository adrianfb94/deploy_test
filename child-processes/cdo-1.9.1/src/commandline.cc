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
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "dmemory.h"

static int    gargc = 0;
static char **gargv;

static char *CDO_CommandLine = NULL;


void freeCommandLine(void)
{
  if ( CDO_CommandLine ) Free(CDO_CommandLine);
}


void initCommandLine(void)
{
  size_t maxlen = 1;
  for ( int iarg = 0; iarg < gargc; iarg++ ) maxlen += strlen(gargv[iarg]) + 1;

  CDO_CommandLine = (char*) Malloc(maxlen);
  atexit(freeCommandLine);
  
  char *pargv;
  size_t offset = 0;
  for ( int iarg = 0; iarg < gargc; iarg++ )
    {
      if ( iarg == 0 )
        {
          pargv = strrchr(gargv[0], '/');
          if ( pargv == 0 ) pargv = gargv[0];
          else              pargv++;
        }
      else
        pargv = gargv[iarg];
      
      size_t len = strlen(pargv);
      if ( offset+len+1 > maxlen ) break;
      memcpy(CDO_CommandLine+offset, pargv, len);
      offset += len;
      CDO_CommandLine[offset] = ' ';
      offset++;
    }

  CDO_CommandLine[offset-1] = '\0';
}

char *commandLine(void)
{
  static bool init = false;

  if ( !init )
    {
      initCommandLine();
      init = true;
    }

  return CDO_CommandLine;
}

void setCommandLine(int argc, char **argv)
{
  gargc = argc;
  gargv = argv;
}
