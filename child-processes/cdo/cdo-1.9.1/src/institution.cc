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
#include "cdo_int.h"


static
int readInstitution(const char *instfile)
{
  int instID = CDI_UNDEFID;
  char line[1024], *pline;
  int lnr = 0;
  int nvar = 0, maxvar = 4;
  char name[1024], longname[1024];
  int center = CDI_UNDEFID, subcenter = CDI_UNDEFID;

  FILE *instfp = fopen(instfile, "r");

  if ( instfp == NULL ) return instID;

  while ( readline(instfp, line, 1024) )
    {
      lnr++;
      if ( line[0] == '#' ) continue;
      if ( nvar == maxvar ) break;
      nvar++;

      pline = line;
      while ( isspace((int) *pline) ) pline++;

      if ( nvar == 1 )
	{
	  if ( isdigit((int) pline[0]) )
	    maxvar = 4;
	  else
	    maxvar = 2;
	}

      if ( nvar == 1 && maxvar == 4 )
	{
	  center = atoi(pline);
	}

      if ( nvar == 2 && maxvar == 4 )
	{
	  if ( ! isdigit((int) pline[0]) )
	    cdoAbort("wrong format in line %d. Missing subcenter!", lnr);

	  subcenter = atoi(pline);
	}

      if ( (nvar == 3 && maxvar == 4) || (nvar == 1 && maxvar == 2) )
	{
	  strcpy(name, pline);
	}

      if ( (nvar == 4 && maxvar == 4) || (nvar == 2 && maxvar == 2) )
	{
	  strcpy(longname, pline);
	}
    }

  fclose(instfp);

  instID = institutInq(center, subcenter, name, longname);
  if ( instID == CDI_UNDEFID )
    instID = institutDef(center, subcenter, name, longname);

  return instID;
}


void defineInstitution(const char *instarg)
{
  const char *instname = instarg;
  int instID = readInstitution(instname);

  if ( instID == CDI_UNDEFID )
    instID = institutInq(0, 0, instname, NULL);

  if ( instID == CDI_UNDEFID )
    cdoAbort("institution <%s> not found", instname);

  cdoDefaultInstID = instID;
}
