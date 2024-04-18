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

void printAtts(int vlistID, int varOrGlobal, int natts, char *argument)
{
  char stdname[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  int found = 0;
  if ( varOrGlobal != CDI_GLOBAL ) 
    {
      vlistInqVarStdname(vlistID, varOrGlobal, stdname);
      vlistInqVarLongname(vlistID, varOrGlobal, longname);
      vlistInqVarUnits(vlistID, varOrGlobal, units);
      double misval = vlistInqVarMissval(vlistID, varOrGlobal);
      if ( argument )
        {
          if ( STR_IS_EQ("standard_name", argument) && stdname[0] )
            {
              fprintf(stdout, "  standard_name = \"%s\"\n", stdname);       
              found++;
            }
          else if ( STR_IS_EQ("long_name", argument) && longname[0] )
            {
              fprintf(stdout, "  long_name = \"%s\"\n", longname);       
              found++;
            }
          else if ( STR_IS_EQ("units", argument) && units[0] )
            {
              fprintf(stdout, "  units = \"%s\"\n", units);       
              found++;
            }
          else if ( STR_IS_EQ("missing_value", argument) )
            {
              fprintf(stdout, "  missing_value = \"%e\"\n", misval);       
              found++;
            }
        }
      else
        {
          if ( stdname[0] )
            fprintf(stdout, "  standard_name = \"%s\"\n", stdname);  
          if ( longname[0] )     
            fprintf(stdout, "  long_name = \"%s\"\n", longname);  
          if ( units[0] )    
            fprintf(stdout, "  units = \"%s\"\n", units);       
          fprintf(stdout, "  missing_value = \"%e\"\n", misval); 
        }      
    }
  if ( found == 0 )
    {
      for ( int i = 0; i < natts; i++ )
        {
          char name[CDI_MAX_NAME];
          char *value = (char *)Malloc(4*CDI_MAX_NAME * sizeof(char));
          int type, len;
          cdiInqAtt(vlistID, varOrGlobal, i, name, &type, &len);
          if ( argument )
            {
              if ( !STR_IS_EQ(argument, name) )
                continue;
              else
                found++;
            }
          switch ( type )
	    {
	    case CDI_DATATYPE_TXT:
	      cdiInqAttTxt(vlistID, varOrGlobal, name, len, value);
	      value[len] = '\0';
	      fprintf(stdout, "  %s = \"%s\"\n", name, value);
	      break;
	    case CDI_DATATYPE_INT32:
	      cdiInqAttInt(vlistID, varOrGlobal, name, len, (int *)value);
	      fprintf(stdout, "  %s = %i\n", name, *(int *)value);
	      break;
	    case CDI_DATATYPE_FLT64:
	      cdiInqAttFlt(vlistID, varOrGlobal, name, len, (double *)value);
	      fprintf(stdout, "  %s = %e\n", name, *(double *)value);
    	      break;
    	    default:
    	      cdoWarning("Unsupported type %i name %s\n", type, name);
    	    }
          Free(value);
        }
    }
  if ( argument && found == 0 )
    {
      if ( varOrGlobal != CDI_GLOBAL )
        cdoAbort("Could not find variable attribute %s in infile.", argument);
      else
        cdoAbort("Could not find global attribute %s in infile.", argument);
    }
}

void check_varname_and_print(int vlistID, int nvars, char *checkvarname, char *attname)
{
  int varID = 0;
  for ( varID = 0; varID < nvars; varID++ )
    {
      char filevarname[CDI_MAX_NAME];
      vlistInqVarName(vlistID, varID, filevarname);  
      if ( !checkvarname || STR_IS_EQ(checkvarname, filevarname) )
        {
          fprintf(stdout, "%s:\n", filevarname);
	  int nfileattsvar;
          cdiInqNatts(vlistID, varID, &nfileattsvar);
          printAtts(vlistID, varID, nfileattsvar, attname);
          break;
        }
    } 
  if ( nvars == varID && checkvarname )
    cdoAbort("Could not find variable %s in infile.", checkvarname);
}

void *Showattribute(void *argument)
{
  const int delim = '@';
  cdoInitialize(argument);

  int SHOWATTRIBUTE = cdoOperatorAdd("showattribute",   0, 0, NULL);
  int SHOWATTSVAR   = cdoOperatorAdd("showattsvar",   0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID = pstreamOpenRead(cdoStreamName(0));
  int vlistID = pstreamInqVlist(streamID);
  int nvars   = vlistNvars(vlistID);

  int natts = operatorArgc();
  if ( natts == 0 && operatorID != SHOWATTSVAR)
    cdoAbort("Parameter missing!");
  else if ( natts == 0 )
    check_varname_and_print(vlistID, nvars, NULL, NULL);
  else
    {
      char **params = operatorArgv();
      char buffer[CDI_MAX_NAME];
      for ( int i = 0; i < natts; i++)
        {
          strcpy(buffer, params[i]);
          char *varname = NULL, *input = NULL;
          char *result = strrchr(buffer, delim);
          input = buffer;
          if ( result == NULL )
            {
              if ( operatorID == SHOWATTRIBUTE )
                {
                  fprintf(stdout, "Global:\n");
                  int nfileatts;
                  cdiInqNatts(vlistID, CDI_GLOBAL, &nfileatts);
                  printAtts(vlistID, CDI_GLOBAL, nfileatts, input);
                }
              else if (operatorID == SHOWATTSVAR )
                check_varname_and_print(vlistID, nvars, input, NULL);
            }
          else
            {
              if ( operatorID == SHOWATTRIBUTE )
                {
                  input = result+1;
                  *result = 0;
                  varname = buffer;
                  check_varname_and_print(vlistID, nvars, varname, input);
                }
              else if ( operatorID == SHOWATTSVAR )
                check_varname_and_print(vlistID, nvars, input, NULL);
            }
        }
    }
  pstreamClose(streamID);

  cdoFinish();

  return 0;
}
