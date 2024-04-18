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

      Filedes    codetab         Parameter code table
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void cdo_print_zaxis(int zaxisID);

void cdo_print_attributes(FILE *fp, int cdiID, int varID, int nblanks)
{
#define MAXATT 8192
  int natts;
  char attname[CDI_MAX_NAME];
  int atttype, attlen;
  char atttxt[MAXATT];
  int attint[MAXATT];
  double attflt[MAXATT];
  char fltstr[128];

  cdiInqNatts(cdiID, varID, &natts);

  for ( int ia = 0; ia < natts; ++ia )
    {
      cdiInqAtt(cdiID, varID, ia, attname, &atttype, &attlen);

      if ( atttype == CDI_DATATYPE_INT8  || atttype == CDI_DATATYPE_UINT8  ||
           atttype == CDI_DATATYPE_INT16 || atttype == CDI_DATATYPE_UINT16 ||
           atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32 )
        {
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  cdiInqAttInt(cdiID, varID, attname, attlen, attint);
          fprintf(fp, "%*s", nblanks, "");
	  fprintf(fp, "%s = ", attname);
	  for ( int i = 0; i < attlen; ++i )
	    {
	      if ( i > 0 ) fprintf(fp, ", ");
	      fprintf(fp, "%d", attint[i]);
	    }
	  fprintf(fp, "\n");
	}
      else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
	{
	  if ( attlen > MAXATT ) attlen = MAXATT;
	  cdiInqAttFlt(cdiID, varID, attname, MAXATT, attflt);
          fprintf(fp, "%*s", nblanks, "");
	  fprintf(fp, "%s = ", attname);
	  for ( int i = 0; i < attlen; ++i )
	    {
	      if ( i > 0 ) fprintf(fp, ", ");
              if ( atttype == CDI_DATATYPE_FLT32 )
                fprintf(fp, "%sf", double_to_attstr(CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
              else
                fprintf(fp, "%s", double_to_attstr(CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
	    }
	  fprintf(fp, "\n");
	}
      else if ( atttype == CDI_DATATYPE_TXT )
	{
	  cdiInqAttTxt(cdiID, varID, attname, sizeof(atttxt), atttxt);
	  atttxt[attlen] = 0;
          fprintf(fp, "%*s", nblanks, "");
	  fprintf(fp, "%s = \"%s\"\n", attname, atttxt);
	}
    }
}

static
void printHistory(FILE *fp, int streamID)
{
  int fileID = pstreamFileID(streamID);
  size_t historysize = (size_t) streamInqHistorySize(fileID);
  if ( historysize > 0 )
    {
      char *history = (char*) Malloc(historysize+1);
      history[historysize] = 0;
      streamInqHistoryString(fileID, history);
      fprintf(fp, "  history = %s\n", history);
      Free(history);
    }
}

static
void printSource(FILE *fp, int vlistID, int varID)
{
  /* institute info */
  const char *instptr = institutInqLongnamePtr(vlistInqVarInstitut(vlistID, varID));
  if ( instptr ) fprintf(fp, "  institution = %s\n", instptr);

  /* source info */
  const char *modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
  if ( modelptr ) fprintf(fp, "  source = %s\n", modelptr);
}

static
void partab(FILE *fp, int streamID, int option)
{
  int vlistID = pstreamInqVlist(streamID);
  int varID, datatype = -1;
  char pstr[32];
  char paramstr[32];
  char varname[CDI_MAX_NAME], varlongname[CDI_MAX_NAME], varstdname[CDI_MAX_NAME], varunits[CDI_MAX_NAME];
      
  int nvars = vlistNvars(vlistID);
  bool linebreak = ( option == 4 ) ? false : true;

  if ( option == 2 )
    {
      int natts;
      cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
      if ( natts > 0 )
	{
	  fprintf(fp, "&parameter\n");
	  fprintf(fp, "  name = _GLOBAL_\n");
          printHistory(fp, streamID);
          printSource(fp, vlistID, 0);
	  cdo_print_attributes(fp, vlistID, CDI_GLOBAL, 2);
	  fprintf(fp, "/\n");
	}
    }

  if ( nvars > 1 )
    {
      datatype = vlistInqVarDatatype(vlistID, 0);
      for ( varID = 1; varID < nvars; varID++ )
	{
	  if ( datatype != vlistInqVarDatatype(vlistID, varID) )
	    {
	      datatype = -1;
	      break;
	    }
	}

      if ( datatype != -1 )
	{
	  fprintf(fp, "&parameter");
	  if ( linebreak ) fprintf(fp, "\n");
	  fprintf(fp, "  name = _default_");
	  if ( linebreak ) fprintf(fp, "\n");
	  if ( datatype2str(datatype, pstr) == 0 )
	    {
	      fprintf(fp, "  datatype = %s", pstr);
	      if ( linebreak ) fprintf(fp, "\n");
	    }
	  fprintf(fp, "/\n");
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      fprintf(fp, "&parameter");
      if ( linebreak ) fprintf(fp, "\n");
      
      varname[0]     = 0;
      varlongname[0] = 0;
      varunits[0]    = 0;
      int param = vlistInqVarParam(vlistID, varID);
      double missval = vlistInqVarMissval(vlistID, varID);
      vlistInqVarName(vlistID, varID, varname);
      vlistInqVarStdname(vlistID, varID, varstdname);
      vlistInqVarLongname(vlistID, varID, varlongname);
      vlistInqVarUnits(vlistID, varID, varunits);
            
      fprintf(fp, "  name = %s", varname);
      if ( linebreak ) fprintf(fp, "\n");
      // if ( code   > 0 ) fprintf(fp, "  code=%d\n", code);
      // if ( tabnum > 0 ) fprintf(fp, "  table=%d\n", tabnum);
      if ( param >= 0 )
	{
	  cdiParamToString(param, paramstr, sizeof(paramstr));
	  fprintf(fp, "  param = %s", paramstr);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varstdname) )
	{
	  fprintf(fp, "  standard_name = %s", varstdname);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varlongname) )
	{
	  fprintf(fp, "  long_name = \"%s\"", varlongname);
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( strlen(varunits) )
	{
	  fprintf(fp, "  units = \"%s\"", varunits);
	  if ( linebreak ) fprintf(fp, "\n");
	}

      if ( datatype == -1 )
	if ( datatype2str(vlistInqVarDatatype(vlistID, varID), pstr) == 0 )
	  {
	    fprintf(fp, "  datatype = %s", pstr);
	    if ( linebreak ) fprintf(fp, "\n");
	  }

      int chunktype = vlistInqVarChunkType(vlistID, varID);
      if ( chunktype == CDI_CHUNK_AUTO )
	{
	  fprintf(fp, "  chunktype = auto");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      else if ( chunktype == CDI_CHUNK_GRID )
	{
	  fprintf(fp, "  chunktype = grid");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      if ( chunktype == CDI_CHUNK_LINES )
	{
	  fprintf(fp, "  chunktype = lines");
	  if ( linebreak ) fprintf(fp, "\n");
	}
      
      if ( option == 2 ) cdo_print_attributes(fp, vlistID, varID, 2);
      if ( option == 2 ) 
	fprintf(fp, "  missing_value = %g\n", missval);
      
      if ( !linebreak ) fprintf(fp, "  ");
      fprintf(fp, "/\n");
    }   
}

static
void filedes(int streamID)
{
  printf("\n");
  int filetype = pstreamInqFiletype(streamID);
  switch ( filetype )
    {
    case CDI_FILETYPE_GRB:  printf("  GRIB data\n"); break;
    case CDI_FILETYPE_GRB2: printf("  GRIB2 data\n"); break;
    case CDI_FILETYPE_NC:   printf("  NetCDF data\n"); break;
    case CDI_FILETYPE_NC2:  printf("  NetCDF2 data\n"); break;
    case CDI_FILETYPE_NC4:  printf("  NetCDF4 data\n"); break;
    case CDI_FILETYPE_NC4C: printf("  NetCDF4 classic data\n"); break;
    case CDI_FILETYPE_NC5:  printf("  NetCDF5 data\n"); break;
    case CDI_FILETYPE_SRV:  printf("  SERVICE data\n"); break;
    case CDI_FILETYPE_EXT:  printf("  EXTRA data\n"); break;
    case CDI_FILETYPE_IEG:  printf("  IEG data\n"); break;
    default: printf("  unsupported filetype %d\n" , filetype);
    }

  switch ( filetype )
    {
    case CDI_FILETYPE_SRV:
    case CDI_FILETYPE_EXT:
    case CDI_FILETYPE_IEG:
      {
        int byteorder = pstreamInqByteorder(streamID);
        switch ( byteorder )
          {
          case CDI_BIGENDIAN:    printf("  byteorder is BIGENDIAN\n"); break;
          case CDI_LITTLEENDIAN: printf("  byteorder is LITTLEENDIAN\n"); break;
          default:  printf("  byteorder %d undefined\n", byteorder); break;
          }
       }
    }  
  
  printf("\n");
}


void *Filedes(void *argument)
{
  cdoInitialize(argument);

  // clang-format off
  int GRIDDES  = cdoOperatorAdd("griddes",   0, 0, NULL);
  int GRIDDES2 = cdoOperatorAdd("griddes2",  0, 0, NULL);
  int ZAXISDES = cdoOperatorAdd("zaxisdes",  0, 0, NULL);
  int VCT      = cdoOperatorAdd("vct",       0, 0, NULL);
  int VCT2     = cdoOperatorAdd("vct2",      0, 0, NULL);
  int CODETAB  = cdoOperatorAdd("codetab",   0, 0, NULL);
  int FILEDES  = cdoOperatorAdd("filedes",   0, 0, NULL);
  int VLIST    = cdoOperatorAdd("vlist",     0, 0, NULL);
  int SPARTAB  = cdoOperatorAdd("spartab",   0, 0, NULL);
  int PARTAB   = cdoOperatorAdd("partab",    0, 0, NULL);
  int PARTAB2  = cdoOperatorAdd("partab2",   0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID = pstreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);

  int nvars  = vlistNvars(vlistID);
  int ngrids = vlistNgrids(vlistID);
  int nzaxis = vlistNzaxis(vlistID);

  if ( operatorID == GRIDDES || operatorID == GRIDDES2 )
    {
      int opt = (operatorID == GRIDDES) ? 1 : 0;
      for ( int index = 0; index < ngrids; index++ )
        {
          printf("#\n" "# gridID %d\n" "#\n", index+1);
          cdo_print_grid(vlistGrid(vlistID, index), opt);
        }
    }
  else if ( operatorID == ZAXISDES )
    {
      for ( int index = 0; index < nzaxis; index++ )
        {
          printf("#\n" "# zaxisID %d\n" "#\n", index+1);
          cdo_print_zaxis(vlistZaxis(vlistID, index));
        }
    }
  else if ( operatorID == VCT || operatorID == VCT2 )
    {
      for ( int index = 0; index < nzaxis; index++ )
	{
	  int zaxisID = vlistZaxis(vlistID, index);
	  int type = zaxisInqType(zaxisID);
	  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
	    {
	      int vctsize = zaxisInqVctSize(zaxisID);
	      const double *vct = zaxisInqVctPtr(zaxisID);
		
	      if ( vctsize%2 == 0 )
		{
		  if ( operatorID == VCT )
		    {
		      fprintf(stdout, "#   k         vct_a(k) [Pa]             vct_b(k) []\n");
		      for ( int i = 0; i < vctsize/2; i++ )
			fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
		    }
		  else
		    {
		      fprintf(stdout, "vctsize   = %d\n", vctsize);
		      int nbyte0 = fprintf(stdout, "vct       = ");
		      int nbyte = nbyte0;
		      for ( int i = 0; i < vctsize; i++ )
			{
			  if ( nbyte > 70 || i == vctsize/2 )
			    {
			      fprintf(stdout, "\n%*s", nbyte0, "");
			      nbyte = nbyte0;
			    }
			  nbyte += fprintf(stdout, "%.9g ", vct[i]);
			}
		      fprintf(stdout, "\n");
		    }
		}
	      else
		for ( int i = 0; i < vctsize; i++ )
		  fprintf(stdout, "%5d %25.17f\n", i, vct[i]);

	      break;
	    }
	}
    }
  else if ( operatorID == VLIST )
    {
      vlistPrint(vlistID);
    }
  else if ( operatorID == CODETAB )
    {
      char varname[CDI_MAX_NAME], varlongname[CDI_MAX_NAME], varunits[CDI_MAX_NAME];

      for ( int varID = 0; varID < nvars; varID++ )
	{
	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  int code = vlistInqVarCode(vlistID, varID);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  vlistInqVarUnits(vlistID, varID, varunits);
	  fprintf(stdout, "%4d  %-12s", code, varname);
	  if ( varlongname[0] )
	    {
	      fprintf(stdout, "  %s", varlongname);
	      if ( varunits[0] )
		fprintf(stdout, " [%s]", varunits);
	    }
	  fprintf(stdout, "\n");
	}   
    }
  else if ( operatorID == PARTAB || operatorID == SPARTAB || operatorID == PARTAB2 )
    {
      int option = 1;
      if ( operatorID == SPARTAB ) option = 4;
      if ( operatorID == PARTAB2 ) option = 2;
      
      partab(stdout, streamID, option);
    }
  else if ( operatorID == FILEDES )
    {
      filedes(streamID);
    }

  pstreamClose(streamID);

  cdoFinish();

  return 0;
}
