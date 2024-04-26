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

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(HAVE_FNMATCH_H)
#include <fnmatch.h>
#endif


#include <stdio.h>
#include <string.h>
#include <ctype.h>   /* tolower */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "util.h"


#if ! defined(VERSION)
#define  VERSION  "0.0.1"
#endif
 

/* refactor: moved here from *.c */

int CDO_opterr = 0;      // refactor: moved here from cdo_getopt.cc
const char *CDO_optarg = NULL; // refactor: moved here from cdo_getopt.cc
int CDO_optind = 1;      // refactor: moved here from cdo_getopt.cc


/* refactor: moved here from cdo.cc */

char *Progname;
const char *CDO_Version = "Climate Data Operators version " VERSION" (http://mpimet.mpg.de/cdo)";

int ompNumThreads = 1;

int stdin_is_tty  = 0;
int stdout_is_tty = 0;
int stderr_is_tty = 0;

char *cdoGridSearchDir   = NULL;

int cdoDefaultFileType   = CDI_UNDEFID;
int cdoDefaultDataType   = CDI_UNDEFID;
int cdoDefaultByteorder  = CDI_UNDEFID;
int cdoDefaultTableID    = CDI_UNDEFID;
int cdoDefaultInstID     = CDI_UNDEFID;     // moved here from institution.cc, was UNDEFID
int cdoDefaultTimeType   = CDI_UNDEFID;

int cdoLockIO            = FALSE;
int cdoCheckDatarange    = FALSE;

int CDO_flt_digits       = 7;
int CDO_dbl_digits       = 15;

int CDO_Color            = FALSE;
int CDO_Use_FFTW         = TRUE;
int CDO_Version_Info     = TRUE;
int CDO_CMOR_Mode        = FALSE;

int cdoDiag              = FALSE;

int CDO_Memtype          = MEMTYPE_DOUBLE;
int CDO_Parallel_Read    = FALSE;

int CDO_Reduce_Dim       = FALSE;
int CDO_Append_History   = TRUE;
int CDO_Reset_History    = FALSE;

int cdoCompType          = CDI_COMPRESS_NONE;  // compression type
int cdoCompLevel         = 0;              // compression level
int cdoDebug             = 0;
int cdoChunkType         = CDI_UNDEFID;
int cdoLogOff            = FALSE;
int cdoSilentMode        = FALSE;
int cdoOverwriteMode     = FALSE;
int cdoBenchmark         = FALSE;
int cdoTimer             = FALSE;
int cdoVerbose           = FALSE;
int cdoCompress          = FALSE;
int cdoInteractive       = FALSE;
int cdoParIO             = FALSE;
int cdoRegulargrid       = FALSE;
int cdoDebugExt          = 0;     //  Debug level for the KNMI extensions
int cdoNumVarnames       = 0;
char **cdoVarnames       = NULL;

char CDO_File_Suffix[32];

int cdoExpMode           = -1;
const char *cdoExpName         = NULL;

int timer_read, timer_write;

const char *cdoComment(void)
{
  static char comment[256];
  static bool init = false;

  if ( ! init )
    {
      init = true;

      int size = strlen(CDO_Version);
      strncat(comment, CDO_Version, size);
      comment[size] = 0;
    }

  return comment;
}


#if defined(HAVE_FNMATCH_H)
int wildcardmatch(const char *pattern, const char *string)
{
  return fnmatch(pattern, string, 0);
}
#else
// The wildcardmatch function checks if two given strings match. 
// The first string may contain wildcard characters
// * --> Matches with 0 or more instances of any character or set of characters.
// ? --> Matches with any one character.
// source code from http://www.geeksforgeeks.org/wildcard-character-matching/
int wildcardmatch(const char *w, const char *s)
{
    // If we reach at the end of both strings, we are done
    if ( *w == '\0' && *s == '\0' ) return 0;
 
    // Make sure that the characters after '*' are present in second string.
    // This function assumes that the first string will not contain two consecutive '*'
    if ( *w == '*' && *(w+1) != '\0' && *s == '\0' ) return 1;
 
    // If the first string contains '?', or current characters of both strings match
    if ( (*w == '?' && *s != '\0') || *w == *s ) return wildcardmatch(w+1, s+1);
 
    // If there is *, then there are two possibilities
    // a) We consider current character of second string
    // b) We ignore current character of second string.
    if ( *w == '*' ) return wildcardmatch(w+1, s) || wildcardmatch(w, s+1);

    return 1;
}
#endif

int cdo_omp_get_thread_num(void)
{
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}


void cdo_omp_set_num_threads(int nthreads)
{
#if defined(_OPENMP)
  if (  omp_get_max_threads() != nthreads ) omp_set_num_threads(nthreads);
#endif
}


char *getProgname(char *string)
{
#if defined(_WIN32)
  /*  progname = strrchr(string, '\\'); */
  char *progname = " cdo";
#else
  char *progname = strrchr(string, '/');
#endif

  if ( progname == NULL ) progname = string;
  else                    progname++;

  return progname;
}

char *getOperator(const char *argument)
{
  char *operatorArg = NULL;

  if ( argument )
    {
      size_t len = 1 + strlen(argument);

      operatorArg = (char*) Malloc(len);

      memcpy(operatorArg, argument, len);
    }

  return operatorArg;
}


const char *getOperatorName(const char *operatorArg)
{
  char *operatorName = NULL;

  if ( operatorArg )
    {
      if ( operatorArg[0] == '-' )
      {
          operatorArg++;
      }
      char *commapos = (char *)strchr(operatorArg, ',');
      size_t len = (commapos != NULL) ? (size_t)(commapos - operatorArg) : strlen(operatorArg);

      operatorName = (char*) Malloc(len+1);

      memcpy(operatorName, operatorArg, len);
      operatorName[len] = '\0';
    }

  /*  return operatorName; */
  if(is_alias(operatorName))
  {
    operatorName = get_original(operatorName);
  }
    return operatorName;
}

char *getOperatorArg(const char *xoperator)
{
  char *operatorArg = NULL;

  if ( xoperator )
    {
      char *commapos = (char *)strchr(xoperator, ',');

      if ( commapos )
        {
          size_t len = strlen(commapos+1);
          if ( len )
            {
              operatorArg = (char*) Malloc(len+1);
              strcpy(operatorArg, commapos+1);
            }
        }
    }

  return operatorArg;
}

char *getFileArg(char *argument)
{
  char *fileArg = NULL;

  if ( argument )
    {
      char *blankpos = strchr(argument, ' ');

      if ( blankpos )
        {
          char *parg = blankpos + 1;
          size_t len = strlen(parg);
          fileArg = (char*) Malloc(len+1);
          strcpy(fileArg, parg);
        }
    }

  return fileArg;
}

static
void trim_flt(char *ss)
{
  char *cp = ss;
  if ( *cp == '-' ) cp++;
  while ( isdigit((int)*cp ) || *cp == '.' ) cp++;
  if ( *--cp == '.' ) return;

  char *ep = cp+1;
  while ( *cp == '0' ) cp--;
  cp++;
  if ( cp == ep ) return;
  while ( *ep ) *cp++ = *ep++;
  *cp = '\0';

  return;
}


char *double_to_attstr(int digits, char *str, size_t len, double value)
{
  int ret = snprintf(str, len, "%#.*g", digits, value);
  assert(ret != -1 && ret < (int)len);
  trim_flt(str);
  return str;
}


void input_int(char *arg, int intarr[], int maxint, int *nintfound)
{
  int nint = 0;

  intarr[nint++] = atoi(arg);

  while ( (arg = strchr(arg, ',')) && (nint < maxint) )
    intarr[nint++] = atoi(++arg);
    
  *nintfound = nint;
}

std::string string2lower(std::string str)
{
    std::string lower_case_string = str;
    for(char c : str)
    {
       c = tolower(c); 
    }
    return lower_case_string;
}
void strtolower(char *str)
{
  if ( str )
    for ( size_t i = 0; str[i]; ++i )
      str[i] = (char)tolower((int)str[i]);
}

void strtoupper(char *str)
{
  if ( str )
    for ( size_t i = 0; str[i]; ++i )
      str[i] = (char)toupper((int)str[i]);
}


const char *parameter2word(const char *string)
{
  size_t len = strlen(string);

  for ( size_t i = 0; i < len; ++i )
    {
      int c = string[i];
      if ( iscntrl(c) || isblank(c) )
        cdoAbort("Word parameter >%s< contains invalid character at position %d!", string, i+1);
    }

  if ( len == 0 ) cdoAbort("Word parameter >%s< is empty!", string);

  return string;
}


bool parameter2bool(const char *str)
{
  size_t len = strlen(str);

  if ( len == 1 )
    {
      if ( *str == 't' || *str == 'T' || *str == '1' ) return true;
      if ( *str == 'f' || *str == 'F' || *str == '0' ) return false;
    }
  else if ( len == 4 && (STR_IS_EQ(str, "true")  || STR_IS_EQ(str, "TRUE")  || STR_IS_EQ(str, "True"))  ) return true;
  else if ( len == 5 && (STR_IS_EQ(str, "false") || STR_IS_EQ(str, "FALSE") || STR_IS_EQ(str, "False")) ) return false;

  cdoAbort("Boolean parameter >%s< contains invalid characters!", str);

  return false;
}


double parameter2double(const char *string)
{
  char *endptr = NULL;
  double fval = strtod(string, &endptr);
  if ( *endptr != 0 )
    cdoAbort("Float parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return fval;
}


int parameter2int(const char *string)
{
  char *endptr = NULL;
  int ival = (int) strtol(string, &endptr, 10);
  if ( *endptr != 0 )
    cdoAbort("Integer parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return ival;
}


int parameter2intlist(const char *string)
{
  char *endptr = NULL;
  int ival = (int) strtol(string, &endptr, 10);
  if ( *endptr != 0 && *endptr != '/' && (endptr - string) == 0 )
    cdoAbort("Integer parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return ival;
}


const char *seas_name_dec[4] = {"DJF", "MAM", "JJA", "SON"};
const char *seas_name_jan[4] = {"JFM", "AMJ", "JAS", "OND"};

static int season_start = START_DEC;

int get_season_start(void)
{
  static bool lgetenv = true;

  if ( lgetenv )
    {
      lgetenv = false;
  
      char *envstr = getenv("CDO_SEASON_START");
      if ( envstr )
        {
          if      ( strcmp(envstr, "DEC") == 0 ) season_start = START_DEC;
          else if ( strcmp(envstr, "JAN") == 0 ) season_start = START_JAN;
      
          if ( cdoVerbose )
            {
              if      ( season_start == START_DEC )
                cdoPrint("Set SEASON_START to December");
              else if ( season_start == START_JAN )
                cdoPrint("Set SEASON_START to January");
            }
        }
    }

  return season_start;
}


void get_season_name(const char *seas_name[])
{
  if ( get_season_start() == START_DEC )
    for ( int i = 0; i < 4; ++i ) seas_name[i] = seas_name_dec[i];
  else
    for ( int i = 0; i < 4; ++i ) seas_name[i] = seas_name_jan[i];
}


int month_to_season(int month)
{
  int season_start = get_season_start();
  int seas = -1;

  if ( month < 0 || month > 16 ) cdoAbort("Month %d out of range!", month);

  if ( season_start == START_DEC )
    {
      if ( month <= 12 )
        seas = (month % 12) / 3;
      else
        seas = month - 13;
    }
  else
    {
      if ( month <= 12 )
        seas = (month - 1) / 3;
      else
        seas = month - 13;
    }

  if ( seas < 0 || seas > 3 ) cdoAbort("Season %d out of range!", seas+1);

  return seas;
}

//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>

bool fileExists(const char *restrict filename)
{
  bool status = false;
  struct stat buf;

  if ( stat(filename, &buf) == 0 )
    {
      if ( S_ISREG(buf.st_mode) && buf.st_size > 0 ) status = true;
    }

  return status;
}


bool userFileOverwrite(const char *restrict filename)
{
  bool status = false;

  if ( !cdoSilentMode && stdin_is_tty && stderr_is_tty )
    {
      fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", filename);
      char line[1024];
      readline(stdin, line, 1024);
      char *pline = line;
      while ( isspace((int) *pline) ) pline++;
      int len = (int) strlen(pline);
      if ( len == 3 )
        {
          if ( pline[0] == 'y' && pline[1] == 'e' && pline[2] == 's' )
            status = true;
          else if ( pline[0] == 'Y' && pline[1] == 'E' && pline[2] == 'S' )
            status = true;
        }
      else if ( len == 1 )
        {
          if ( pline[0] == 'y' || pline[0] == 'Y' ) status = true;
        }
    }

  return status;
}


bool ps_lhead = false;
int ps_nch   = 0;
int ps_cval  = -1;

void progressInit(void)
{
  ps_lhead = false;
  ps_nch   = 0;
  ps_cval  = -1;
}


void progressStatus(double offset, double refval, double curval)
{
  if ( cdoSilentMode ) return;
  if ( !stdout_is_tty ) return;

  offset = offset < 0 ? 0: offset;
  offset = offset > 1 ? 1: offset;
  refval = refval < 0 ? 0: refval;
  refval = refval > 1 ? 1: refval;
  curval = curval < 0 ? 0: curval;
  curval = curval > 1 ? 1: curval;

  int ival = (offset + refval*curval)*100;

  if ( ps_cval == -1 )
    {
      ps_nch = fprintf(stdout, "%s: %3d%%", processInqPrompt(), 0);
      fflush(stdout);
      ps_lhead = true;
    }

  if ( ival != ps_cval )
    {
      ps_cval = ival;
      fprintf(stdout, "\b\b\b\b%3d%%", ps_cval);
      fflush(stdout);
    }

  if ( ps_cval == 100 && ps_lhead )
    {
      ps_lhead = false;
      while ( ps_nch-- ) fprintf(stdout, "\b \b");
      fflush(stdout);
    }
}


int datatype2str(int datatype, char *datatypestr)
{
  int status = 0;

  if      ( datatype == CDI_DATATYPE_PACK   ) strcpy(datatypestr, "P0");
  else if ( datatype > 0 && datatype <= 32  ) snprintf(datatypestr, 4, "P%d", datatype);
  else if ( datatype == CDI_DATATYPE_CPX32  ) strcpy(datatypestr, "C32");
  else if ( datatype == CDI_DATATYPE_CPX64  ) strcpy(datatypestr, "C64");
  else if ( datatype == CDI_DATATYPE_FLT32  ) strcpy(datatypestr, "F32");
  else if ( datatype == CDI_DATATYPE_FLT64  ) strcpy(datatypestr, "F64");
  else if ( datatype == CDI_DATATYPE_INT8   ) strcpy(datatypestr, "I8");
  else if ( datatype == CDI_DATATYPE_INT16  ) strcpy(datatypestr, "I16");
  else if ( datatype == CDI_DATATYPE_INT32  ) strcpy(datatypestr, "I32");
  else if ( datatype == CDI_DATATYPE_UINT8  ) strcpy(datatypestr, "U8");
  else if ( datatype == CDI_DATATYPE_UINT16 ) strcpy(datatypestr, "U16");
  else if ( datatype == CDI_DATATYPE_UINT32 ) strcpy(datatypestr, "U32");
  else                                  { strcpy(datatypestr, "-1"); status = -1;}

  return status;
}


int str2datatype(const char *datatypestr)
{
  int datatype = -1;
  size_t len = strlen(datatypestr);

  if ( len > 1 )
    {
      int ilen = atoi(datatypestr+1);
      if      ( strncmp(datatypestr, "P0",  len) == 0 ) datatype = CDI_DATATYPE_PACK;
      else if ( strncmp(datatypestr, "P",     1) == 0 &&
                ilen > 0 && ilen <= 32 )                datatype = atoi(datatypestr+1);
      else if ( strncmp(datatypestr, "C32", len) == 0 ) datatype = CDI_DATATYPE_CPX32;
      else if ( strncmp(datatypestr, "C64", len) == 0 ) datatype = CDI_DATATYPE_CPX64;
      else if ( strncmp(datatypestr, "F32", len) == 0 ) datatype = CDI_DATATYPE_FLT32;
      else if ( strncmp(datatypestr, "F64", len) == 0 ) datatype = CDI_DATATYPE_FLT64;
      else if ( strncmp(datatypestr, "I8",  len) == 0 ) datatype = CDI_DATATYPE_INT8;
      else if ( strncmp(datatypestr, "I16", len) == 0 ) datatype = CDI_DATATYPE_INT16;
      else if ( strncmp(datatypestr, "I32", len) == 0 ) datatype = CDI_DATATYPE_INT32;
      else if ( strncmp(datatypestr, "U8",  len) == 0 ) datatype = CDI_DATATYPE_UINT8;
      else if ( strncmp(datatypestr, "U16", len) == 0 ) datatype = CDI_DATATYPE_UINT16;
      else if ( strncmp(datatypestr, "U32", len) == 0 ) datatype = CDI_DATATYPE_UINT32;
      else if ( strncmp(datatypestr, "real",   len) == 0 ) datatype = CDI_DATATYPE_FLT32;
      else if ( strncmp(datatypestr, "double", len) == 0 ) datatype = CDI_DATATYPE_FLT64;
    }

  return datatype;
}


off_t fileSize(const char *restrict filename)
{
  off_t filesize = 0;

  if ( filename[0] == '(' && filename[1] == 'p' )
    {
    }
  else
    {
      struct stat buf;
      if ( stat(filename, &buf) == 0 ) filesize = buf.st_size;
    }
  
  return filesize;
}


/* 
 * Return the filetype extension (const char) for a given filetype (int)
 * TODO: handle lists of extensions i.e. grb and grb2 for GRIB2-format
 */
const char *filetypeext(int filetype)
{
  switch ( filetype )
    {
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2: return ".grb";
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC5:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C: return ".nc";
    case CDI_FILETYPE_SRV:  return ".srv";
    case CDI_FILETYPE_EXT:  return ".ext";
    case CDI_FILETYPE_IEG:  return ".ieg";
    default:                return "";
    }
}


/*
 * Remove file extension:
 * -------------------------------------------------
 * Remove file extension if it is the expected one
 * Do nothing otherwise
 */
void rm_filetypeext(char *file, const char *ext)
{
  // length of filename
  int namelen = (int) strlen(file);
  // length of the original file extension
  int extlen =  (int) strlen(ext);

  // delete original extension if it is the expected one
  if ( strcmp(&file[namelen-extlen], ext) == 0 )
      file[namelen-extlen] = 0;
}


/*
 * Replace or just add file extension:
 * -------------------------------------------------
 * Replace file extension with new one
 * or just add the new file extension 
 * if the original extension is not the expected one
 */
void repl_filetypeext(char file[], const char *oldext, const char *newext)
{
  // delete original extension if it is the expected one
  rm_filetypeext(file, oldext);

  // add new file extension
  strcat(file, newext);
}


void cdoGenFileSuffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname)
{
  if ( strncmp(CDO_File_Suffix, "NULL", 4) != 0 )
    {
      if ( CDO_File_Suffix[0] != 0 )
        {
          strncat(filesuffix, CDO_File_Suffix, maxlen-1);
        }
      else
        {
          bool lready = false;
          bool lcompsz = false;
          
          if ( filetype == cdoDefaultFileType && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )
            {
              size_t len = 0;
              if ( refname != NULL && *refname != 0 && *refname != '-' && *refname != '.' ) len = strlen(refname);

              if ( len > 2 )
                {
                  char *result = (char *)strrchr(refname, '.');
                  if ( result != NULL && result[1] != 0 )
                    {
                      int firstchar = tolower(result[1]);
                      switch (firstchar)
                        {
                        case 'g':
                          if ( cdoDefaultFileType == CDI_FILETYPE_GRB || cdoDefaultFileType == CDI_FILETYPE_GRB2 ) lready = true;
                          break;
                        case 'n':
                          if ( cdoDefaultFileType == CDI_FILETYPE_NC || cdoDefaultFileType == CDI_FILETYPE_NC2 ||
                               cdoDefaultFileType == CDI_FILETYPE_NC4 || cdoDefaultFileType == CDI_FILETYPE_NC4C ||
                               cdoDefaultFileType == CDI_FILETYPE_NC5 ) lready = true;
                          break;
                        case 's':
                          if ( cdoDefaultFileType == CDI_FILETYPE_SRV ) lready = true;
                          break;
                        case 'e':
                          if ( cdoDefaultFileType == CDI_FILETYPE_EXT ) lready = true;
                          break;
                        case 'i':
                          if ( cdoDefaultFileType == CDI_FILETYPE_IEG ) lready = true;
                          break;
                        }
                    }

                  //if ( lready )  strncat(filesuffix, result, maxlen-1);
		  if ( lready && ((len=strlen(result)) < (maxlen-1)) )
		    {
		      while ( len-- )
			{
			  if ( *result == '.' || isalnum(*result) ) 
			    strncat(filesuffix, result, 1);
			  result++;
			}
		    }
                }
            }

          if ( !lready )
            {
              strncat(filesuffix, streamFilesuffix(cdoDefaultFileType), maxlen-1);
              if ( cdoDefaultFileType == CDI_FILETYPE_GRB && vlistIsSzipped(vlistID) ) lcompsz = true;
            }

          if ( cdoDefaultFileType == CDI_FILETYPE_GRB && cdoCompType == CDI_COMPRESS_SZIP ) lcompsz = true;
          if ( lcompsz ) strncat(filesuffix, ".sz", maxlen-1);
        }
    }
}


int cdoFiletype(void)
{
  if ( cdoDefaultFileType == CDI_UNDEFID )
    {
      cdoDefaultFileType = CDI_FILETYPE_GRB;
      if ( ! cdoSilentMode )
        cdoPrint("Set default filetype to GRIB");
    }

  return cdoDefaultFileType;
}


void cdoSetNAN(double missval, size_t gridsize, double *array)
{
  if ( DBL_IS_NAN(missval) )
    {
      double newmissval = -9e33;
      for ( size_t i = 0; i < gridsize; ++i )
        if ( DBL_IS_EQUAL(array[i], missval) )
          array[i] = newmissval;
    }
}


void minmaxval(long nvals, double *array, int *imiss, double *minval, double *maxval)
{
  double xmin =  DBL_MAX;
  double xmax = -DBL_MAX;

  if ( imiss )
    {
      for ( long i = 0; i < nvals; ++i )
	{
	  if ( ! imiss[i] )
	    {
	      if      ( array[i] > xmax ) xmax = array[i];
	      else if ( array[i] < xmin ) xmin = array[i];
	    }
	}
    }
  else
    {
      xmin = array[0];
      xmax = array[0];
      for ( long i = 1; i < nvals; ++i )
	{
	  if      ( array[i] > xmax ) xmax = array[i];
	  else if ( array[i] < xmin ) xmin = array[i];
	}
    }

  *minval = xmin;
  *maxval = xmax;
}


void cdo_check_round(void)
{
  static bool checked = false;
  if ( !checked )
    {
      checked = true;
#define NMAX 3
      double vals2[NMAX];
      double vals[NMAX] = {2783.333, 1.45678921, 1.54678921};
      double rvals[NMAX] = {2783, 1, 2};
      for ( int i = 0; i < NMAX; ++i )
        {
          vals2[i] = round(vals[i]);
          if ( IS_NOT_EQUAL(vals2[i], rvals[i]) )
            cdoAbort("Function round() produces wrong results!");
        }
    }
}


