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

#if defined (HAVE_EXECINFO_H)
#include <execinfo.h>
#endif

#include <signal.h>
#include <fenv.h>
/*#include <malloc.h>*/ /* mallopt and malloc_stats */
#include <sys/stat.h>
#if defined(HAVE_GETRLIMIT)
#if defined(HAVE_SYS_RESOURCE_H)
#include <sys/time.h>       /* getrlimit */
#include <sys/resource.h>   /* getrlimit */
#endif
#endif
#include <unistd.h>         /* sysconf, gethostname */

#include <thread>

#if defined(SX)
#define RLIM_T  long long
#else
#define RLIM_T  rlim_t
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "cdo_task.h"

#include "cdo_getopt.h"

#if defined(HAVE_LIBPTHREAD)
#include "pstream_int.h"
#include "pthread_debug.h"
#endif

#include "modules.h"
#include "error.h"
#include "grid_proj.h"

#if defined(_OPENMP)
#  include <omp.h>
#endif

#if ! defined(VERSION)
#  define  VERSION  "0.0.1"
#endif

#define MAX_NUM_VARNAMES 256

#include <string>

static int Debug = 0;
static int Version = 0;
static int Help = 0;
static int DebugLevel = 0;
static int numThreads = 0;
static int timer_total;
static int CDO_netcdf_hdr_pad = 0;
static int CDO_Rusage = 0;
static const char *username;

#ifdef __cplusplus
extern "C" {
#endif
void streamGrbDefDataScanningMode(int scanmode);
#if defined (__cplusplus)
}
#endif

void gridsearch_set_method(const char *methodstr);

#define PRINT_RLIMIT(resource) \
      { \
        int status; \
        struct rlimit rlim; \
        status = getrlimit(resource, &rlim); \
        if ( status == 0 ) \
          { \
            if ( sizeof(RLIM_T) > sizeof(long) ) \
              { \
                fprintf(stderr, "CUR %-15s = %llu\n", #resource, (long long) rlim.rlim_cur); \
                fprintf(stderr, "MAX %-15s = %llu\n", #resource, (long long) rlim.rlim_max); \
              } \
            else \
              { \
                fprintf(stderr, "CUR %-15s = %lu\n", #resource, (long) rlim.rlim_cur); \
                fprintf(stderr, "MAX %-15s = %lu\n", #resource, (long) rlim.rlim_max); \
              } \
          } \
      }

#define ITSME  (strcmp(username, "\x6d\x32\x31\x34\x30\x30\x33") == 0)

static
void cdo_stackframe(void)
{
#if defined HAVE_EXECINFO_H && defined HAVE_BACKTRACE
  void *callstack[32];
  int frames = backtrace(callstack, 32);
  char **messages = backtrace_symbols(callstack, frames);

  fprintf(stderr, "[bt] Execution path:\n");
  if ( messages ) {
    for ( int i = 0; i < frames; ++i )
      fprintf(stderr, "[bt] %s\n", messages[i]);
    free(messages);
  }
#endif
}

static
int cdo_feenableexcept(int excepts)
{
#if defined HAVE_FEENABLEEXCEPT
  int feenableexcept(int);
  int old_excepts = feenableexcept(excepts);
  return old_excepts;
#else
  static fenv_t fenv;
  unsigned new_excepts = ((unsigned)excepts) & FE_ALL_EXCEPT;
  int old_excepts = -1;  // previous masks

  if ( fegetenv(&fenv) ) return -1;
#if defined(HAVE_FENV_T___CONTROL) && defined(HAVE_FENV_T___MXCSR)
  old_excepts = (int) (fenv.__control & FE_ALL_EXCEPT);

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);
#endif

  return ( fesetenv(&fenv) ? -1 : (int)old_excepts );
#endif
}

static
void cdo_sig_handler(int signo)
{
  if ( signo == SIGFPE )
    {
      cdo_stackframe();
      cdoAbort("floating-point exception!");
    }
}

static
void cdo_set_digits(const char *optarg)
{
  char *ptr1 = 0;
  if ( optarg != 0 && (int) strlen(optarg) > 0 && optarg[0] != ',' )
    CDO_flt_digits = (int)strtol(optarg, &ptr1, 10);

  if ( CDO_flt_digits < 1 || CDO_flt_digits > 20 )
    cdoAbort("Unreasonable value for float significant digits: %d", CDO_flt_digits);

  if ( ptr1 && *ptr1 == ',' )
    {
      char *ptr2 = 0;
      CDO_dbl_digits = (int)strtol(ptr1+1, &ptr2, 10);
      if  ( ptr2 == ptr1+1 || CDO_dbl_digits < 1 || CDO_dbl_digits > 20 )
        cdoAbort("Unreasonable value for double significant digits: %d", CDO_dbl_digits);
    }
}

static
void cdo_version(void)
{
  const int   filetypes[] = {CDI_FILETYPE_SRV, CDI_FILETYPE_EXT, CDI_FILETYPE_IEG, CDI_FILETYPE_GRB, CDI_FILETYPE_GRB2, CDI_FILETYPE_NC, CDI_FILETYPE_NC2, CDI_FILETYPE_NC4, CDI_FILETYPE_NC4C, CDI_FILETYPE_NC5};
  const char* typenames[] = {        "srv",        "ext",        "ieg",       "grb1",        "grb2",       "nc1",        "nc2",        "nc4",        "nc4c",        "nc5"};

  fprintf(stderr, "%s\n", CDO_Version);
#if defined(USER_NAME) && defined(HOST_NAME) && defined(SYSTEM_TYPE)
  fprintf(stderr, "Compiled: by %s on %s (%s) %s %s\n", USER_NAME, HOST_NAME, SYSTEM_TYPE, __DATE__, __TIME__);
#endif
#if defined(CXX_COMPILER)
  fprintf(stderr, "CXX Compiler: %s\n", CXX_COMPILER);
#endif
#if defined(CXX_VERSION)
  fprintf(stderr, "CXX version : %s\n", CXX_VERSION);
#endif
#if defined(C_COMPILER)
  fprintf(stderr, "C Compiler: %s\n", C_COMPILER);
#endif
#if defined(C_VERSION)
  fprintf(stderr, "C version : %s\n", C_VERSION);
#endif

  printFeatures();
  printLibraries();

  fprintf(stderr, "Filetypes: ");
  set_text_color(stderr, BRIGHT, GREEN);
  for ( size_t i = 0; i < sizeof(filetypes)/sizeof(int); ++i )
    if ( cdiHaveFiletype(filetypes[i]) ) fprintf(stderr, "%s ", typenames[i]);
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  cdiPrintVersion();
  fprintf(stderr, "\n");
}

static
void cdo_usage(void)
{
  const char *name;

  /*  fprintf(stderr, "%s\n", CDO_Version);*/
  /*  fprintf(stderr, "\n");*/
  fprintf(stderr, "usage : cdo  [Options]  Operator1  [-Operator2  [-OperatorN]]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  set_text_color(stderr, RESET, BLUE);
  fprintf(stderr, "    -a             Generate an absolute time axis\n");
  fprintf(stderr, "    -b <nbits>     Set the number of bits for the output precision\n");
  fprintf(stderr, "                   (I8/I16/I32/F32/F64 for nc1/nc2/nc4/nc4c/nc5; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb1/grb2)\n");
  fprintf(stderr, "                   Add L or B to set the byteorder to Little or Big endian\n");
  fprintf(stderr, "    --cmor         CMOR conform NetCDF output\n");
  fprintf(stderr, "    -C, --color    Colorized output messages\n");
  fprintf(stderr, "    --enableexcept <except>\n");
  fprintf(stderr, "                   Set individual floating-point traps (DIVBYZERO, INEXACT, INVALID, OVERFLOW, UNDERFLOW, ALL_EXCEPT)\n");
  fprintf(stderr, "    -f, --format <format>\n");
  fprintf(stderr, "                   Format of the output file. (grb1/grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg)\n");
  fprintf(stderr, "    -g <grid>      Set default grid name or file. Available grids: \n");
  fprintf(stderr, "                   n<N>, t<RES>, tl<RES>, global_<DXY>, r<NX>x<NY>, g<NX>x<NY>, gme<NI>, lon=<LON>/lat=<LAT>\n");
  fprintf(stderr, "    -h, --help     Help information for the operators\n");
  fprintf(stderr, "    --history      Do not append to NetCDF \"history\" global attribute\n");
  fprintf(stderr, "    --netcdf_hdr_pad, --hdr_pad, --header_pad <nbr>\n");
  fprintf(stderr, "                   Pad NetCDF output header with nbr bytes\n");
  /*
  fprintf(stderr, "    -i <inst>      Institution name/file\n");
  fprintf(stderr, "                   Predefined instituts: ");
  for ( int id = 0; id < institutInqNumber; id++ )
    if ( (name = institutInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");
  */
  /* fprintf(stderr, "    -l <level>     Level file\n"); */
  fprintf(stderr, "    -k <chunktype> NetCDF4 chunk type: auto, grid or lines\n");
  fprintf(stderr, "    -L             Lock IO (sequential access)\n");
  fprintf(stderr, "    -M             Switch to indicate that the I/O streams have missing values\n");
  fprintf(stderr, "    -m <missval>   Set the missing value of non NetCDF files (default: %g)\n", cdiInqMissval());
  fprintf(stderr, "    --no_warnings  Inhibit warning messages\n");
  fprintf(stderr, "    -O             Overwrite existing output file, if checked\n");
  fprintf(stderr, "    --operators    List of all operators\n");
#if defined(_OPENMP)
  fprintf(stderr, "    -P <nthreads>  Set number of OpenMP threads\n");
#endif
  fprintf(stderr, "    --percentile <method>\n");
  fprintf(stderr, "                   Percentile method: nrank, nist, numpy, numpy_lower, numpy_higher, numpy_nearest\n");
  fprintf(stderr, "    --precision <float_digits[,double_digits]>\n");
  fprintf(stderr, "                   Precision to use in displaying floating-point data (default: 7,15)\n");
  fprintf(stderr, "    --reduce_dim   Reduce NetCDF dimensions\n");
  if ( ITSME )
    fprintf(stderr, "    --remap_genweights 0/1\n");
  fprintf(stderr, "    -R, --regular  Convert GRIB1 data from reduced to regular grid (cgribex only)\n");
  fprintf(stderr, "    -r             Generate a relative time axis\n");
  fprintf(stderr, "    -S             Create an extra output stream for the module TIMSTAT. This stream\n");
  fprintf(stderr, "                   contains the number of non missing values for each output period.\n");
  fprintf(stderr, "    -s, --silent   Silent mode\n");
  fprintf(stderr, "    --sortname     Alphanumeric sorting of NetCDF parameter names\n");
  fprintf(stderr, "    -t <codetab>   Set GRIB1 default parameter code table name or file (cgribex only)\n");
  fprintf(stderr, "                   Predefined tables: ");
  for ( int id = 0; id < tableInqNumber(); id++ )
    if ( (name = tableInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    --timestat_date <srcdate>\n");
  fprintf(stderr, "                   Target timestamp (time statistics): first, middle, midhigh or last source timestep.\n");
  fprintf(stderr, "    -V, --version  Print the version number\n");
  fprintf(stderr, "    -v, --verbose  Print extra details for some operators\n");
  fprintf(stderr, "    -W             Print extra warning messages\n");
  fprintf(stderr, "    -z szip        SZIP compression of GRIB1 records\n");
  fprintf(stderr, "       jpeg        JPEG compression of GRIB2 records\n");
  fprintf(stderr, "        zip[_1-9]  Deflate compression of NetCDF4 variables\n");
#ifdef HIRLAM_EXTENSIONS
  fprintf(stderr, "    --Dkext <debLev>   Setting debugLevel for extensions\n");
  fprintf(stderr, "    --outputGribDataScanningMode <mode>   Setting grib scanning mode for data in output file <0, 64, 96>; Default is 64\n");
#endif // HIRLAM_EXTENSIONS
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "  Operators:\n");
  fprintf(stderr, "    Use option --operators for a list of all operators.\n");
  /*
  set_text_color(stderr, RESET, GREEN);
  operatorPrintAll();
  reset_text_color(stderr);
  */

  fprintf(stderr, "\n");
  fprintf(stderr, "  CDO version %s, Copyright (C) 2003-2017 Uwe Schulzweida\n", VERSION);
  //  fprintf(stderr, "  Available from <http://mpimet.mpg.de/cdo>\n");
  fprintf(stderr, "  This is free software and comes with ABSOLUTELY NO WARRANTY\n");
  fprintf(stderr, "  Report bugs to <http://mpimet.mpg.de/cdo>\n");
}

static
void cdo_init_is_tty(void)
{
  struct stat statbuf;
  fstat(0, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdin_is_tty = 1;
  fstat(1, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;
  fstat(2, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stderr_is_tty = 1;
}

static
void cdoPrintHelp(std::vector<std::string> help/*, char *xoperator*/)
{
  if (help.empty())
    fprintf(stderr, "No help available for this operator!\n");
  else
    {
      bool lprint;
      for(unsigned long i =  0; i < help.size(); i++)
        {
          lprint = !(help[i][0] == '\0'  && help[i+1][0] == ' ');
          
          if ( lprint )
            {
              if ( COLOR_STDOUT )
                {
                  if ( (help[i].compare( "NAME")        == 0) ||
                       (help[i].compare( "SYNOPSIS")    == 0) ||
                       (help[i].compare( "DESCRIPTION") == 0) ||
                       (help[i].compare( "OPERATORS")   == 0) ||
                       (help[i].compare( "NAMELIST")    == 0) ||
                       (help[i].compare( "PARAMETER")   == 0) ||
                       (help[i].compare( "ENVIRONMENT") == 0) ||
                       (help[i].compare( "NOTE")        == 0) ||
                       (help[i].compare( "EXAMPLES")    == 0) )
                    {
                      set_text_color(stdout, BRIGHT, BLACK);
                      fprintf(stdout, "%s", help[i].c_str());
                      reset_text_color(stdout);
                      fprintf(stdout, "\n");
                    }
                  else
                    fprintf(stdout, "%s\n", help[i].c_str());
                }
              else
                {
                  fprintf(stdout, "%s\n", help[i].c_str());
                }
            }
        }
    }
}
static
void cdoSetDebug(int level)
{
  /*
    level   0: off
    level   1: on
    level   2: cdi
    level   4: memory
    level   8: file
    level  16: format
    level  32: cdo
    level  64: stream
    level 128: pipe
    level 256: pthread
   */
  cdiDebug(level);

  if ( level == 1 || (level &  32) ) cdoDebug = 1;
  if ( level == 1 || (level &  64) ) pstreamDebug(1);
#if defined(HAVE_LIBPTHREAD)
  if ( level == 1 || (level & 128) ) pipeDebug(1);
  if ( level == 1 || (level & 256) ) Pthread_debug(1);
#endif
}

#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

static
void setDefaultDataType(const char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int nbits = -1;
  enum {D_UINT, D_INT, D_FLT, D_CPX};
  int dtype = -1;

  int datatype = tolower(*datatypestr);
  if      ( datatype == 'i' ) { dtype = D_INT;  datatypestr++; }
  else if ( datatype == 'u' ) { dtype = D_UINT; datatypestr++; }
  else if ( datatype == 'f' ) { dtype = D_FLT;  datatypestr++; }
  else if ( datatype == 'c' ) { dtype = D_CPX;  datatypestr++; }
  else if ( datatype == 'p' ) { datatypestr++; }

  if ( isdigit((int) *datatypestr) )
    {
      nbits = atoi(datatypestr);
      datatypestr += 1;
      if ( nbits >= 10 ) datatypestr += 1;

      if ( dtype == -1 )
        {
          if      ( nbits > 0 && nbits < 32 ) cdoDefaultDataType = nbits;
          else if ( nbits == 32 )
            {
              if ( cdoDefaultFileType == CDI_FILETYPE_GRB )
                cdoDefaultDataType = CDI_DATATYPE_PACK32;
              else
                cdoDefaultDataType = CDI_DATATYPE_FLT32;
            }
          else if ( nbits == 64 ) cdoDefaultDataType = CDI_DATATYPE_FLT64;
          else
            {
              fprintf(stderr, "Unsupported number of bits %d!\n", nbits);
              fprintf(stderr, "Use I8/I16/I32/F32/F64 for nc1/nc2/nc4/nc4c/nc5; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb1/grb2.\n");
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if ( dtype == D_INT )
            {
              if      ( nbits ==  8 ) cdoDefaultDataType = CDI_DATATYPE_INT8;
              else if ( nbits == 16 ) cdoDefaultDataType = CDI_DATATYPE_INT16;
              else if ( nbits == 32 ) cdoDefaultDataType = CDI_DATATYPE_INT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype INT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_UINT )
            {
              if      ( nbits ==  8 ) cdoDefaultDataType = CDI_DATATYPE_UINT8;
              else if ( nbits == 16 ) cdoDefaultDataType = CDI_DATATYPE_UINT16;
              else if ( nbits == 32 ) cdoDefaultDataType = CDI_DATATYPE_UINT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype UINT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_FLT )
            {
              if      ( nbits == 32 ) cdoDefaultDataType = CDI_DATATYPE_FLT32;
              else if ( nbits == 64 ) cdoDefaultDataType = CDI_DATATYPE_FLT64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype FLT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_CPX )
            {
              if      ( nbits == 32 ) cdoDefaultDataType = CDI_DATATYPE_CPX32;
              else if ( nbits == 64 ) cdoDefaultDataType = CDI_DATATYPE_CPX64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype CPX!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
        }
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
        {
          if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;
        }
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
        {
          if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: >%s< !\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}
/*
static
void setDefaultDataTypeByte(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int datatype = -1;

  if ( isdigit((int) *datatypestr) )
    {
      datatype = atoi(datatypestr);
      datatypestr++;

      if      ( datatype == 1 ) cdoDefaultDataType = CDI_DATATYPE_PACK8;
      else if ( datatype == 2 ) cdoDefaultDataType = CDI_DATATYPE_PACK16;
      else if ( datatype == 3 ) cdoDefaultDataType = CDI_DATATYPE_PACK24;
      else if ( datatype == 4 ) cdoDefaultDataType = CDI_DATATYPE_FLT32;
      else if ( datatype == 8 ) cdoDefaultDataType = CDI_DATATYPE_FLT64;
      else
        {
          fprintf(stderr, "Unsupported datatype %d!\n", datatype);
          fprintf(stderr, "Use 4/8 for filetype nc/srv/ext/ieg and 1/2/3 for grb1/grb2.\n");
          exit(EXIT_FAILURE);
        }
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
        {
          if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;setDefaultDataTypeByte
        }
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
        {
          if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: %s!\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}
*/
static
void setDefaultFileType(const char *filetypestr, int labort)
{
  if ( filetypestr )
    {
      const char *ftstr = filetypestr;
      size_t len;

      // clang-format off
      if      ( cmpstrlen(filetypestr, "grb2", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB2;}
      else if ( cmpstrlen(filetypestr, "grb1", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "grb",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "nc2",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC2; }
      else if ( cmpstrlen(filetypestr, "nc4c", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC4C;}
      else if ( cmpstrlen(filetypestr, "nc4",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC4; }
      else if ( cmpstrlen(filetypestr, "nc5",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC5; }
      else if ( cmpstrlen(filetypestr, "nc1",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC;  }
      else if ( cmpstrlen(filetypestr, "nc",   len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC2; }
      else if ( cmpstrlen(filetypestr, "srv",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_SRV; }
      else if ( cmpstrlen(filetypestr, "ext",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_EXT; }
      else if ( cmpstrlen(filetypestr, "ieg",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_IEG; }
      else
        {
          if ( labort )
            {
              fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
              fprintf(stderr, "Available filetypes: grb1/grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg\n");
              exit(EXIT_FAILURE);
            }
          else
            {
              return;
            }
        }
      // clang-format on

      if ( cdoDefaultFileType != CDI_UNDEFID && *ftstr != 0 )
        {
          if ( *ftstr == '_' )
            {
              ftstr++;

              setDefaultDataType(ftstr);
            }
          else
            {
              fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
              fprintf(stderr, "Use format[_nbits] with:\n");
              fprintf(stderr, "    format = grb1, grb2, nc1, nc2, nc4, nc4c, nc5, srv, ext or ieg\n");
              fprintf(stderr, "    nbits  = 32/64 for grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg; 1 - 24 for grb1/grb2\n");
              exit(EXIT_FAILURE);
            }
        }
    }
}

#define NTESTS 11
#include <inttypes.h>
static
int getMemAlignment(void)
{
  int ma = -1;
  double *ptr[NTESTS];
  int64_t iptr;
  size_t tsize[NTESTS] = {1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025};
  size_t ma_check[4] = {8, 16, 32, 64};
  int ma_result[4] = {1, 1, 1, 1};

  for ( int i = 0; i < NTESTS; ++i )
    {
      ptr[i] = (double*) malloc(tsize[i]);
      iptr = (int64_t) ptr[i];
      for ( int k = 0; k < 4; ++k ) if ( iptr%ma_check[k] ) ma_result[k] = 0; 
    }
  for ( int i = 0; i < NTESTS; ++i ) free(ptr[i]);

  for ( int i = NTESTS-1; i >= 0; i-- )
    {
      ptr[i] = (double*) malloc(tsize[i]+5);
      iptr = (int64_t) ptr[i];
      for ( int k = 0; k < 4; ++k ) if ( iptr%ma_check[k] ) ma_result[k] = 0; 
    }
  for ( int i = 0; i < NTESTS; ++i ) free(ptr[i]);

  for ( int k = 0; k < 4; ++k ) if ( ma_result[k] ) ma = ma_check[k];

  return ma;
}


static
void defineCompress(const char *arg)
{
  size_t len = strlen(arg);

  if      ( strncmp(arg, "szip", len) == 0 )
    {
      cdoCompType  = CDI_COMPRESS_SZIP;
      cdoCompLevel = 0;
    }
  else if ( strncmp(arg, "jpeg", len) == 0 )
    {
      cdoCompType = CDI_COMPRESS_JPEG;
      cdoCompLevel = 0;
    }
  else if ( strncmp(arg, "gzip", len) == 0 )
    {
      cdoCompType  = CDI_COMPRESS_GZIP;
      cdoCompLevel = 6;
    }
  else if ( strncmp(arg, "zip", 3) == 0 )
    {
      cdoCompType  = CDI_COMPRESS_ZIP;
      if ( len == 5 && arg[3] == '_' && isdigit(arg[4]) )
        cdoCompLevel = atoi(&arg[4]);
      else
        cdoCompLevel = 1;
    }
  else
    {
      fprintf(stderr, "Compression type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static
void defineChunktype(const char *arg)
{
  if      ( strcmp("auto",  arg)   == 0 ) cdoChunkType = CDI_CHUNK_AUTO;
  else if ( strcmp("grid",  arg)   == 0 ) cdoChunkType = CDI_CHUNK_GRID;
  else if ( strcmp("lines", arg)   == 0 ) cdoChunkType = CDI_CHUNK_LINES;
  else
    {
      fprintf(stderr, "Chunk type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static
void defineVarnames(const char *arg)
{
  size_t len = strlen(arg);
  size_t istart = 0;
  while ( istart < len && (arg[istart] == ' ' || arg[istart] == ',') ) istart++;

  len -= istart;

  if ( len )
    {      
      cdoVarnames = (char **) Malloc(MAX_NUM_VARNAMES*sizeof(char *));

      char *pbuf = strdup(arg+istart);
      cdoVarnames[cdoNumVarnames++] = pbuf;    

      char *commapos = pbuf;
      while ( (commapos = strchr(commapos, ',')) != NULL )
        {
          *commapos++ = '\0';
          if ( strlen(commapos) )
            {
              if ( cdoNumVarnames >= MAX_NUM_VARNAMES )
                cdoAbort("Too many variable names (limit=%d)!", MAX_NUM_VARNAMES);

              cdoVarnames[cdoNumVarnames++] = commapos;
            }
        }
      /*
      for ( int i = 0; i < cdoNumVarnames; ++i )
        printf("varname %d: %s\n", i+1, cdoVarnames[i]);
      */
    }
}

static
void get_env_vars(void)
{
  username = getenv("LOGNAME");
  if ( username == NULL )
    {
      username = getenv("USER");
      if ( username == NULL ) username = "unknown";
    }

  char *envstr = getenv("CDO_GRID_SEARCH_DIR");
  if ( envstr )
    {
      size_t len = strlen(envstr);
      if ( len > 0 )
        {
          len += 2;
          cdoGridSearchDir = (char*) Malloc(len);
          memcpy(cdoGridSearchDir, envstr, len-1);
          if ( cdoGridSearchDir[len-3] != '/' )
            {
              cdoGridSearchDir[len-2] = '/';
              cdoGridSearchDir[len-1] = 0;
            }
        }
    }

  envstr = getenv("CDO_LOG_OFF");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          cdoLogOff = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_LOG_OFF         = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DISABLE_HISTORY");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          CDO_Reset_History = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DISABLE_HISTORY = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_RESET_HISTORY");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          CDO_Reset_History = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_RESET_HISTORY = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_HISTORY_INFO");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Append_History = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_HISTORY_INFO = %s\n", envstr);
        }
    }

  CDO_File_Suffix[0] = 0;

  envstr = getenv("CDO_FILE_SUFFIX");
  if ( envstr )
    {
      if ( envstr[0] )
        {
          strncat(CDO_File_Suffix, envstr, sizeof(CDO_File_Suffix)-1);
          if ( cdoVerbose )
            fprintf(stderr, "CDO_FILE_SUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DISABLE_FILESUFFIX");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          strcat(CDO_File_Suffix, "NULL");
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DISABLE_FILESUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DIAG");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          cdoDiag = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DIAG = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_USE_FFTW");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Use_FFTW = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_Use_FFTW = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_VERSION_INFO");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Version_Info = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_Version_Info = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_COLOR");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Color = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_COLOR = %s\n", envstr);
        }
    }
  else if ( CDO_Color == FALSE && ITSME ) CDO_Color = TRUE;
}

static
void print_system_info()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "CDO_Color           = %d\n", CDO_Color);
  fprintf(stderr, "CDO_Reset_History   = %d\n", CDO_Reset_History);
  fprintf(stderr, "CDO_File_Suffix     = %s\n", CDO_File_Suffix);
  fprintf(stderr, "cdoDefaultFileType  = %d\n", cdoDefaultFileType);
  fprintf(stderr, "cdoDefaultDataType  = %d\n", cdoDefaultDataType);
  fprintf(stderr, "cdoDefaultByteorder = %d\n", cdoDefaultByteorder);
  fprintf(stderr, "cdoDefaultTableID   = %d\n", cdoDefaultTableID);
  fprintf(stderr, "\n");

  const char *envstr;
  envstr = getenv("HOSTTYPE");
  if ( envstr ) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if ( envstr ) fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if ( envstr ) fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if ( envstr ) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  fprintf(stderr, "Predefined: __SSE2__\n");
#endif 
  fprintf(stderr, "\n");

  fprintf(stderr, "mem alignment       = %d\n\n", getMemAlignment());

#if defined(HAVE_MMAP)
  fprintf(stderr, "HAVE_MMAP\n");
#endif
#if defined(HAVE_MEMORY_H)
  fprintf(stderr, "HAVE_MEMORY_H\n");
#endif
  fprintf(stderr, "\n");

#if defined(_OPENACC)
  fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
  /* OPENMP 3:  201107 */
  /* OPENMP 4:  201307 gcc 4.9 */
#if defined(_OPENMP)
  fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
#if defined(__cplusplus)
  fprintf(stderr, "__cplusplus         = %ld\n", __cplusplus);
#endif
#if defined(__GNUC__)
  fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#if defined(__GNUC_MINOR__)
  fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#if defined(__ICC)
  fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#if defined(__STDC__)
  fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#if defined(__STD_VERSION__)
  fprintf(stderr, "STD VERSION         = %ld\n", __STD_VERSION__);
#endif
#if defined(__STDC_VERSION__)
  fprintf(stderr, "STDC VERSION        = %ld\n", __STDC_VERSION__);
#endif
#if defined(__STD_HOSTED__)
  fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#if defined(FLT_EVAL_METHOD)
  fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#if defined(FP_FAST_FMA)
  fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#if defined(__FAST_MATH__)
  fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  fprintf(stderr, "\n");

#if defined(_SC_VERSION)
  fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#if defined(_SC_ARG_MAX)
  fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#if defined(_SC_CHILD_MAX)
  fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#if defined(_SC_STREAM_MAX)
  fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#if defined(_SC_OPEN_MAX)
  fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#if defined(_SC_PAGESIZE)
  fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  fprintf(stderr, "\n");

#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_FSIZE)
  PRINT_RLIMIT(RLIMIT_FSIZE);
#endif
#if defined(RLIMIT_NOFILE)
  PRINT_RLIMIT(RLIMIT_NOFILE);
#endif
#if defined(RLIMIT_STACK)
  PRINT_RLIMIT(RLIMIT_STACK);
#endif
#endif
  fprintf(stderr, "\n");
}


static
void check_stacksize()
{
#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_STACK)
  {
    struct rlimit rlim;
    int status = getrlimit(RLIMIT_STACK, &rlim);
    if ( status == 0 )
      {
#define  MIN_STACK_SIZE  67108864L  /* 64MB */
        RLIM_T min_stack_size = MIN_STACK_SIZE;
        if ( min_stack_size > rlim.rlim_max ) min_stack_size = rlim.rlim_max;
        if ( rlim.rlim_cur < min_stack_size )
          {
            rlim.rlim_cur = min_stack_size;

            status = setrlimit(RLIMIT_STACK, &rlim);
            if ( Debug )
              {
                if ( status == 0 )
                  {
                    fprintf(stderr, "Set stack size to %ld\n", (long) min_stack_size);
                    PRINT_RLIMIT(RLIMIT_STACK);
                  }
                else
                  fprintf(stderr, "Set stack size to %ld failed!\n", (long) min_stack_size);
                fprintf(stderr, "\n");
              }
          }
      }
  }
#endif
#endif
}

static
void cdo_set_options(void)
{
  if ( Debug )
    {
      fprintf(stderr, "CDO_CMOR_Mode       = %d\n", CDO_CMOR_Mode);
      fprintf(stderr, "CDO_netcdf_hdr_pad  = %d\n", CDO_netcdf_hdr_pad);
      fprintf(stderr, "\n");
    }
  
  if ( CDO_CMOR_Mode )          cdiDefGlobal("CMOR_MODE", CDO_CMOR_Mode);
  if ( CDO_Reduce_Dim )         cdiDefGlobal("REDUCE_DIM", CDO_Reduce_Dim);
  if ( CDO_netcdf_hdr_pad > 0 ) cdiDefGlobal("NETCDF_HDR_PAD", CDO_netcdf_hdr_pad);  
}


static
long str_to_int(const char *intstring)
{
  long intval = -1;

  if ( intstring )
    {
      long fact = 1;
      int len = (int) strlen(intstring);
      for ( int loop = 0; loop < len; loop++ )
        {
          if ( ! isdigit((int) intstring[loop]) )
            {
              switch ( tolower((int) intstring[loop]) )
                {
                case 'k':  fact = 1024;        break;
                case 'm':  fact = 1048576;     break;
                case 'g':  fact = 1073741824;  break;
                default:   fact = 0;           break;
                }
              break;
            }
        }

      if ( fact ) intval = fact*atol(intstring);
    }

  return intval;
}

static
int parse_options_long(int argc, char *argv[])
{
  int c;
  int lnetcdf_hdr_pad;
  int luse_fftw;
  int lgridsearchnn;
  int lgridsearchradius;
  int lremap_genweights;
  int lprecision;
  int lpercentile;
  int lprintoperatorsno = 0;
  int lprintoperators = 0;
  int lenableexcept;
  int ltimestat_date;
  int ltimestat_bounds;
  int lsortname;
  int lsortparam;
  int ldebLevel;
  int lscmode;

  // clang-format off
  struct cdo_option opt_long[] =
    {
      { "precision",         required_argument,        &lprecision,   1  },
      { "percentile",        required_argument,        &lpercentile,  1  },
      { "netcdf_hdr_pad",    required_argument,    &lnetcdf_hdr_pad,  1  },
      { "header_pad",        required_argument,    &lnetcdf_hdr_pad,  1  },
      { "hdr_pad",           required_argument,    &lnetcdf_hdr_pad,  1  },
      { "use_fftw",          required_argument,          &luse_fftw,  1  },
      { "gridsearchnn",      required_argument,      &lgridsearchnn,  1  },
      { "gridsearchradius",  required_argument,  &lgridsearchradius,  1  },
      { "remap_genweights",  required_argument,  &lremap_genweights,  1  },
      { "enableexcept",      required_argument,      &lenableexcept,  1  },
      { "timestat_date",     required_argument,     &ltimestat_date,  1  },
      { "timestat_bounds",         no_argument,   &ltimestat_bounds,  1  },
      { "cmor",                    no_argument,      &CDO_CMOR_Mode,  1  },
      { "reduce_dim",              no_argument,     &CDO_Reduce_Dim,  1  },
      { "float",                   no_argument,        &CDO_Memtype,  MEMTYPE_FLOAT  },
      { "rusage",                  no_argument,         &CDO_Rusage,  1  },
      { "operators_no_output",     no_argument,  &lprintoperatorsno,  1  },
      { "operators",               no_argument,    &lprintoperators,  1  },
      { "no_warnings",             no_argument,           &_Verbose,  0  },
      { "color",                   no_argument,                NULL, 'C' },
      { "format",            required_argument,                NULL, 'f' },
      { "help",                    no_argument,                NULL, 'h' },
      { "history",                 no_argument, &CDO_Append_History,  0  },
      { "no_history",              no_argument, &CDO_Append_History,  0  },
      { "regular",                 no_argument,                NULL, 'R' },
      { "silent",                  no_argument,                NULL, 's' },
      { "sort",                    no_argument,                NULL, 'Q' },
      { "sortname",                no_argument,          &lsortname,  1  },
      { "sortparam",               no_argument,         &lsortparam,  1  },
      { "table",             required_argument,                NULL, 't' },
      { "verbose",                 no_argument,                NULL, 'v' },
      { "version",                 no_argument,                NULL, 'V' },
      { "Dkext",             required_argument,          &ldebLevel,  1  },
      { "outputGribDataScanningMode", required_argument,  &lscmode,   1  },
      { NULL,                                0,                NULL,  0  }
    };
  // clang-format on

  CDO_opterr = 1;

  while ( 1 )
    {
      // IMPORTANT: BY EVERY OPTION that takes arguments you MUST set its trigger variable to ZERO;
      // otherwise the parameters of other options get wrongly assigned.
      lprecision = 0;
      lpercentile = 0;
      lnetcdf_hdr_pad = 0;
      luse_fftw = 0;
      lgridsearchnn = 0;
      lgridsearchradius = 0;
      lremap_genweights = 0;
      lenableexcept = 0;
      ltimestat_date = 0;
      ltimestat_bounds = 0;
      lsortname = 0;
      lsortparam = 0;
      ldebLevel = 0;
      lscmode = 0;

      c = cdo_getopt_long(argc, argv, "f:b:e:P:g:i:k:l:m:n:t:D:z:aBCcdhLMOpQRrsSTuVvWXZ", opt_long, NULL);
      if ( c == -1 ) break;

      switch (c)
        {
        case '?':
          //cdo_usage();
          //fprintf(stderr, "Illegal option!\n");
          return -1;
          // break;
        case ':':
          //cdo_usage();
          //fprintf(stderr, "Option requires an argument!\n");
          return -1;
          // break;
        case 0:
          if ( lnetcdf_hdr_pad )
            {
              int netcdf_hdr_pad = str_to_int(CDO_optarg);
              if ( netcdf_hdr_pad >= 0 ) CDO_netcdf_hdr_pad = netcdf_hdr_pad;
            }
          else if ( lprecision )
            {
              cdo_set_digits(CDO_optarg);
            }
          else if ( lpercentile )
            {
              percentile_set_method(CDO_optarg);
            }
          else if ( lenableexcept )
            {
              int except = -1;
              if      ( strcmp(CDO_optarg, "DIVBYZERO")  == 0 ) except = FE_DIVBYZERO;
              else if ( strcmp(CDO_optarg, "INEXACT")    == 0 ) except = FE_INEXACT;
              else if ( strcmp(CDO_optarg, "INVALID")    == 0 ) except = FE_INVALID;
              else if ( strcmp(CDO_optarg, "OVERFLOW")   == 0 ) except = FE_OVERFLOW;
              else if ( strcmp(CDO_optarg, "UNDERFLOW")  == 0 ) except = FE_UNDERFLOW;
              else if ( strcmp(CDO_optarg, "ALL_EXCEPT") == 0 ) except = FE_ALL_EXCEPT;
              if ( except < 0 ) cdoAbort("option --%s: unsupported argument: %s", "enableexcept", CDO_optarg);
              cdo_feenableexcept((unsigned)except);
              if ( signal(SIGFPE, cdo_sig_handler) == SIG_ERR ) cdoWarning("can't catch SIGFPE!");
            }
          else if ( ltimestat_date )
            {
              int timestatdate = -1;
              if      ( strcmp(CDO_optarg, "first")   == 0 ) timestatdate = TIMESTAT_FIRST;
              else if ( strcmp(CDO_optarg, "last")    == 0 ) timestatdate = TIMESTAT_LAST;
              else if ( strcmp(CDO_optarg, "middle")  == 0 ) timestatdate = TIMESTAT_MEAN;
              else if ( strcmp(CDO_optarg, "midhigh") == 0 ) timestatdate = TIMESTAT_MIDHIGH;
              if ( timestatdate < 0 ) cdoAbort("option --%s: unsupported argument: %s", "timestat_date", CDO_optarg);
              extern int CDO_Timestat_Date;
              CDO_Timestat_Date = timestatdate;
            }
          else if ( ltimestat_bounds )
            {
              extern bool CDO_Timestat_Bounds;
              CDO_Timestat_Bounds = true;
            }
          else if ( luse_fftw )
            {
              int intarg = parameter2int(CDO_optarg);
              if ( intarg != 0 && intarg != 1 )
                cdoAbort("Unsupported value for option --use_fftw=%d [range: 0-1]", intarg);
              CDO_Use_FFTW = intarg;
            }
          else if ( lgridsearchnn )
            {
              gridsearch_set_method(CDO_optarg);
            }
          else if ( lgridsearchradius )
            {
              extern double gridsearch_radius;
              double fval = radius_str_to_deg(CDO_optarg);
              if ( fval < 0 || fval > 180 ) cdoAbort("%s=%g out of bounds (0-180 deg)!", "gridsearchradius", fval);
              gridsearch_radius = fval;
            }
          else if ( lremap_genweights )
            {
              int intarg = parameter2int(CDO_optarg);
              if ( intarg != 0 && intarg != 1 )
                cdoAbort("Unsupported value for option --remap_genweights %d [0/1]", intarg);
              remap_genweights = intarg;
            }
          else if ( lsortname )
            {
              cdiDefGlobal("SORTNAME", TRUE);
            }
          else if ( lsortparam )
            {
              cdiDefGlobal("SORTPARAM", TRUE);
            }
#ifdef HIRLAM_EXTENSIONS
          else if ( ldebLevel )
            {
              int newDebLevelVal = parameter2int(CDO_optarg);
              if ( newDebLevelVal > 0 )
                {
                  extern int cdiDebugExt;
                  cdoDebugExt = newDebLevelVal;
                  cdiDebugExt = newDebLevelVal;
                }
            }
          else if ( lscmode )
            {
              int scanningModeValue = atoi(CDO_optarg);
              if ( cdoDebugExt ) printf("scanningModeValue=%d\n", scanningModeValue);
              
              if ( (scanningModeValue==0) || (scanningModeValue==64) || (scanningModeValue==96) )
                {
                  streamGrbDefDataScanningMode(scanningModeValue); // -1: not used; allowed modes: <0, 64, 96>; Default is 64
                }
              else
                {
                  cdoAbort("Warning: %d not in allowed modes: <0, 64, 96>; Using default: 64\n", scanningModeValue);
                  streamGrbDefDataScanningMode(64);
                }
            }
#endif
          break;
        case 'a':
          cdoDefaultTimeType = TAXIS_ABSOLUTE;
          break;
        case 'b':
          setDefaultDataType(CDO_optarg);
          break;
        case 'B':
          cdoBenchmark = TRUE;
          break;
        case 'C':
          CDO_Color = TRUE;
          break;
        case 'c':
          cdoCheckDatarange = TRUE;
          break;
        case 'd':
          Debug = 1;
          break;
        case 'D':
          Debug = 1;
          DebugLevel = atoi(CDO_optarg);
          break;
        case 'e':
          {
#if defined(HAVE_GETHOSTNAME)
          char host[1024];
          gethostname(host, sizeof(host));
          cdoExpName = CDO_optarg;
          /* printf("host: %s %s\n", host, cdoExpName); */
          if ( strcmp(host, cdoExpName) == 0 )
            cdoExpMode = CDO_EXP_REMOTE;
          else
            cdoExpMode = CDO_EXP_LOCAL;
#else
          fprintf(stderr, "Function gethostname not available!\n");
          exit(EXIT_FAILURE);
#endif
          break;
          }
        case 'f':
          setDefaultFileType(CDO_optarg, 1);
          break;
        case 'g':
          cdo_set_grids(CDO_optarg);
          break;
        case 'h':        
          Help = 1;
          break;
        case 'i':
          defineInstitution(CDO_optarg);
          break;
        case 'k':
          defineChunktype(CDO_optarg);
          break;
        case 'L':        
          cdoLockIO = TRUE;
          break;
        case 'l':
          defineZaxis(CDO_optarg);
          break;
        case 'm':
          cdiDefMissval(atof(CDO_optarg));
          break;
        case 'M':
          cdiDefGlobal("HAVE_MISSVAL", TRUE);
          break;
        case 'n':
          defineVarnames(CDO_optarg);
          break;
        case 'O':
          cdoOverwriteMode = TRUE;
          break;
        case 'P':
          if ( *CDO_optarg < '1' || *CDO_optarg > '9' )
            {
              fprintf(stderr, "Unexpected character in number of OpenMP threads (-P <nthreads>): %s!\n", CDO_optarg);
              exit(EXIT_FAILURE);
            }
          numThreads = atoi(CDO_optarg);
          break;
        case 'p':
          CDO_Parallel_Read = TRUE;
          CDO_task = true;
          break;
        case 'Q':
          cdiDefGlobal("SORTNAME", TRUE);
          break;
        case 'R':
          cdoRegulargrid = TRUE;
          cdiDefGlobal("REGULARGRID", TRUE);
          break;
        case 'r':
          cdoDefaultTimeType = TAXIS_RELATIVE;
          break;
        case 'S':
          cdoDiag = TRUE;
          break;
        case 's':
          cdoSilentMode = TRUE;
          break;
        case 'T':
          cdoTimer = TRUE;
          break;
        case 't':
          cdoDefaultTableID = defineTable(CDO_optarg);
          break;
        case 'u':
          cdoInteractive = TRUE;
          break;
        case 'V':
          Version = 1;
          break;
        case 'v':
          cdoVerbose = TRUE;
          _Verbose = 1;
          break;
        case 'W': /* Warning messages */
          _Verbose = 1;
          break;
        case 'X': /* multi threaded I/O */
          cdoParIO = TRUE;
          break;
        case 'Z':
          cdoCompress = TRUE;
          break;
        case 'z':
          defineCompress(CDO_optarg);
          break;
        }
    }

  if ( lprintoperators || lprintoperatorsno )
    {
      set_text_color(stderr, RESET, GREEN);
      bool print_no_output = lprintoperatorsno > 0;
      operatorPrintList(print_no_output);
      //operatorPrintAll();
      reset_text_color(stderr);
      return 1;
    }
 
  return 0;
}

static
void cdo_rusage(void)
{
#if defined(HAVE_SYS_RESOURCE_H) && defined(RUSAGE_SELF)
  struct rusage ru;
  int status = getrusage(RUSAGE_SELF, &ru);

  if ( status == 0 )
    {
      double ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
      double st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

      fprintf(stderr, "  User time:     %.3f seconds\n", ut);
      fprintf(stderr, "  System time:   %.3f seconds\n", st);
      fprintf(stderr, "  Total time:    %.3f seconds\n", ut+st);
      fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss/(1024.*1024.));
      fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
      fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
      fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
      fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
      fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
    }
#endif
}


int main(int argc, char *argv[])
{
  int lstop = FALSE;
  int noff = 0;
  int status = 0;
  const char *operatorArg = NULL;
  argument_t *argument = NULL;

  cdo_init_is_tty();

  memExitOnError();

  _Verbose = 1;
  CDO_Reduce_Dim = 0;

  /* mallopt(M_MMAP_MAX, 0); */
 
  setCommandLine(argc, argv);

  Progname = getProgname(argv[0]);

  if ( strncmp(Progname, "cdo", 3) == 0 && strlen(Progname) > 3 ) noff = 3;

  if ( noff ) setDefaultFileType(Progname+noff, 0);

  get_env_vars();
  init_modules();
  status = parse_options_long(argc, argv);

  if ( status != 0 ) return -1;

  cdo_set_options();

  if ( Debug || Version ) cdo_version();

  if ( Debug )
    {
      fprintf(stderr, "stdin_is_tty:   %d\n", stdin_is_tty);
      fprintf(stderr, "stdout_is_tty:  %d\n", stdout_is_tty);
      fprintf(stderr, "stderr_is_tty:  %d\n", stderr_is_tty);
    }

  if ( Debug ) print_system_info();

  check_stacksize();

  if ( Debug ) print_pthread_info();

  if ( Debug )
    {
      //      fprintf(stderr, "C++ max thread      = %u\n", std::thread::hardware_concurrency());
    }

#if defined(_OPENMP)
  if ( numThreads <= 0 ) numThreads = 1;
  omp_set_num_threads(numThreads);

  if ( Debug )
    {
      fprintf(stderr, "OMP num procs       = %d\n", omp_get_num_procs());
      fprintf(stderr, "OMP max threads     = %d\n", omp_get_max_threads());
      fprintf(stderr, "OMP num threads     = %d\n", omp_get_num_threads());
#if defined(HAVE_OPENMP3)
      fprintf(stderr, "OMP thread limit    = %d\n", omp_get_thread_limit());
      omp_sched_t kind;
      int modifer;
      omp_get_schedule(&kind, &modifer);
      fprintf(stderr, "OMP schedule        = %d (1:static; 2:dynamic; 3:guided; 4:auto)\n", (int) kind);
#endif
#if defined(HAVE_OPENMP4)
      fprintf(stderr, "OMP proc bind       = %d (0:false; 1:true; 2:master; 3:close; 4:spread)\n", (int) omp_get_proc_bind());
#if !defined(__ICC)
      fprintf(stderr, "OMP num devices     = %d\n", omp_get_num_devices());
#endif
#endif
    }

  ompNumThreads = omp_get_max_threads();
  if ( omp_get_max_threads() > omp_get_num_procs() )
    fprintf(stderr, "Warning: Number of OMP threads is greater than number of Cores=%d!\n", omp_get_num_procs());
  if ( ompNumThreads < numThreads )
    fprintf(stderr, "Warning: omp_get_max_threads() returns %d!\n", ompNumThreads);
  if ( cdoVerbose )
    {
      fprintf(stderr, " OpenMP:  num_procs = %d  max_threads = %d", omp_get_num_procs(), omp_get_max_threads());
#if defined(HAVE_OPENMP4)
#if !defined(__ICC)
      fprintf(stderr, "  num_devices = %d", omp_get_num_devices());
#endif
#endif
      fprintf(stderr, "\n");
    }
#else
  if ( numThreads > 0 )
    {
      fprintf(stderr, "Option -P failed, OpenMP support not compiled in!\n");
      return -1;
    }
#endif


  if ( CDO_optind < argc )
    {
      operatorArg = argv[CDO_optind];
      argument = argument_new(argc-CDO_optind, 0);
      argument_fill(argument, argc-CDO_optind, &argv[CDO_optind]);
    }
  else
    {
      if ( ! Version && ! Help )
        {
          fprintf(stderr, "\nNo operator given!\n\n");
          cdo_usage();
          status = 1;
        }

      if ( Help ) cdo_usage();
      lstop = TRUE;
    }

  if ( lstop ) return status;

  if ( cdoDefaultTableID != CDI_UNDEFID ) cdiDefTableID(cdoDefaultTableID);

  extern int (*proj_lonlat_to_lcc_func)();
  proj_lonlat_to_lcc_func = (int (*)()) proj_lonlat_to_lcc;
  extern int (*proj_lcc_to_lonlat_func)();
  proj_lcc_to_lonlat_func = (int (*)()) proj_lcc_to_lonlat;

  const char *operatorName = getOperatorName(operatorArg);

  if ( Help )
    {
      cdoPrintHelp(operatorHelp(operatorName));
    }
  else if ( cdoExpMode == CDO_EXP_LOCAL )
    {
      exp_run(argc, argv, cdoExpName);
    }
  else
    {
      if ( Debug )
        {
          if ( DebugLevel == 0 ) DebugLevel = 1;
          cdoSetDebug(DebugLevel);
        }

      timer_total  = timer_new("total");
      timer_read   = timer_new("read");
      timer_write  = timer_new("write");

      timer_start(timer_total);

#ifdef CUSTOM_MODULES
      load_custom_modules("custom_modules");
      operatorModule(operatorName)(argument);
      close_library_handles();
#else
      operatorModule(operatorName)(argument);

#endif

      timer_stop(timer_total);

      if ( cdoTimer ) timer_report();
    }

  if ( argument ) argument_free(argument);

  if ( cdoVarnames )
    {
      if ( cdoNumVarnames ) Free(cdoVarnames[0]);
      Free(cdoVarnames);
    }

  /* problems with alias!!! if ( operatorName ) Free(operatorName); */ 

  /* malloc_stats(); */

  if ( cdoGridSearchDir ) Free(cdoGridSearchDir);

  if ( CDO_Rusage ) cdo_rusage();

  return status;
}
