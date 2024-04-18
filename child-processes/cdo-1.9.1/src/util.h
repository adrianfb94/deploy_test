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

#ifndef _UTIL_H
#define _UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include "percentiles.h"

/* dummy use of unused parameters to silence compiler warnings */
#define  UNUSED(x) (void)x

#undef   TRUE
#define  TRUE   1
#undef   FALSE
#define  FALSE  0

#undef   MIN
#define  MIN(a,b)  ((a) < (b) ? (a) : (b))
#undef   MAX
#define  MAX(a,b)  ((a) > (b) ? (a) : (b))

#undef   SQR
#define  SQR(a)    ((a)*(a))


#define  ADD_PLURAL(n)  ((n)!=1 ? "s" : "")

#define  UNCHANGED_RECORD  (processSelf().m_ID == 0 && cdoStreamName(0)->argv[0][0] != '-' && cdoRegulargrid == FALSE && cdoDefaultFileType == -1 && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )

#include <string>
extern char *Progname;
extern char *cdoGridSearchDir;
extern int CDO_Reduce_Dim;
extern int CDO_Memtype;
extern int CDO_Parallel_Read;
extern int CDO_Append_History;
extern int CDO_Reset_History;
extern int timer_read, timer_write; // refactor: both pstream.cc and CDIread.cc CDIwrite.cc defined in cdo.cc

extern int CDO_optind;
extern const char *CDO_optarg;
extern int CDO_opterr;

extern int CDO_flt_digits;
extern int CDO_dbl_digits;

extern int remap_genweights;

extern const char *cdoExpName;
extern int ompNumThreads;

extern int stdin_is_tty;
extern int stdout_is_tty;
extern int stderr_is_tty;

extern int cdoDefaultFileType;
extern int cdoDefaultDataType;
extern int cdoDefaultByteorder;
extern int cdoDefaultTableID;
extern int cdoDefaultInstID;
extern int cdoDefaultTimeType;
extern int cdoLogOff;

extern int cdoLockIO;
extern int cdoCheckDatarange;

extern int cdoSilentMode;
extern int cdoOverwriteMode;
extern int cdoRegulargrid;
extern int cdoBenchmark;
extern int cdoTimer;
extern int cdoVerbose;
extern int cdoDebug;
extern int cdoCompress;
extern int cdoInteractive;
extern int cdoParIO;
extern int cdoDebugExt;

extern int cdoCompType;
extern int cdoCompLevel;

extern int cdoChunkType;

extern int cdoExpMode;

extern int CDO_Color;
extern int CDO_Use_FFTW;
extern int CDO_Version_Info;
extern int CDO_CMOR_Mode;
extern int cdoDiag;

extern int cdoNumVarnames;
extern char **cdoVarnames;
extern char CDO_File_Suffix[32]; // refactor: added keyword extern

extern const char *CDO_Version;




char *getProgname(char *string);
char *GetOperator(const char *argument);
const char *getOperatorName(const char *xoperator);
char *getOperatorArg(const char *xoperator);
const char *cdoComment(void);

char *getFileArg(char *argument);

enum {START_DEC, START_JAN};
int get_season_start(void);
void get_season_name(const char *seas_name[]);
int month_to_season(int month);

void init_is_tty(void);

char *double_to_attstr(int digits, char *str, size_t len, double value);

void progressInit(void);
void progressStatus(double offset, double refval, double curval);

bool fileExists(const char *filename);
bool userFileOverwrite(const char *filename);

/* convert a CDI datatype to string */
int datatype2str(int datatype, char *datatypestr);
int str2datatype(const char *datatypestr);

/* filename manipulation */
const char *filetypeext(int filetype);
void rm_filetypeext(char *file, const char *ext);
void repl_filetypeext(char file[], const char *oldext, const char *newext);


/* moved here from cdo.h */
void    cdiOpenError(int cdiErrno, const char *fmt, const char *path);
void    cdoAbort(const char *fmt, ...);
void    cdoWarning(const char *fmt, ...);
void    cdoPrint(const char *fmt, ...);
void    cdoPrintBlue(const char *fmt, ...);
void    cdoPrintRed(const char *fmt, ...);

int  timer_new(const char *text);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);
double timer_val(int it);


void    operatorInputArg(const char *enter);
int     operatorArgc(void);
char  **operatorArgv(void);
void    operatorCheckArgc(int numargs);


void    cdoInitialize(void *argument);
void    cdoFinish(void);

int     cdoStreamNumber(void);
int     cdoStreamCnt(void);
int     cdoOperatorAdd(const char *name, int func, int intval, const char *enter);
int     cdoOperatorID(void);
int     cdoOperatorF1(int operID);
int     cdoOperatorF2(int operID);
const char *cdoOperatorName(int operID);
const char *cdoOperatorEnter(int operID);

int     cdoFiletype(void);

void cdoSetNAN(double missval, size_t gridsize, double *array);

void    cdoInqHistory(int fileID);
void    cdoDefHistory(int fileID, char *histstring);
void cdo_def_tracking_id(int vlistID, const char *uuid_attribute);
void cdo_def_creation_date(int vlistID);

int     cdoDefineGrid(const char *gridfile);
int     cdoDefineZaxis(const char *zaxisfile);

int     vlistInqNWPV(int vlistID, int varID);
int     vlistIsSzipped(int vlistID);
int     vlist_check_gridsize(int vlistID);
int     vlist_get_psvarid(int vlistID, int zaxisID);
double *vlist_read_vct(int vlistID, int *rzaxisIDh, int *rnvct, int *rnhlev, int *rnhlevf, int *rnhlevh);
void vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2);

void cdoGenFileSuffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname);

void writeNCgrid(const char *gridfile, int gridID, int *imask);
void defineZaxis(const char *zaxisarg);

int grid_from_name(const char *gridname);
int zaxisFromName(const char *zaxisname);

/* refactor: moved here from cdo.h */
int cdo_omp_get_thread_num(void);
void cdo_omp_set_num_threads(int nthreads);
std::string string2lower(std::string str);
void strtolower(char *str);
void strtoupper(char *str);

/* refactor: moved here from cdo.cc */
void exp_run(int argc, char *argv[], const char *cdoExpName); // job.cc
void printFeatures(void); // features.cc
void printLibraries(void);  // features.cc  

int wildcardmatch(const char *w, const char *s);

void cdo_check_round(void);

#endif  /* _UTIL_H */
