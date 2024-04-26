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

#ifndef _PROCESS_H
#define _PROCESS_H

#include <sys/types.h> /* off_t */
#include <vector>
#include "util.h"
#include "pstream.h"
#include "modules.h"

#include <vector>
#include <iostream>

constexpr int MAX_PROCESS  =   128;
constexpr int MAX_STREAM   =    64;
constexpr int MAX_OPERATOR =   128;
constexpr int MAX_OARGC    =  4096;
constexpr int MAX_FILES    = 65536;


typedef struct {
  int         f1;
  int         f2;
  const char *name;
  const char *enter;
}
oper_t;

class process_t {
    public:
   int m_ID;
#if defined(HAVE_LIBPTHREAD)
  pthread_t   threadID;
  int         l_threadID;
#endif
  short       nchild;
  std::vector<pstream_t*>       inputStreams;
  std::vector<pstream_t*>       outputStreams;
  double      s_utime;
  double      s_stime;
  double      a_utime;
  double      a_stime;
  double      cputime;

  off_t       nvals;
  short       nvars;
  int         ntimesteps;
  short       streamCnt;
  std::vector<argument_t> streamNames;
  char       *xoperator;
  const char *operatorName;
  char       *operatorArg;
  int         oargc;
  std::vector<char *> oargv;
  char        prompt[64];
  short       noper;
  oper_t      oper[MAX_OPERATOR];

  modules_t module;

  int getInStreamCnt();
  int getOutStreamCnt();
  void initProcess();
  void print_process();
  process_t(int ID);
 private: 
  process_t();
  void OpenRead(int p_input_idx);
  void OpenWrite(int p_input_idx);
  void OpenAppend(int p_input_idx);
};

  pstream_t*  processInqInputStream(int streamindex);
  pstream_t*  processInqOutputStream(int streamindex);
  process_t&  processSelf(void);
int  processCreate(void);
void processDelete(void);
int  processInqTimesteps(void);
void processDefTimesteps(int streamID);
int  processInqVarNum(void);
int processInqInputStreamNum(void);
int processInqOutputStreamNum(void);
void processAddInputStream(pstream_t *p_pstream_ptr);
void processAddOutputStream(pstream_t *p_pstream_ptr);
void processDelStream(int streamID);
void processDefVarNum(int nvars);
void processDefArgument(void *vargument);

void processStartTime(double *utime, double *stime);
void processEndTime(double *utime, double *stime);
void processAccuTime(double utime, double stime);

void processDefCputime(int processID, double cputime);
double processInqCputime(int processID);

void processAddNvals(off_t nvals);
off_t processInqNvals(int processID);
int processNums(void);

int  processInqChildNum(void);

const char *processOperatorArg(void);
const char *processInqOpername(void);
const char *processInqOpername2(int processID);
const char *processInqPrompt(void);

const argument_t *cdoStreamName(int cnt);

#endif  /* _PROCESS_H */
