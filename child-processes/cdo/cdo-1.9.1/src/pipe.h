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

#ifndef PIPE_H
#define PIPE_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdbool.h>
#include <sys/types.h>

#if defined(HAVE_LIBPTHREAD)

#include <pthread.h>
#include "pthread_debug.h"
//#include "pstream.h"
#include <iostream>

#endif

#if defined(HAVE_LIBPTHREAD)

#include <condition_variable>
#include <mutex>

struct pipe_t
{

public:
  pipe_t();
  int pipeInqVlist(int &vlistID);
  void pipe_init();
  void pipeDefRecord(int p_varId, int p_levelID);
  void pipeDefTimestep(int p_vlistID, int tsID);
  void pipeDefVlist(int &target_vlistID, int new_vlistID);

  int pipeInqTimestep(int p_tsID);
  int pipeInqRecord(int *varID, int *levelID);

  void pipeWriteRecord(double *p_data, int p_nmiss);
  void pipeReadRecord(int p_vlistID, double *data, int *nmiss);
  void pipeReadPipeRecord(double *data, int vlistID, int *p_nmiss);

  bool EOP;
  bool usedata;
  bool hasdata;
  int nrecs;
  int varID, levelID;
  int recIDr, recIDw, tsIDr, tsIDw;
  int nmiss;
  double *data;
  // pstream_t *pstreamptr_in;
  /* unsigned long */ off_t nvals;

  std::mutex m_mutex;
  std::condition_variable tsDef, tsInq, vlistDef, isclosed;
  std::condition_variable recDef, recInq;
  std::condition_variable writeCond, readCond;

  std::string name;
};

void pipeDebug(int debug);

#endif

#endif /* PIPE_H */
