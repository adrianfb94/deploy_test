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

#include <stdio.h>
#include <string.h>
#include <time.h>  // time()
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "dmemory.h"
#include "pipe.h"

#if defined(HAVE_LIBPTHREAD)

static int PipeDebug = 0;

pipe_t::pipe_t() { pipe_init(); }

void
pipe_t::pipe_init()
{
/*  pthread_mutexattr_t m_attr;
  pthread_condattr_t c_attr;

  pthread_mutexattr_init(&m_attr);
  pthread_condattr_init(&c_attr);
  */
  /*
#if defined(_POSIX_THREAD_PROCESS_SHARED)
  if ( PipeDebug )
    {
      Message("setpshared mutexattr to PTHREAD_PROCESS_SHARED");
      Message("setpshared condattr to PTHREAD_PROCESS_SHARED");
    }

  pthread_mutexattr_setpshared(&m_attr, PTHREAD_PROCESS_SHARED);
  pthread_condattr_setpshared(&c_attr, PTHREAD_PROCESS_SHARED);

  if ( PipeDebug )
    {
      int pshared;
      pthread_mutexattr_getpshared(&m_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
        Message("getpshared mutexattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
        Message("getpshared mutexattr is PTHREAD_PROCESS_PRIVATE");

      pthread_condattr_getpshared(&c_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
        Message("getpshared condattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
        Message("getpshared condattr is PTHREAD_PROCESS_PRIVATE");
    }
#else
  if ( PipeDebug )
    Message("_POSIX_THREAD_PROCESS_SHARED undefined");
#endif
  */
  EOP = false;

  recIDr = -1;
  recIDw = -1;
  tsIDr = -1;
  tsIDw = -1;

  nvals = 0;
  nmiss = 0;
  data = NULL;
  hasdata = false;
  usedata = true;
  // pstreamptr_in = 0;
/*
  mutex = (pthread_mutex_t *) Malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(mutex, &m_attr);

  tsDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(tsDef, &c_attr);
  tsInq = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(tsInq, &c_attr);

  recDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(recDef, &c_attr);
  recInq = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(recInq, &c_attr);

  vlistDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(vlistDef, &c_attr);
  isclosed = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(isclosed, &c_attr);

  writeCond = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(writeCond, &c_attr);

  readCond = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(readCond, &c_attr);

  pthread_mutexattr_destroy(&m_attr);
  pthread_condattr_destroy(&c_attr);
  */
}
int
pipe_t::pipeInqTimestep(int p_tsID)
{
  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  usedata = false;
  recIDr = -1;
  if (p_tsID != tsIDr + 1)
    {
      if (!(p_tsID == tsIDr && tsIDr == tsIDw && recIDr == -1))
        Error("%s unexpected tsID %d %d %d", name.c_str(), p_tsID, tsIDr + 1, tsIDw);
    }

  tsIDr = p_tsID;
  while (tsIDw != p_tsID)
    {
      if (EOP)
        {
          if (PipeDebug)
            Message("%s EOP", name.c_str());
          break;
        }
      if (hasdata)
        {
          if (PipeDebug)
            Message("%s has data", name.c_str());
          hasdata = false;
          data = NULL;
          readCond.notify_all();
        }
      else if (PipeDebug)
        Message("%s has no data", name.c_str());

      recInq.notify_all(); /* o.k. ??? */

      if (PipeDebug)
        Message("%s wait of tsDef", name.c_str());
      tsDef.wait(locked_mutex);
    }

  int numrecs = EOP ? 0 : nrecs;

  locked_mutex.unlock();
  // UNLOCK

  tsInq.notify_all();

  return numrecs;
}

void
pipe_t::pipeDefVlist(int &target_vlistID, int new_vlistID)
{

  // LOCK
  m_mutex.lock();
  target_vlistID = new_vlistID;
  m_mutex.unlock();
  // UNLOCK

  // lets the program know that the vlist is now defined
  vlistDef.notify_all();
}

//#define TIMEOUT 1  // wait 1 seconds
constexpr std::chrono::milliseconds TIMEOUT = std::chrono::milliseconds(1000);
#define MIN_WAIT_CYCLES 10
#define MAX_WAIT_CYCLES 3600
int processNumsActive(void);

int
pipe_t::pipeInqVlist(int &p_vlistID)
{
  int vlistID = -1;
  std::chrono::milliseconds time_to_wait(0);
  std::cv_status retcode = std::cv_status::timeout;
  int nwaitcycles = 0;

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  while (p_vlistID == -1 && retcode == std::cv_status::timeout)
    {
      time_to_wait += TIMEOUT;
      // fprintf(stderr, "tvsec %g\n", (double) time_to_wait.tv_sec);
      if (PipeDebug)
        Message("%s wait of vlistDef", name.c_str());
      // pthread_cond_wait(pipe->vlistDef, pipe->mutex);
      retcode = vlistDef.wait_for(locked_mutex, time_to_wait);
      // fprintf(stderr, "self %d retcode %d %d %d\n", pstreamptr->self, retcode, processNumsActive(),
      // vlistID);
      //if (retcode != 0 && nwaitcycles++ < MAX_WAIT_CYCLES)
      if (retcode != std::cv_status::timeout && nwaitcycles++ < MAX_WAIT_CYCLES)
        {
          if (processNumsActive() > 1 || (processNumsActive() == 1 && nwaitcycles < MIN_WAIT_CYCLES))
            retcode = std::cv_status::timeout;
        }
    }

  if (retcode == std::cv_status::timeout)
    vlistID = p_vlistID;
  else if (PipeDebug)
    Message("%s timeout!", name.c_str());

  // UNLOCK

  return vlistID;
}

void
pipe_t::pipeDefTimestep(int p_vlistID, int p_tsID)
{
  int numrecs;

  // LOCK
  m_mutex.lock();
  recIDw = -1;
  tsIDw++;
  if (p_tsID != tsIDw)
    Error("unexpected p_tsID %d(%d) for %s", p_tsID, tsIDw, name.c_str());

  if (p_tsID == 0)
    numrecs = vlistNrecs(p_vlistID);
  else
    {
      int vlistID = p_vlistID;
      numrecs = 0;
      for (int varID = 0; varID < vlistNvars(vlistID); varID++)
        {
          if (vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT)
            {
              numrecs += zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
            }
        }
      if (PipeDebug)
        Message("%s numrecs = %d nvars = %d", name.c_str(), numrecs, vlistNvars(vlistID));
    }

  nrecs = numrecs;
  if (PipeDebug)
    Message("%s numrecs %d p_tsID %d %d %d", name.c_str(), numrecs, p_tsID, tsIDw, tsIDr);
  if (numrecs == 0)
    EOP = true;
  m_mutex.unlock();
  // UNLOCK

  tsDef.notify_all();
  // sleep(1);

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  while (tsIDr < p_tsID)
    {
      if (EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of tsInq (p_tsID %d %d)", name.c_str(), p_tsID, tsIDr);
      tsInq.wait(locked_mutex);
    }
  // UNLOCK
}

int
pipe_t::pipeInqRecord(int *p_varID, int *p_levelID)
{
  bool condSignal = false;

  // if (PipeDebug)

  // LOCK
  m_mutex.lock();
  if (PipeDebug)
    Message("%s has no data %d %d", name.c_str(), recIDr, recIDw);
  if (hasdata || usedata)
    {
      hasdata = false;
      data = NULL;
      usedata = false;
      condSignal = true;
    }
  m_mutex.unlock();
  // UNLOCK

  if (condSignal)
    readCond.notify_all();

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  usedata = true;
  recIDr++;

  if (PipeDebug)
    Message("%s recID %d %d", name.c_str(), recIDr, recIDw);

  while (recIDw != recIDr)
    {
      if (EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of recDef", name.c_str());
      recDef.wait(locked_mutex);
    }

  if (EOP)
    {
      *p_varID = -1;
      *p_levelID = -1;
    }
  else
    {
      *p_varID = varID;
      *p_levelID = levelID;
    }

  locked_mutex.unlock();
  // UNLOCK

  recInq.notify_all();

  return 0;
}

void
pipe_t::pipeDefRecord(int p_varID, int p_levelID)
{
  bool condSignal = false;

  // LOCK
  m_mutex.lock();
  if (PipeDebug)
    Message("%s has data %d %d", name.c_str(), recIDr, recIDw);
  if (hasdata)
    {
      hasdata = false;
      data = NULL;
      condSignal = true;
    }
  m_mutex.unlock();
  // UNLOCK

  if (condSignal)
    readCond.notify_all();

  // LOCK
  m_mutex.lock();
  usedata = true;
  recIDw++;
  varID = p_varID;
  levelID = p_levelID;
  if (PipeDebug)
    Message("%s recID %d %d", name.c_str(), recIDr, recIDw);
  m_mutex.unlock();
  // UNLOCK

  recDef.notify_all();

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  while (recIDr < recIDw)
    {
      if (tsIDw != tsIDr)
        break;
      if (EOP)
        break;
      if (PipeDebug)
        Message("%s wait of recInq %d", name.c_str(), recIDr);
      recInq.wait(locked_mutex);
    }
  // UNLOCK
}

/***
 * copys data from a pipe to data
 *
 * @param data destination for the record data
 * @param pipe pipe that has the wanted data
 */
void
pipe_t::pipeReadPipeRecord(double *p_data, int vlistID, int *p_nmiss)
{
  if (!p_data)
    Error("No data pointer for %s", name.c_str());

  int datasize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  nvals += datasize;
  if (vlistNumber(vlistID) != CDI_REAL)
    datasize *= 2;
  memcpy(p_data, data, datasize * sizeof(double));
  *p_nmiss = nmiss;
}

/*
void
pipeGetReadTarget(pstream_t *pstreamptr, pstream_t *pstreamptr_in)
{

  pstreamptr_in = pstreamptr->pipe->pstreamptr_in;
  pstreamptr = pstreamptr_in;
  while (pstreamptr_in->ispipe)
    {
      if (PipeDebug)
        fprintf(stderr, "%s: istream %d is pipe\n", __func__, pstreamptr_in->self);
      pstreamptr = pstreamptr_in;
      pstreamptr_in = pstreamptr_in->pipe->pstreamptr_in;
      if (pstreamptr_in == 0)
        break;
    }

  if (pstreamptr->pipe->hasdata == false)
    {
      Error("Internal problem! istream undefined");
    }
  else if (!pstreamptr->pipe->data)
    {
      Error("No data pointer for %s", pstreamptr->pipe->name.c_str());
    }
}
*/
void
pipe_t::pipeReadRecord(int p_vlistID, double *data, int *nmiss)
{
  *nmiss = 0;

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  while (!hasdata)
    {
      if (PipeDebug)
        Message("%s wait of writeCond", name.c_str());
      writeCond.wait(locked_mutex);
    }

  if (hasdata)
    {
      pipeReadPipeRecord(data, p_vlistID, nmiss);
    }
  else
    {
      Error("data type %d not implemented", hasdata);
    }

  if (PipeDebug)
    Message("%s read record %d", name.c_str(), recIDr);

  hasdata = false;
  data = NULL;
  locked_mutex.unlock();
  // UNLOCK

  readCond.notify_all();
}

void
pipe_t::pipeWriteRecord(double *p_data, int p_nmiss)
{
  /*
  if ( ! usedata ) return;
  */
  // LOCK
  m_mutex.lock();
  hasdata = true; /* data pointer */
  data = p_data;
  nmiss = p_nmiss;
  m_mutex.unlock();
  // UNLOCK

  writeCond.notify_all();

  if (PipeDebug)
    Message("%s write record %d", name.c_str(), recIDw);

  // LOCK
  std::unique_lock<std::mutex> locked_mutex(m_mutex);
  while (hasdata)
    {
      if (!usedata)
        break;
      /*
      printf("ts ids %d %d\n", tsIDw, tsIDr);
      printf("rec ids %d %d\n", recIDw, recIDr);
      */
      if (recIDw != recIDr)
        break;

      if (EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of readCond", name.c_str());
      readCond.wait(locked_mutex);
    }
  // UNLOCK
}

void
pipeDebug(int debug)
{
  PipeDebug = debug;
}

#endif
