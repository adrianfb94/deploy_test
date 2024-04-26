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


#include <thread>
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h> /* stat */

FILE *popen(const char *command, const char *type);
int pclose(FILE *stream);

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "pstream.h"
#include "pstream_int.h"
#include "util.h"
#include "pipe.h"
#include "error.h"

static int PSTREAM_Debug = 0;

//#define MAX_PSTREAMS 4096

//static int _pstream_max = MAX_PSTREAMS;


//static void pstream_initialize(void);

//static bool _pstream_init = false;

#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>
#include "pthread_debug.h"

// TODO: make threadsafe
static int pthreadScope = 0;

static pthread_mutex_t streamOpenReadMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamOpenWriteMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamMutex = PTHREAD_MUTEX_INITIALIZER;

//static pthread_once_t _pstream_init_thread = PTHREAD_ONCE_INIT;
//static pthread_mutex_t _pstream_mutex;

static std::mutex _pstream_map_mutex;
#define PSTREAM_LOCK() _pstream_map_mutex.lock();
#define PSTREAM_UNLOCK() _pstream_map_mutex.unlock();
/*
#define PSTREAM_LOCK() pthread_mutex_lock(&_pstream_mutex)
#define PSTREAM_UNLOCK() pthread_mutex_unlock(&_pstream_mutex)
#define PSTREAM_INIT() \
  if (!_pstream_init)  \
  pthread_once(&_pstream_init_thread, pstream_initialize)

*/
#else

#define PSTREAM_LOCK()
#define PSTREAM_UNLOCK()
#define PSTREAM_INIT() \
  if (!_pstream_init)  \
  pstream_initialize()

#endif

/*
typedef struct _pstreamPtrToIdx
{
  int idx;
  pstream_t *ptr;
  struct _pstreamPtrToIdx *next;
} pstreamPtrToIdx;
*/
/*
static pstreamPtrToIdx *_pstreamList = NULL;
static pstreamPtrToIdx *_pstreamAvail = NULL;
*/
static std::map<int,pstream_t> _pstream_map;
static int next_pstream_id = 1;
/*
static void
pstream_list_new(void)
{
  assert(_pstreamList == NULL);

  _pstreamList = (pstreamPtrToIdx *) Malloc(_pstream_max * sizeof(pstreamPtrToIdx));
}
static void
pstream_list_delete(void)
{
  if (_pstreamList)
    Free(_pstreamList);
}

static void
pstream_init_pointer(void)
{
  for (int i = 0; i < _pstream_max; ++i)
    {
      _pstreamList[i].next = _pstreamList + i + 1;
      _pstreamList[i].idx = i;
      _pstreamList[i].ptr = 0;
    }

  _pstreamList[_pstream_max - 1].next = 0;

  _pstreamAvail = _pstreamList;
}
*/

static pstream_t *create_pstream()
{
    PSTREAM_LOCK();
    auto new_entry  = _pstream_map.insert(
            std::make_pair(next_pstream_id, pstream_t(next_pstream_id))
            );
    next_pstream_id++;
    PSTREAM_UNLOCK();

    return &new_entry.first->second;
}

static pstream_t * pstream_to_pointer(int idx)
{
    PSTREAM_LOCK();
    auto pstream_iterator = _pstream_map.find(idx);
    PSTREAM_UNLOCK();
    if(pstream_iterator == _pstream_map.end())
    {
        Error("pstream index %d undefined!", idx);
    }

    return &pstream_iterator->second;
}
/*
static pstream_t *
pstream_to_pointer(int idx)
{
  pstream_t *pstreamptr = NULL;

  PSTREAM_INIT();

  if (idx >= 0 && idx < _pstream_max)
    //{
      PSTREAM_LOCK();

      pstreamptr = _pstreamList[idx].ptr;

      PSTREAM_UNLOCK();
    }
  else
    Error("pstream index %d undefined!", idx);

  return pstreamptr;
}
*/
/* Create an index from a pointer */
/*
static int
pstream_from_pointer(pstream_t *ptr)
{
  int idx = -1;

  if (ptr)
    {
      PSTREAM_LOCK();

      if (_pstreamAvail)
        {
          pstreamPtrToIdx *newptr = _pstreamAvail;
          _pstreamAvail = _pstreamAvail->next;
          newptr->next = 0;
          idx = newptr->idx;
          newptr->ptr = ptr;

          if (PSTREAM_Debug)
            Message("Pointer %p has idx %d from pstream list", ptr, idx);
        }
      else
        Error("Too many open pstreams (limit is %d)!", _pstream_max);

      PSTREAM_UNLOCK();
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return idx;
}
*/

void pstream_t::init()
{
  isopen = true;
  ispipe = false;
  m_fileID = -1;
  m_vlistID = -1;
  tsID = -1;
  m_filetype = -1;
  m_name = "";
  tsID0 = 0;
  mfiles = 0;
  nfiles = 0;
  varID = -1;
  m_varlist = NULL;
#if defined(HAVE_LIBPTHREAD)
  argument = NULL;
  pipe = NULL;
//  pstreamptr->rthreadID  = 0;
//  pstreamptr->wthreadID  = 0;
#endif
}
pstream_t::pstream_t(int p_id) : self(p_id) { init(); }
pstream_t::~pstream_t(){
  //vlistDestroy(m_vlistID);
}

static void
pstream_delete_entry(pstream_t *pstreamptr)
{
  int idx = pstreamptr->self;

  PSTREAM_LOCK();

  _pstream_map.erase(idx);

  PSTREAM_UNLOCK();

  if (PSTREAM_Debug)
    Message("Removed idx %d from pstream list", idx);
}
/*
static void
pstream_initialize(void)
{
#if defined(HAVE_LIBPTHREAD)
  // initialize global API mutex lock 
  pthread_mutex_init(&_pstream_mutex, NULL);
#endif

  char *env = getenv("PSTREAM_DEBUG");
  if (env)
    PSTREAM_Debug = atoi(env);

  env = getenv("PSTREAM_MAX");
  if (env)
    _pstream_max = atoi(env);

  if (PSTREAM_Debug)
    Message("PSTREAM_MAX = %d", _pstream_max);

  pstream_list_new();
  atexit(pstream_list_delete);

  pstream_init_pointer();

  _pstream_init = true;
}
*/
static int pstreamFindID(const char *p_name)
{
    std::string cur_name;
    for(auto map_pair :  _pstream_map)
    {
        cur_name = map_pair.second.m_name;
        if(!(cur_name.empty())){
            if(cur_name.compare(p_name) == 0)
            {
                return map_pair.first;
            }
        }
    }
    return -1;
}
/*
static int
pstreamFindID(const char *name)
{
  pstream_t *pstreamptr;
  int pstreamID;

  for (pstreamID = 0; pstreamID < _pstream_max; ++pstreamID)
    {
      pstreamptr = pstream_to_pointer(pstreamID);

      if (pstreamptr)
        if (pstreamptr->m_name)
          if (strcmp(pstreamptr->m_name, m_name) == 0)
            break;
    }

  if (pstreamID == _pstream_max)
    pstreamID = -1;

  return pstreamID;
}
*/
bool
pstream_t::isPipe()
{
  return ispipe;
}

pthread_t
pCreateReadThread(argument_t *argument)
{
  pthread_attr_t attr;
  int status = pthread_attr_init(&attr);
  if (status)
    SysError("pthread_attr_init failed for '%s'", argument->operatorName.c_str());
  status = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (status)
    SysError("pthread_attr_setdetachstate failed for '%s'", argument->operatorName.c_str());
  /*
    param.sched_priority = 0;
    status = pthread_attr_setschedparam(&attr, &param);
    if ( status ) SysError("pthread_attr_setschedparam failed for '%s'", newarg+1);
  */
  /* status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED); */
  /* if ( status ) SysError("pthread_attr_setinheritsched failed for '%s'", newarg+1); */

  pthread_attr_getscope(&attr, &pthreadScope);

  /* status = pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS); */
  /* if ( status ) SysError("pthread_attr_setscope failed for '%s'", newarg+1); */
  /* If system scheduling scope is specified, then the thread is scheduled against all threads in the system */
  /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */

  size_t stacksize = 0;
  status = pthread_attr_getstacksize(&attr, &stacksize);
  if (stacksize < 2097152)
    {
      stacksize = 2097152;
      pthread_attr_setstacksize(&attr, stacksize);
    }

  pthread_t thrID;
  int rval = pthread_create(&thrID, &attr, operatorModule(argument->operatorName.c_str()), argument);
  if (rval != 0)
    {
      errno = rval;
      SysError("pthread_create failed for '%s'", argument->operatorName.c_str());
    }
  return thrID;
}

void
pstream_t::pstreamOpenReadPipe(const char *pipename)
{
#if defined(HAVE_LIBPTHREAD)
  // int pstreamID = pstreamptr->self;

  ispipe = true;
  m_name = pipename;
  rthreadID = pthread_self();
  pipe = new pipe_t();
  pipe->name = std::string(pipename);

    /* Free(operatorName); */
  /*      pipeInqInfo(pstreamID); */
  if (PSTREAM_Debug)
    Message("pipe %s", pipename);
#else
  cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
}

void pstream_t::createFilelist(const char * p_args)
{
  size_t i;
  size_t len = strlen(p_args);

  for (i = 0; i < len; i++)
    if (p_args[i] == ':')
      break;

  if (i < len)
    {
      int nfiles = 1, j;

      const char *pch = &p_args[i + 1];
      len -= (i + 1);
      if (len && (strncmp(p_args, "filelist:", 9) == 0 || strncmp(p_args, "flist:", 6) == 0))
        {
          for (i = 0; i < len; i++)
            if (pch[i] == ',')
              nfiles++;

          if (nfiles == 1)
            {
              char line[4096];
              FILE *fp, *fp2;
              fp = fopen(pch, "r");
              if (fp == NULL)
                {
                  cdoAbort("Open failed on %s", pch);
                }
              if (cdoVerbose)
                {
                  cdoPrint("Reading file names from %s", pch);
                }
              /* find number of files */
              nfiles = 0;
              while (readline(fp, line, 4096))
                {
                  if (line[0] == '#' || line[0] == '\0' || line[0] == ' ')
                    continue;

                  fp2 = fopen(line, "r");
                  if (fp2 == NULL)
                    cdoAbort("Open failed on %s", line);
                  fclose(fp2);
                  nfiles++;
                  if (cdoVerbose)
                    cdoPrint("File number %d is %s", nfiles, line);
                }

              if (nfiles == 0)
                cdoAbort("No imput file found in %s", pch);

              mfiles = nfiles;
              m_mfnames.resize(nfiles);
              rewind(fp);

              nfiles = 0;
              while (readline(fp, line, 4096))
                {
                  if (line[0] == '#' || line[0] == '\0' || line[0] == ' ')
                    continue;

                  m_mfnames[nfiles] = line;
                  nfiles++;
                }

              fclose(fp);
            }
          else
            {
              char line[65536];

              mfiles = nfiles;
              m_mfnames.resize(nfiles);

              strcpy(line, pch);
              for (i = 0; i < len; i++)
                {
                  if (line[i] == ',')
                    {
                      line[i] = 0;
                    }
                }
              i = 0;
              for (j = 0; j < nfiles; j++)
                {
                  m_mfnames[j] = line[i];
                  i += strlen(&line[i]) + 1;
                }
            }
        }
      else if (len && strncmp(p_args, "ls:", 3) == 0)
        {
          char line[4096];
          char command[4096];
          char *fnames[16384];
          FILE *pfp;

          strcpy(command, "ls ");
          strcat(command, pch);

          pfp = popen(command, "r");
          if (pfp == 0)
            SysError("popen %s failed", command);

          nfiles = 0;
          while (readline(pfp, line, 4096))
            {
              if (nfiles >= 16384)
                cdoAbort("Too many input files (limit: 16384)");
              fnames[nfiles++] = strdupx(line);
            }

          pclose(pfp);

          mfiles = nfiles;
          m_mfnames.resize(nfiles);

          for (j = 0; j < nfiles; j++)
            m_mfnames[j] = std::string(fnames[j]);
        }
    }
}

void
pstream_t::pstreamOpenReadFile(const char* p_args)
{
  createFilelist(p_args);

  std::string filename; 

  if (mfiles)
    {
      filename = m_mfnames[0];
      nfiles = 1;
    }
  else
    {
     filename = std::string(p_args);
    }

  if (PSTREAM_Debug)
    Message("file %s", filename.c_str());

#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenReadMutex);
#endif
  int fileID = streamOpenRead(filename.c_str());
  if (fileID < 0)
    {
      isopen = false;
      cdiOpenError(fileID, "Open failed on >%s<", filename.c_str());
    }

  if (cdoDefaultFileType == CDI_UNDEFID)
    cdoDefaultFileType = streamInqFiletype(fileID);
  /*
    if ( cdoDefaultInstID == CDI_UNDEFID )
    cdoDefaultInstID = streamInqInstID(fileID);
  */
  cdoInqHistory(fileID);
#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenReadMutex);
#endif

  mode = 'r';
  m_name = filename;
  m_fileID = fileID;
}

void createPipeName(char *pipename, int pnlen)
{

  snprintf(pipename, pnlen, "(pipe%d.%d)", processSelf().m_ID + 1, processInqChildNum() + 1);
}

int
pstreamOpenRead(const argument_t *argument)
{

  pstream_t *pstreamptr = create_pstream();
  if (!pstreamptr)
    Error("No memory");

  int pstreamID = pstreamptr->self;

  int ispipe = argument->args[0] == '-';
  /*
  printf("pstreamOpenRead: args >%s<\n", argument->args);
  for ( int i = 0; i < argument->argc; ++i )
    printf("pstreamOpenRead: arg %d >%s<\n", i, argument->argv[i]);
  */
  if (ispipe)
    {
      size_t pnlen = 16;
      char *pipename = (char *) Malloc(pnlen);
      createPipeName(pipename, pnlen);
      argument_t * newargument = pipe_argument_new(argument, pipename, pnlen);
      pstreamptr->pstreamOpenReadPipe(pipename);
      pCreateReadThread(newargument);
      if (!cdoSilentMode)
      {
        cdoPrint("Started child process \"%s\".", newargument->args + 1);
      }

      processAddInputStream(pstreamptr);
    }
  else
    {
      pstreamptr->pstreamOpenReadFile(argument->args);
    }

  if (pstreamID < 0)
    cdiOpenError(pstreamID, "Open failed on >%s<", argument->args);

  return pstreamID;
}

static void
query_user_exit(const char *argument)
{
/* modified code from NCO */
#define USR_RPL_MAX_LNG 10 /* Maximum length for user reply */
#define USR_RPL_MAX_NBR 10 /* Maximum number of chances for user to reply */
  char usr_rpl[USR_RPL_MAX_LNG];
  int usr_rpl_int;
  short nbr_itr = 0;
  size_t usr_rpl_lng = 0;

  /* Initialize user reply string */
  usr_rpl[0] = 'z';
  usr_rpl[1] = '\0';

  while (!(usr_rpl_lng == 1 && (*usr_rpl == 'o' || *usr_rpl == 'O' || *usr_rpl == 'e' || *usr_rpl == 'E')))
    {
      if (nbr_itr++ > USR_RPL_MAX_NBR)
        {
          (void) fprintf(stdout,
                         "\n%s: ERROR %d failed attempts to obtain valid interactive input.\n",
                         processInqPrompt(),
                         nbr_itr - 1);
          exit(EXIT_FAILURE);
        }

      if (nbr_itr > 1)
        (void) fprintf(stdout, "%s: ERROR Invalid response.\n", processInqPrompt());
      (void) fprintf(stdout,
                     "%s: %s exists ---`e'xit, or `o'verwrite (delete existing file) (e/o)? ",
                     processInqPrompt(),
                     argument);
      (void) fflush(stdout);
      if (fgets(usr_rpl, USR_RPL_MAX_LNG, stdin) == NULL)
        continue;

      /* Ensure last character in input string is \n and replace that with \0 */
      usr_rpl_lng = strlen(usr_rpl);
      if (usr_rpl_lng >= 1)
        if (usr_rpl[usr_rpl_lng - 1] == '\n')
          {
            usr_rpl[usr_rpl_lng - 1] = '\0';
            usr_rpl_lng--;
          }
    }

  /* Ensure one case statement for each exit condition in preceding while loop */
  usr_rpl_int = (int) usr_rpl[0];
  switch (usr_rpl_int)
    {
    case 'E':
    case 'e': exit(EXIT_SUCCESS); break;
    case 'O':
    case 'o': break;
    default: exit(EXIT_FAILURE); break;
    } /* end switch */
}

static int
pstreamOpenWritePipe(const argument_t *argument, int filetype)
{
  int pstreamID = -1;

#if defined(HAVE_LIBPTHREAD)
  if (PSTREAM_Debug)
    {
      Message("pipe %s", argument->args);
    }
  pstreamID = pstreamFindID(argument->args);
  if (pstreamID == -1)
    {
      Error("%s is not open!", argument->args);
    }

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  pstreamptr->wthreadID = pthread_self();
  pstreamptr->m_filetype = filetype;
  processAddOutputStream(pstreamptr);
#endif

  return pstreamID;
}

static void
set_comp(int fileID, int filetype)
{
  if (cdoCompress)
    {
      if (filetype == CDI_FILETYPE_GRB)
        {
          cdoCompType = CDI_COMPRESS_SZIP;
          cdoCompLevel = 0;
        }
      else if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C)
        {
          cdoCompType = CDI_COMPRESS_ZIP;
          cdoCompLevel = 1;
        }
    }

  if (cdoCompType != CDI_COMPRESS_NONE)
    {
      streamDefCompType(fileID, cdoCompType);
      streamDefCompLevel(fileID, cdoCompLevel);

      if (cdoCompType == CDI_COMPRESS_SZIP
          && (filetype != CDI_FILETYPE_GRB && filetype != CDI_FILETYPE_GRB2 && filetype != CDI_FILETYPE_NC4
              && filetype != CDI_FILETYPE_NC4C))
        cdoWarning("SZIP compression not available for non GRIB/NetCDF4 data!");

      if (cdoCompType == CDI_COMPRESS_JPEG && filetype != CDI_FILETYPE_GRB2)
        cdoWarning("JPEG compression not available for non GRIB2 data!");

      if (cdoCompType == CDI_COMPRESS_ZIP && (filetype != CDI_FILETYPE_NC4 && filetype != CDI_FILETYPE_NC4C))
        cdoWarning("Deflate compression not available for non NetCDF4 data!");
    }
}

static int
pstreamOpenWriteFile(const argument_t *argument, int filetype)
{
  char *filename = (char *) Malloc(strlen(argument->args) + 1);

  pstream_t *pstreamptr = create_pstream();
  if (!pstreamptr)
    Error("No memory");

  int pstreamID = pstreamptr->self;

  if (PSTREAM_Debug)
    Message("file %s", argument->args);

  if (filetype == CDI_UNDEFID)
    filetype = CDI_FILETYPE_GRB;

  if (cdoInteractive)
    {
      struct stat stbuf;

      int rstatus = stat(argument->args, &stbuf);
      /* If permanent file already exists, query user whether to overwrite or exit */
      if (rstatus != -1)
        query_user_exit(argument->args);
    }

  if (processNums() == 1 && ompNumThreads == 1)
    timer_start(timer_write);

#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenWriteMutex);
#endif

  int fileID = streamOpenWrite(argument->args, filetype);

#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenWriteMutex);
#endif

  if (processNums() == 1 && ompNumThreads == 1)
    timer_stop(timer_write);
  if (fileID < 0)
    cdiOpenError(fileID, "Open failed on >%s<", argument->args);

  cdoDefHistory(fileID, commandLine());

  if (cdoDefaultByteorder != CDI_UNDEFID)
    streamDefByteorder(fileID, cdoDefaultByteorder);

  set_comp(fileID, filetype);
  /*
    if ( cdoDefaultInstID != CDI_UNDEFID )
    streamDefInstID(fileID, cdoDefaultInstID);
  */
  strcpy(filename, argument->args);

  pstreamptr->mode = 'w';
  pstreamptr->m_name = filename;
  pstreamptr->m_fileID = fileID;
  pstreamptr->m_filetype = filetype;

  return pstreamID;
}

int
pstreamOpenWrite(const argument_t *argument, int filetype)
{
  int pstreamID = -1;

  //PSTREAM_INIT();

  int ispipe = strncmp(argument->args, "(pipe", 5) == 0;

  if (ispipe)
    {
      pstreamID = pstreamOpenWritePipe(argument, filetype);
    }
  else
    {
      pstreamID = pstreamOpenWriteFile(argument, filetype);
    }

  return pstreamID;
}

int
pstreamOpenAppend(const argument_t *argument)
{
  int ispipe = strncmp(argument->args, "(pipe", 5) == 0;

  if (ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("pipe %s", argument->args);
        }
      cdoAbort("this operator doesn't work with pipes!");
    }

  pstream_t *pstreamptr = create_pstream();

  if (!pstreamptr)
    Error("No memory");

  if (PSTREAM_Debug)
    Message("file %s", argument->args);

  pstreamptr->openAppend(argument->args);

  return pstreamptr->self;
}
void
pstream_t::openAppend(const char *p_filename)
{
  if (processNums() == 1 && ompNumThreads == 1)
    {
      timer_start(timer_write);
    }
#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    {
      pthread_mutex_lock(&streamMutex);
    }
  else
    {
      pthread_mutex_lock(&streamOpenReadMutex);
    }
#endif

  int fileID = streamOpenAppend(p_filename);

#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    {
      pthread_mutex_unlock(&streamMutex);
    }
  else
    {
      pthread_mutex_unlock(&streamOpenReadMutex);
    }
#endif
  if (processNums() == 1 && ompNumThreads == 1)
    {
      timer_stop(timer_write);
    }
  if (fileID < 0)
    {
      cdiOpenError(fileID, "Open failed on >%s<", p_filename);
    }
  /*
  cdoInqHistory(fileID);
  cdoDefHistory(fileID, commandLine());
  */
  int filetype = streamInqFiletype(fileID);
  set_comp(fileID, filetype);

  m_name = p_filename;
  mode = 'a';
  m_fileID = fileID;
}
void
pstreamCloseChildStream(pstream_t *pstreamptr)
{
  pipe_t *pipe = pstreamptr->pipe;
  pthread_mutex_lock(pipe->m_mutex);
  pipe->EOP = true;
  if (PSTREAM_Debug)
    Message("%s read closed", pstreamptr->m_name.c_str());
  pthread_mutex_unlock(pipe->m_mutex);
  pthread_cond_signal(pipe->tsDef);
  pthread_cond_signal(pipe->tsInq);

  pthread_cond_signal(pipe->recInq);

  pthread_mutex_lock(pipe->m_mutex);
  pstreamptr->isopen = false;
  pthread_mutex_unlock(pipe->m_mutex);
  pthread_cond_signal(pipe->isclosed);

  pthread_join(pstreamptr->wthreadID, NULL);

  pthread_mutex_lock(pipe->m_mutex);
    if (pstreamptr->argument)
    {
      argument_t *argument = (argument_t *) (pstreamptr->argument);
      if (argument->args)
        Free(argument->args);
      delete (argument);
    }
  pthread_mutex_unlock(pipe->m_mutex);

  processAddNvals(pipe->nvals);
}
void
pstreamCloseParentStream(pstream_t *pstreamptr)
{

  pipe_t *pipe = pstreamptr->pipe;
  pthread_mutex_lock(pipe->m_mutex);
  pipe->EOP = true;
  if (PSTREAM_Debug)
    Message("%s write closed", pstreamptr->m_name.c_str());
  pthread_mutex_unlock(pipe->m_mutex);
  pthread_cond_signal(pipe->tsDef);
  pthread_cond_signal(pipe->tsInq);

  std::unique_lock<std::mutex> locked_mutex(pipe->m_mutex);
  while (pstreamptr->isopen)
    {
      if (PSTREAM_Debug)
        Message("wait of read close");
      pthread_cond_wait(pipe->isclosed, locked_mutex);
    }
  locked_mutex.unlock();
}

void
pstreamClose(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  if (pstreamptr == NULL)
    Error("Internal problem, stream %d not open!", pstreamID);

  pstreamptr->close();

  if(!pstreamptr->ispipe)
  {
    pstream_delete_entry(pstreamptr);
  }
}

void pstream_t::close(){
  if (ispipe)
    {
#if defined(HAVE_LIBPTHREAD)
      pthread_t threadID = pthread_self();

      if (pthread_equal(threadID, rthreadID))
        pstreamCloseChildStream(this);
      else if (pthread_equal(threadID, wthreadID))
        pstreamCloseParentStream(this);
      else
        Error("Internal problem! Close pipe %s", m_name.c_str());

     // processDelStream(pstreamID);
#else
      cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
    }
  else
    {
      if (PSTREAM_Debug)
        Message("%s fileID %d", m_name.c_str(), m_fileID);

      if (mode == 'r')
        {
          processAddNvals(streamNvals(m_fileID));
        }

#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamClose(m_fileID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif

      if (cdoExpMode == CDO_EXP_REMOTE)
        {
          if (mode == 'w')
            {
              extern const char *cdojobfiles;
              FILE *fp = fopen(cdojobfiles, "a");
              fprintf(fp, "%s\n", m_name.c_str());
              fclose(fp);
            }
        }

      if (m_varlist)
        {
          Free(m_varlist);
          m_varlist = NULL;
        }

    }
}

int
pstreamInqVlist(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  return pstreamptr->inqVlist();
}

int
pstream_t::inqVlist()
{
  int vlistID = -1;

#if defined(HAVE_LIBPTHREAD)
  // read from pipe
  if (ispipe)
    {
      vlistID = pipe->pipeInqVlist(m_vlistID);
      if (vlistID == -1)
        cdoAbort("Couldn't read data from input stream %s!", m_name.c_str());
    }
  // read from file through cdi streamInqVlist
  else
#endif
    {
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      vlistID = streamInqVlist(m_fileID);
    if (vlistID == -1){
        cdoAbort("Couldn't read data from input fileID %d!", m_fileID);
    }

#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_read);

      int nsubtypes = vlistNsubtypes(vlistID);
      if (nsubtypes > 1)
        cdoWarning("Subtypes are unsupported, the processing results are possibly wrong!");

      if (cdoDefaultTimeType != CDI_UNDEFID)
        taxisDefType(vlistInqTaxis(vlistID), cdoDefaultTimeType);

      m_vlistID = vlistID;
    }

  if (vlistNumber(vlistID) == CDI_COMP && cdoStreamNumber() == CDI_REAL)
    cdoAbort("Complex fields are not supported by this operator!");

  if (vlistNumber(vlistID) == CDI_REAL && cdoStreamNumber() == CDI_COMP)
    cdoAbort("This operator needs complex fields!");

  processDefVarNum(vlistNvars(vlistID));

  return vlistID;
}

void pstream_t::defVarList(int p_vlistID)
{
  int filetype = m_filetype;

  if (m_vlistID != -1)
    cdoAbort("Internal problem, vlist already defined!");

  if (m_varlist != NULL)
    cdoAbort("Internal problem, varlist already allocated!");

  int nvars = vlistNvars(p_vlistID);
  assert(nvars > 0);

  varlist_t *varlist = (varlist_t *) Malloc(nvars * sizeof(varlist_t));

  for (int varID = 0; varID < nvars; ++varID)
    {
      varlist[varID].gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
      varlist[varID].datatype = vlistInqVarDatatype(p_vlistID, varID);
      varlist[varID].missval = vlistInqVarMissval(p_vlistID, varID);
      varlist[varID].addoffset = vlistInqVarAddoffset(p_vlistID, varID);
      varlist[varID].scalefactor = vlistInqVarScalefactor(p_vlistID, varID);

      varlist[varID].check_datarange = false;

      int laddoffset = IS_NOT_EQUAL(varlist[varID].addoffset, 0);
      int lscalefactor = IS_NOT_EQUAL(varlist[varID].scalefactor, 1);

      int datatype = varlist[varID].datatype;

      if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4
          || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NC5)
        {
          if (datatype == CDI_DATATYPE_UINT8 &&
              (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT16;
              varlist[varID].datatype = datatype;
            }

          if (datatype == CDI_DATATYPE_UINT16 &&
              (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT32;
              varlist[varID].datatype = datatype;
            }

          if (laddoffset || lscalefactor)
            {
              if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
                  || datatype == CDI_DATATYPE_UINT16)
                varlist[varID].check_datarange = true;
            }
          else if (cdoCheckDatarange)
            {
              varlist[varID].check_datarange = true;
            }
        }
    }

  m_varlist = varlist;
  m_vlistID = p_vlistID; /* used for -r/-a */
}

void
pstreamDefVlist(int pstreamID, int vlistID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->defVlist(vlistID);
}

void pstream_t::defVlist(int p_vlistID){
#if defined(HAVE_LIBPTHREAD)
  if (ispipe)
    {
      if (PSTREAM_Debug)
        Message("%s pstreamID %d", m_name.c_str(), self);
      int vlistIDcp = vlistDuplicate(p_vlistID);
      /*    pipeDefVlist(pstreamptr, p_vlistID);*/
      pipe->pipeDefVlist(m_vlistID, vlistIDcp);
    }
  else
#endif
    {
      if (cdoDefaultDataType != CDI_UNDEFID)
        {
          int varID, nvars = vlistNvars(p_vlistID);

          for (varID = 0; varID < nvars; ++varID)
            vlistDefVarDatatype(p_vlistID, varID, cdoDefaultDataType);

          if (cdoDefaultDataType == CDI_DATATYPE_FLT64 || cdoDefaultDataType == CDI_DATATYPE_FLT32)
            {
              for (varID = 0; varID < nvars; varID++)
                {
                  vlistDefVarAddoffset(p_vlistID, varID, 0.0);
                  vlistDefVarScalefactor(p_vlistID, varID, 1.0);
                }
            }
        }

      if (cdoChunkType != CDI_UNDEFID)
        {
          int varID, nvars = vlistNvars(p_vlistID);

          for (varID = 0; varID < nvars; ++varID)
            vlistDefVarChunkType(p_vlistID, varID, cdoChunkType);
        }

      if (CDO_CMOR_Mode)
        {
          cdo_def_tracking_id(p_vlistID, "tracking_id");
          cdo_def_creation_date(p_vlistID);
        }

      if (CDO_Version_Info)
        cdiDefAttTxt(p_vlistID, CDI_GLOBAL, "CDO", (int) strlen(cdoComment()), cdoComment());

#if defined(_OPENMP)
      if (ompNumThreads > 1)
        cdiDefAttInt(p_vlistID, CDI_GLOBAL, "cdo_openmp_thread_number", CDI_DATATYPE_INT32, 1, &ompNumThreads);
#endif
      defVarList(p_vlistID);

      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamDefVlist(m_fileID, p_vlistID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_write);
    }
}

int
pstreamInqRecord(int pstreamID, int *varID, int *levelID)
{
  pstream_t *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pstreamptr->pipe->name.c_str(), pstreamptr->self);
        }
      pstreamptr->pipe->pipeInqRecord(varID, levelID);
    }

  else
#endif
    {
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamInqRecord(pstreamptr->m_fileID, varID, levelID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_read);
    }

  return 0;
}

void
pstreamDefRecord(int pstreamID, int varID, int levelID)
{
  pstream_t *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

  pstreamptr->varID = varID;

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {

      if (PSTREAM_Debug)
        {
          Message("%s pstreamid %d", pstreamptr->m_name.c_str(), pstreamptr->self);
        }
      pstreamptr->pipe->pipeDefRecord(varID, levelID);
    }
  else
#endif
    {
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamDefRecord(pstreamptr->m_fileID, varID, levelID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_write);
    }
}

void
pstreamReadRecord(int pstreamID, double *data, int *nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pstreamptr->pipe->name.c_str(), pstreamptr->self);
        }
      pstreamptr->pipe->pipeReadRecord(pstreamptr->m_vlistID, data, nmiss);
    }
  else
#endif
    {
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamReadRecord(pstreamptr->m_fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_read);
    }
}

void
pstreamReadRecordF(int pstreamID, float *data, int *nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      cdoAbort("pipeReadRecord not implemented for memtype float!");
      // pipeReadRecord(pstreamptr, data, nmiss);
    }
  else
#endif
    {
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamReadRecordF(pstreamptr->m_fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_read);
    }
}

void
pstreamCheckDatarange(pstream_t *pstreamptr, int varID, double *array, int nmiss)
{
  long i;
  long gridsize = pstreamptr->m_varlist[varID].gridsize;
  int datatype = pstreamptr->m_varlist[varID].datatype;
  double missval = pstreamptr->m_varlist[varID].missval;
  double addoffset = pstreamptr->m_varlist[varID].addoffset;
  double scalefactor = pstreamptr->m_varlist[varID].scalefactor;

  long ivals = 0;
  double arrmin = 1.e300;
  double arrmax = -1.e300;
  if (nmiss > 0)
    {
      for (i = 0; i < gridsize; ++i)
        {
          if (!DBL_IS_EQUAL(array[i], missval))
            {
              if (array[i] < arrmin)
                arrmin = array[i];
              if (array[i] > arrmax)
                arrmax = array[i];
              ivals++;
            }
        }
    }
  else
    {
      for (i = 0; i < gridsize; ++i)
        {
          if (array[i] < arrmin)
            arrmin = array[i];
          if (array[i] > arrmax)
            arrmax = array[i];
        }
      ivals = gridsize;
    }

  if (ivals > 0)
    {
      double smin = (arrmin - addoffset) / scalefactor;
      double smax = (arrmax - addoffset) / scalefactor;

      if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
          || datatype == CDI_DATATYPE_UINT16)
        {
          smin = (int) lround(smin);
          smax = (int) lround(smax);
        }

      double vmin = 0, vmax = 0;

      // clang-format off
      if      ( datatype == CDI_DATATYPE_INT8   ) { vmin =        -128.; vmax =        127.; }
      else if ( datatype == CDI_DATATYPE_UINT8  ) { vmin =           0.; vmax =        255.; }
      else if ( datatype == CDI_DATATYPE_INT16  ) { vmin =      -32768.; vmax =      32767.; }
      else if ( datatype == CDI_DATATYPE_UINT16 ) { vmin =           0.; vmax =      65535.; }
      else if ( datatype == CDI_DATATYPE_INT32  ) { vmin = -2147483648.; vmax = 2147483647.; }
      else if ( datatype == CDI_DATATYPE_UINT32 ) { vmin =           0.; vmax = 4294967295.; }
      else if ( datatype == CDI_DATATYPE_FLT32  ) { vmin = -3.40282e+38; vmax = 3.40282e+38; }
      else                                        { vmin =     -1.e+300; vmax =     1.e+300; }
      // clang-format on

      if (smin < vmin || smax > vmax)
        cdoWarning("Some data values (min=%g max=%g) are outside the\n"
                   "    valid range (%g - %g) of the used output precision!\n"
                   "    Use the CDO option%s -b 64 to increase the output precision.",
                   smin,
                   smax,
                   vmin,
                   vmax,
                   (datatype == CDI_DATATYPE_FLT32) ? "" : " -b 32 or");
    }
}

void
pstreamWriteRecord(int pstreamID, double *data, int nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pstreamptr->pipe->name.c_str(), pstreamptr->self);
        }
      pstreamptr->pipe->pipeWriteRecord(data, nmiss);
    }
  else
#endif
    {
      int varID = pstreamptr->varID;
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_write);

      if (pstreamptr->m_varlist)
        if (pstreamptr->m_varlist[varID].check_datarange)
          pstreamCheckDatarange(pstreamptr, varID, data, nmiss);

#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecord(pstreamptr->m_fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif

      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_write);
    }
}

void
pstreamWriteRecordF(int pstreamID, float *data, int nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      cdoAbort("pipeWriteRecord not implemented for memtype float!");
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pstreamptr->pipe->name.c_str(), pstreamptr->self);
        }
      // pipeWriteRecord(pstreamptr, data, nmiss);
    }
  else
#endif
    {
      // int varID = pstreamptr->varID;
      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_write);
/*
if ( pstreamptr->m_varlist )
  if ( pstreamptr->m_varlist[varID].check_datarange )
    pstreamCheckDatarange(pstreamptr, varID, data, nmiss);
*/
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecordF(pstreamptr->m_fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_write);
    }
}

int
pstreamInqTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int nrecs = 0;

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pstreamptr->pipe->name.c_str(), pstreamptr->self);
        }
      nrecs = pstreamptr->pipe->pipeInqTimestep(tsID);
    }
  else
#endif
    {
      if (pstreamptr->mfiles)
        tsID -= pstreamptr->tsID0;

      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      nrecs = streamInqTimestep(pstreamptr->m_fileID, tsID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_read);

      if (nrecs == 0 && pstreamptr->mfiles && (pstreamptr->nfiles < pstreamptr->mfiles))
        {
          int nfile = pstreamptr->nfiles;
          std::string filename; 
          int fileID;
          int vlistIDold, vlistIDnew;

          pstreamptr->tsID0 += tsID;

          vlistIDold = vlistDuplicate(streamInqVlist(pstreamptr->m_fileID));
          streamClose(pstreamptr->m_fileID);

          filename = pstreamptr->m_mfnames[nfile];
          pstreamptr->nfiles++;

#if defined(HAVE_LIBPTHREAD)
          if (cdoLockIO)
            pthread_mutex_lock(&streamMutex);
          else
            pthread_mutex_lock(&streamOpenReadMutex);
#endif
          if (cdoVerbose)
            cdoPrint("Continuation file: %s", filename.c_str());

          if (processNums() == 1 && ompNumThreads == 1)
            timer_start(timer_read);
          fileID = streamOpenRead(filename.c_str());
          vlistIDnew = streamInqVlist(fileID);
          if (processNums() == 1 && ompNumThreads == 1)
            timer_stop(timer_read);

          vlistCompare(vlistIDold, vlistIDnew, CMP_HRD);
          vlistDestroy(vlistIDold);
#if defined(HAVE_LIBPTHREAD)
          if (cdoLockIO)
            pthread_mutex_unlock(&streamMutex);
          else
            pthread_mutex_unlock(&streamOpenReadMutex);
#endif
          if (fileID < 0)
            cdiOpenError(fileID, "Open failed on >%s<", filename.c_str());

          pstreamptr->m_name = filename;
          pstreamptr->m_fileID = fileID;

          if (processNums() == 1 && ompNumThreads == 1)
            timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
          if (cdoLockIO)
            pthread_mutex_lock(&streamMutex);
#endif
          nrecs = streamInqTimestep(pstreamptr->m_fileID, 0);
#if defined(HAVE_LIBPTHREAD)
          if (cdoLockIO)
            pthread_mutex_unlock(&streamMutex);
#endif
          if (processNums() == 1 && ompNumThreads == 1)
            timer_stop(timer_read);
        }

      if (tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID)
        taxisDefType(vlistInqTaxis(pstreamptr->m_vlistID), cdoDefaultTimeType);
    }

  if (nrecs && tsID != pstreamptr->tsID)
    {
      processDefTimesteps(pstreamID);
      pstreamptr->tsID = tsID;
    }

  return nrecs;
}

void
pstreamDefTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->defTimestep(tsID);
}

void
pstream_t::defTimestep(int p_tsID)
{
#if defined(HAVE_LIBPTHREAD)
  if (ispipe)
    {
      if (PSTREAM_Debug)
        {
          Message("%s pstreamID %d", pipe->name.c_str(), self);
        }
      pipe->pipeDefTimestep(m_vlistID, p_tsID);
    }
  else
#endif
    {
      if (p_tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID)
        {
          int taxisID, vlistID;
          vlistID = m_vlistID;
          taxisID = vlistInqTaxis(vlistID);
          taxisDefType(taxisID, cdoDefaultTimeType);
        }

      if (processNums() == 1 && ompNumThreads == 1)
        timer_start(timer_write);
/* don't use sync -> very slow on GPFS */
//  if ( p_tsID > 0 ) streamSync(fileID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_lock(&streamMutex);
#endif
      streamDefTimestep(m_fileID, p_tsID);
#if defined(HAVE_LIBPTHREAD)
      if (cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
#endif
      if (processNums() == 1 && ompNumThreads == 1)
        timer_stop(timer_write);
    }
}

void
pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc)
{
  if (PSTREAM_Debug)
    Message("pstreamIDdest = %d  pstreamIDsrc = %d", pstreamIDdest, pstreamIDsrc);

  pstream_t *pstreamptr_dest = pstream_to_pointer(pstreamIDdest);
  pstream_t *pstreamptr_src = pstream_to_pointer(pstreamIDsrc);

  if (pstreamptr_dest->ispipe || pstreamptr_src->ispipe)
    cdoAbort("This operator can't be combined with other operators!");

#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_lock(&streamMutex);
#endif
  streamCopyRecord(pstreamptr_dest->m_fileID, pstreamptr_src->m_fileID);
#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
#endif
}

void
pstreamDebug(int debug)
{
  PSTREAM_Debug = debug;
}

void
cdoInitialize(void *argument)
{
#if defined(_OPENMP)
  omp_set_num_threads(ompNumThreads); /* Have to be called for every module (pthread)! */
#endif

  processCreate();

#if defined(HAVE_LIBPTHREAD)
  if (PSTREAM_Debug)
    Message("process %d  thread %ld", processSelf().m_ID, pthread_self());
#endif

  processDefArgument(argument);
}

void
pstreamCloseAll()
{
  for (auto pstream_iter : _pstream_map)
    {
      if ( pstream_iter.second.m_fileID != CDI_UNDEFID )
        {
          if (PSTREAM_Debug)
            Message("Close file %s id %d", pstream_iter.second.m_name.c_str(), pstream_iter.second.m_fileID);
          streamClose(pstream_iter.second.m_fileID);
        }
    }
}
/*
void
pstreamCloseAll(void)
{
  if (_pstreamList == NULL)
    return;

  for (int i = 0; i < _pstream_max; i++)
    {
      pstream_t *pstreamptr = _pstreamList[i].ptr;
      if (pstreamptr && pstreamptr->isopen)
        {
          if (!pstreamptr->ispipe && pstreamptr->m_fileID != CDI_UNDEFID)
            {
              if (PSTREAM_Debug)
                Message("Close file %s id %d", pstreamptr->m_name.c_str(), pstreamptr->m_fileID);
              streamClose(pstreamptr->m_fileID);
            }
        }
    }
}
*/


int
pstreamInqFiletype(int pstreamID)
{
  return pstream_to_pointer(pstreamID)->inqFileType();
}
int
pstream_t::inqFileType()
{
  int filetype;

#if defined(HAVE_LIBPTHREAD)
  if (ispipe)
    {
      filetype = m_filetype;
    }
  else
#endif
    {
      filetype = streamInqFiletype(m_fileID);
    }
  return filetype;
}

int
pstreamInqByteorder(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int byteorder;

#if defined(HAVE_LIBPTHREAD)
  if (pstreamptr->ispipe)
    byteorder = pstreamptr->m_filetype;
  else
#endif
    byteorder = streamInqByteorder(pstreamptr->m_fileID);

  return byteorder;
}

void
pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  streamInqGRIBinfo(pstreamptr->m_fileID, intnum, fltnum, bignum);
}

int
pstreamFileID(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  return pstreamptr->m_fileID;
}

void
cdoVlistCopyFlag(int vlistID2, int vlistID1)
{
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&streamMutex);
#endif

  vlistCopyFlag(vlistID2, vlistID1);

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&streamMutex);
#endif
}

void
openLock(void)
{
#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenReadMutex);
#endif
}

void
openUnlock(void)
{
#if defined(HAVE_LIBPTHREAD)
  if (cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenReadMutex);
#endif
}

//TODO remove when processes create the new threads
const int &getPthreadScope()
{
    return pthreadScope;
}

