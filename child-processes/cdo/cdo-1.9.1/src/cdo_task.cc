#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>
#endif

#include "cdo_task.h"

bool CDO_task = false;

enum cdo_pt_state {
  SETUP,
  IDLE,
  JOB,
  DIE
};

typedef struct cdo_task_info {
  void *(*routine)(void *);
  void *arg;
  void *result;
  enum cdo_pt_state state;
#if defined(HAVE_LIBPTHREAD)
  pthread_t thread;
  pthread_cond_t work_cond;
  pthread_mutex_t work_mtx;
  pthread_cond_t boss_cond;
  pthread_mutex_t boss_mtx;
#endif
} cdo_task_t;


#if defined(HAVE_LIBPTHREAD)
static
void *cdo_task(void *task)
{
  cdo_task_t *task_info = (cdo_task_t *) task;

  // cond_wait mutex must be locked before we can wait
  pthread_mutex_lock(&(task_info->work_mtx));

  //printf("<worker> start\n");

  // ensure boss is waiting
  pthread_mutex_lock(&(task_info->boss_mtx));

  // signal to boss that setup is complete
  task_info->state = IDLE;

  // wake-up signal
  pthread_cond_signal(&(task_info->boss_cond));
  pthread_mutex_unlock(&(task_info->boss_mtx));

  while ( 1 )
    {
      pthread_cond_wait(&(task_info->work_cond), &(task_info->work_mtx));

      if ( DIE == task_info->state ) break; // kill thread

      if ( IDLE == task_info->state ) continue; // accidental wake-up
      
      // do blocking task
      //printf("<worker> JOB start\n");
      task_info->result = task_info->routine(task_info->arg);
      //printf("<worker> JOB end\n");

      // ensure boss is waiting
      pthread_mutex_lock(&(task_info->boss_mtx));

      // indicate that job is done
      task_info->state = IDLE;

      // wake-up signal
      pthread_cond_signal(&(task_info->boss_cond));
      pthread_mutex_unlock(&(task_info->boss_mtx));
    }

  pthread_mutex_unlock(&(task_info->work_mtx));
  pthread_exit(NULL);

  return NULL;
}
#endif

void cdo_task_start(void *task, void *(*task_routine)(void *), void *task_arg)
{
  if ( !task ) return;

  cdo_task_t *task_info = (cdo_task_t *) task;

  // ensure worker is waiting
#if defined(HAVE_LIBPTHREAD)
  if ( CDO_task ) pthread_mutex_lock(&(task_info->work_mtx));
#endif

  // set job information & state
  task_info->routine = task_routine;
  task_info->arg = task_arg;
  task_info->state = JOB;

  bool run_task = !CDO_task;
#if !defined(HAVE_LIBPTHREAD)
  run_task = true;
#endif
  if ( run_task ) task_info->result = task_info->routine(task_info->arg);

  // wake-up signal
#if defined(HAVE_LIBPTHREAD)
  if ( CDO_task )
    {
      pthread_cond_signal(&(task_info->work_cond));
      pthread_mutex_unlock(&(task_info->work_mtx));
    }
#endif
}


void *cdo_task_wait(void *task)
{
  if ( !task ) return NULL;

  cdo_task_t *task_info = (cdo_task_t *) task;

#if defined(HAVE_LIBPTHREAD)
  if ( CDO_task )
    {
      while (1)
        {
          if ( IDLE == task_info->state ) break;

          pthread_cond_wait(&(task_info->boss_cond), &(task_info->boss_mtx));

          // if ( IDLE == task_info->state ) break;
        }
    }
#endif

  return task_info->result;
}


void *cdo_task_new()
{
  cdo_task_t *task_info = NULL;

  task_info = (cdo_task_t *) malloc(sizeof(cdo_task_t));
  task_info->routine = NULL;
  task_info->arg = NULL;
  task_info->result = NULL;
  task_info->state = SETUP;

#if defined(HAVE_LIBPTHREAD)
  if ( CDO_task )
    {
      pthread_attr_t attr;
      size_t stacksize;
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      pthread_attr_getstacksize(&attr, &stacksize);
      if ( stacksize < 2097152 )
        {
          stacksize = 2097152;
          pthread_attr_setstacksize(&attr, stacksize);
        }

      pthread_cond_init(&(task_info->work_cond), NULL);
      pthread_mutex_init(&(task_info->work_mtx), NULL);
      pthread_cond_init(&(task_info->boss_cond), NULL);
      pthread_mutex_init(&(task_info->boss_mtx), NULL);

      pthread_mutex_lock(&(task_info->boss_mtx));

      pthread_create(&(task_info->thread), &attr, cdo_task, (void *)task_info);

      cdo_task_wait(task_info);
    }
#endif

  return (void *)task_info;
}


void cdo_task_delete(void *task)
{
  cdo_task_t *task_info = (cdo_task_t *) task;

#if defined(HAVE_LIBPTHREAD)
  if ( CDO_task )
    {
      // ensure the worker is waiting
      pthread_mutex_lock(&(task_info->work_mtx));

      //printf("cdo_task_delete: send DIE to <worker>\n");
      task_info->state = DIE;

      // wake-up signal
      pthread_cond_signal(&(task_info->work_cond));
      pthread_mutex_unlock(&(task_info->work_mtx));

      // wait for thread to exit
      pthread_join(task_info->thread, NULL);

      pthread_mutex_destroy(&(task_info->work_mtx));
      pthread_cond_destroy(&(task_info->work_cond));

      pthread_mutex_unlock(&(task_info->boss_mtx));
      pthread_mutex_destroy(&(task_info->boss_mtx));
      pthread_cond_destroy(&(task_info->boss_cond));
    }
#endif

  if ( task_info ) free(task_info);
}

#ifdef TEST_CDO_TASK
// gcc -DTEST_CDO_TASK -DHAVE_LIBPTHREAD cdo_task.c

void *mytask(void *arg)
{
  printf("run mytask\n");
}


int main(int argc, char **argv)
{
  CDO_task = true;

  void *task = cdo_task_new();

  printf("Init done\n");
  void *myarg = NULL;
  void *myresult;

  cdo_task_start(task, mytask, myarg);

  myresult = cdo_task_wait(task);
  
  cdo_task_start(task, mytask, myarg);

  myresult = cdo_task_wait(task);

  cdo_task_delete(task);

  return 0;
}
#endif
