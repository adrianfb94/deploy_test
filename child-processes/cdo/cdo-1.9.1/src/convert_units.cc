#include <errno.h>
#include "convert_units.h"

#if defined(HAVE_UDUNITS2)

static void udunitsInitialize(void);
static int udunitsInit = 0;

#if defined(HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  udunitsInitThread = PTHREAD_ONCE_INIT;
static pthread_mutex_t udunitsMutex;

#  define UDUNITS_LOCK()         pthread_mutex_lock(&udunitsMutex)
#  define UDUNITS_UNLOCK()       pthread_mutex_unlock(&udunitsMutex)
#  define UDUNITS_INIT()         pthread_once(&udunitsInitThread, udunitsInitialize)

#else

#  define UDUNITS_LOCK()
#  define UDUNITS_UNLOCK()
#  define UDUNITS_INIT()         if ( !udunitsInit ) udunitsInitialize();

#endif

static ut_system *ut_read = NULL;

static
void udunitsInitialize(void)
{
#if defined(HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&udunitsMutex, NULL);
#endif

  udunitsInit = 1;
}

static
void *get_converter(char *src_unit_str, char *tgt_unit_str, int *rstatus)
{
  ut_unit *src_unit, *tgt_unit;
  cv_converter *ut_units_converter = NULL;
  int status;

  *rstatus = -1;

  if ( ut_read == NULL )
    {
      ut_set_error_message_handler(ut_ignore);

      errno = 0;
      ut_read = ut_read_xml(NULL);
      status = ut_get_status();
      if ( status == UT_PARSE )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Couldn't parse unit database!");
	}
      if ( status == UT_OPEN_ENV || status == UT_OPEN_DEFAULT || status == UT_OS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: %s", strerror(errno));
	}
      errno = 0;
      if ( status != UT_SUCCESS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Error reading units system!");
	  return NULL;
	}
    }

  ut_trim(src_unit_str, UT_ASCII);
  src_unit = ut_parse(ut_read, src_unit_str, UT_ASCII);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error parsing units: [%s]", src_unit_str);
      return NULL;
    }

  ut_trim(tgt_unit_str, UT_ASCII);
  tgt_unit = ut_parse(ut_read, tgt_unit_str, UT_ASCII);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error parsing units: [%s]", tgt_unit_str);
      return NULL;
    }

  status = ut_compare(src_unit, tgt_unit);
  if ( status == 0 ) *rstatus = -2;

  if ( *rstatus == -1 )
    {
      status = ut_are_convertible(src_unit, tgt_unit);
      if ( status == 0 ) *rstatus = -3;
    }

  if ( *rstatus == -1 )
    {
      ut_units_converter = ut_get_converter(src_unit, tgt_unit);
      if ( ut_units_converter == NULL || ut_get_status() != UT_SUCCESS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Error getting converter from [%s] to [%s]", src_unit_str, tgt_unit_str);
	}
      else
	*rstatus = 0;
    }

  ut_free(src_unit);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error freeing units [%s]", src_unit_str);
      return NULL;
    }
     
  ut_free(tgt_unit);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error freeing units [%s]", tgt_unit_str);
      return NULL;
    }

  return (void *) ut_units_converter;
}


void cdoConvertFree(void *ut_converter)
{
  UDUNITS_LOCK();
  cv_free((cv_converter*)ut_converter);
  UDUNITS_UNLOCK();
}


void cdoConvertDestroy()
{
  UDUNITS_LOCK();
  if ( ut_read )
    { 
      ut_free_system(ut_read);
      ut_read = NULL;
    }
  UDUNITS_UNLOCK();
}
#endif

void cdoConvertUnits(void **ut_converter, bool *changeunits, char *units, char *units_old, const char *name)
{
  if ( *changeunits )
    {
#if defined(HAVE_UDUNITS2)
      int status;
      UDUNITS_INIT();
      UDUNITS_LOCK();
      *ut_converter = get_converter(units_old, units, &status);
      UDUNITS_UNLOCK();
      if ( *ut_converter == NULL )
	{
	  if ( status == -2 )
	    {
	      if ( cdoVerbose )
		cdoPrint("%s - not converted from  [%s] to [%s], units are equal!", name, units_old, units);
	    }
	  else if ( status == -3 )
	    {
	      cdoWarning("%s - converting units from [%s] to [%s] failed, not convertible!", name, units_old, units);
	    }
	  else
	    cdoWarning("%s - converting units from [%s] to [%s] failed!", name, units_old, units);
	  *changeunits = false;
	}
      else
	{
	  // if ( cdoVerbose )
	    {
	      char buf[64];
	      cv_get_expression((const cv_converter*)*ut_converter, buf, sizeof(buf), name);
	      cdoPrint("%s - convert units from [%s] to [%s] (expression: %s).", name, units_old, units, buf);
	    }
	}
#else
      static bool lwarn_udunits = true;
      if ( lwarn_udunits )
	{
	  cdoWarning("%s - converting units from [%s] to [%s] failed, UDUNITS2 support not compiled in!", name, units_old, units);
	  *changeunits = false;
	  lwarn_udunits = false;
	}
#endif
    }
}
