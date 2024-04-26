#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBGRIB_API
#include <limits.h>
#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "gribapi_utilities.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"
#include "calendar.h"
#include "subtype.h"


#include "cgribex.h"      /* gribGetSize, gribRead, gribGetZip, GRIB1_LTYPE_99 */
#include "gribapi.h"

#include <grib_api.h>

extern int cdiInventoryMode;

static const var_tile_t dummy_tiles = { 0, -1, -1, -1, -1, -1 };

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
#ifdef HIRLAM_EXTENSIONS
    // NOTE: tsteptype MUST be part of attributes used to compare variables!
    // Modern NWP models (HARMONIE, HIRLAM) use timeRangeIndicator to specify
    // if the field is instantanous or accumulated.
    // Both types are typically in the same GRIB-file.
    // (181; 105, 0, timeRangeIndicator=0) .. instantanous rain
    // (181; 105, 0, timeRangeIndicator=4) .. accumulated rain  .. both can be in the same grib file
#endif // HIRLAM_EXTENSIONS
  char name[32];

  var_tile_t tiles;

} compvar2_t;


static
int gribapiGetZaxisType(long editionNumber, int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  if ( editionNumber <= 1 )
    {
      zaxistype = grib1ltypeToZaxisType(grib_ltype);
    }
  else
    {
      zaxistype = grib2ltypeToZaxisType(grib_ltype);
    }

  return zaxistype;
}

static
int getTimeunits(long unitsOfTime)
{
  int timeunits = -1;

  switch (unitsOfTime)
    {
    case 13:  timeunits = TUNIT_SECOND;  break;
    case  0:  timeunits = TUNIT_MINUTE;  break;
    case  1:  timeunits = TUNIT_HOUR;    break;
    case 10:  timeunits = TUNIT_3HOURS;  break;
    case 11:  timeunits = TUNIT_6HOURS;  break;
    case 12:  timeunits = TUNIT_12HOURS; break;
    case  2:  timeunits = TUNIT_DAY;     break;
    default:  timeunits = TUNIT_HOUR;    break;
    }

  return timeunits;
}

static
double timeunit_factor(int tu1, int tu2)
{
  double factor = 1;

  if ( tu2 == TUNIT_HOUR )
    {
      switch (tu1)
        {
        case TUNIT_SECOND:  factor = 3600;   break;
        case TUNIT_MINUTE:  factor = 60;     break;
        case TUNIT_HOUR:    factor = 1;      break;
        case TUNIT_3HOURS:  factor = 1./3;   break;
        case TUNIT_6HOURS:  factor = 1./6;   break;
        case TUNIT_12HOURS: factor = 1./12;  break;
        case TUNIT_DAY:     factor = 1./24;  break;
        }
    }

  return factor;
}

static
int gribapiGetTimeUnits(grib_handle *gh)
{
  int timeunits = -1;
  long unitsOfTime = -1;

  grib_get_long(gh, "indicatorOfUnitOfTimeRange", &unitsOfTime);

  GRIB_CHECK(my_grib_set_long(gh, "stepUnits", unitsOfTime), 0);

  timeunits = getTimeunits(unitsOfTime);

  return timeunits;
}

static
void gribapiGetSteps(grib_handle *gh, int timeunits, int *startStep, int *endStep)
{
  int timeunits2 = timeunits;
  long unitsOfTime;
  int status = grib_get_long(gh, "stepUnits", &unitsOfTime);
  if ( status == 0 ) timeunits2 = getTimeunits(unitsOfTime);
  //timeunits2 = gribapiGetTimeUnits(gh);

  long lpar;
  status = grib_get_long(gh, "forecastTime", &lpar);
  if ( status == 0 ) *startStep = (int) lpar;
  else
    {
      status = grib_get_long(gh, "startStep", &lpar);
      if ( status == 0 )
        *startStep = (int) (((double)lpar * timeunit_factor(timeunits, timeunits2)) + 0.5);
    }

  *endStep = *startStep;
  status = grib_get_long(gh, "endStep", &lpar);
  if ( status == 0 )
    *endStep = (int) (((double)lpar * timeunit_factor(timeunits, timeunits2)) + 0.5);
  // printf("%d %d %d %d %d %g\n", *startStep, *endStep, lpar, timeunits, timeunits2, timeunit_factor(timeunits, timeunits2));
}

static
void gribapiGetDataDateTime(grib_handle *gh, int *datadate, int *datatime)
{
  long lpar;

  GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
  *datadate = (int) lpar;
  GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);  //FIXME: This looses the seconds in GRIB2 files.
  *datatime = (int) lpar*100;
}

static
void gribapiSetDataDateTime(grib_handle *gh, int datadate, int datatime)
{
  GRIB_CHECK(my_grib_set_long(gh, "dataDate", datadate), 0);
  GRIB_CHECK(my_grib_set_long(gh, "dataTime", datatime/100), 0);
}

static
int gribapiGetValidityDateTime(grib_handle *gh, int *vdate, int *vtime)
{
  int rdate, rtime;
  int timeUnits, startStep = 0, endStep;
  int tstepRange = 0;
  int range;
  long sigofrtime = 3;

  if ( gribEditionNumber(gh) > 1 )
    {
      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);
    }
  else
    {
      GRIB_CHECK(grib_get_long(gh, "timeRangeIndicator", &sigofrtime), 0);
    }

  if ( sigofrtime == 3 )        //XXX: This looks like a bug to me, because timeRangeIndicator == 3 does not seem to have the same meaning as significanceOfReferenceTime == 3. I would recommend replacing this condition with `if(!gribapiTimeIsFC())`.
    {
      gribapiGetDataDateTime(gh, vdate, vtime);
    }
  else
    {
      gribapiGetDataDateTime(gh, &rdate, &rtime);

      timeUnits = gribapiGetTimeUnits(gh);
      gribapiGetSteps(gh, timeUnits, &startStep, &endStep);

      range = endStep - startStep;

      if ( range > 0 )
	{
	  if ( startStep == 0 ) tstepRange = -1;
	  else                  tstepRange =  1;
	}

      {
	static bool lprint = true;
	extern int grib_calendar;
	int ryear, rmonth, rday, rhour, rminute, rsecond;
	int julday, secofday;
	int64_t time_period = endStep;
        int64_t addsec;

	cdiDecodeDate(rdate, &ryear, &rmonth, &rday);
	cdiDecodeTime(rtime, &rhour, &rminute, &rsecond);

        if ( rday > 0 )
          {
            encode_caldaysec(grib_calendar, ryear, rmonth, rday, rhour, rminute, rsecond, &julday, &secofday);

            addsec = 0;
            switch ( timeUnits )
              {
              case TUNIT_SECOND:  addsec =         time_period; break;
              case TUNIT_MINUTE:  addsec =    60 * time_period; break;
              case TUNIT_HOUR:    addsec =  3600 * time_period; break;
              case TUNIT_3HOURS:  addsec = 10800 * time_period; break;
              case TUNIT_6HOURS:  addsec = 21600 * time_period; break;
              case TUNIT_12HOURS: addsec = 43200 * time_period; break;
              case TUNIT_DAY:     addsec = 86400 * time_period; break;
              default:
                if ( lprint )
                  {
                    Warning("Time unit %d unsupported", timeUnits);
                    lprint = false;
                  }
                break;
              }

            julday_add_seconds(addsec, &julday, &secofday);

            decode_caldaysec(grib_calendar, julday, secofday, &ryear, &rmonth, &rday, &rhour, &rminute, &rsecond);
          }

	*vdate = cdiEncodeDate(ryear, rmonth, rday);
	*vtime = cdiEncodeTime(rhour, rminute, rsecond);
      }
    }

  return tstepRange;
}

static
void grib1GetLevel(grib_handle *gh, int *leveltype, int *lbounds, int *level1, int *level2)
{
  *leveltype = 0;
  *lbounds   = 0;
  *level1    = 0;
  *level2    = 0;

  long lpar;
  if ( !grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar) )       //1 byte
    {
      *leveltype = (int) lpar;

      switch (*leveltype)
	{
	case GRIB1_LTYPE_SIGMA_LAYER:
	case GRIB1_LTYPE_HYBRID_LAYER:
	case GRIB1_LTYPE_LANDDEPTH_LAYER:
	  { *lbounds = 1; break; }
	}

      if ( *lbounds )
	{
	  GRIB_CHECK(grib_get_long(gh, "topLevel", &lpar), 0);  //1 byte
	  *level1 = (int)lpar;
	  GRIB_CHECK(grib_get_long(gh, "bottomLevel", &lpar), 0);       //1 byte
	  *level2 = (int)lpar;
	}
      else
	{
          double dlevel;
	  GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0); //2 byte
	  if ( *leveltype == GRIB1_LTYPE_ISOBARIC ) dlevel *= 100;
	  if ( dlevel < -2.e9 || dlevel > 2.e9 ) dlevel = 0;
	  if ( *leveltype == GRIB1_LTYPE_99 || *leveltype == GRIB1_LTYPE_ISOBARIC_PA ) *leveltype = GRIB1_LTYPE_ISOBARIC;

	  *level1 = (int) dlevel;
	  *level2 = 0;
	}
    }
}

static
double grib2ScaleFactor(long factor)
{
  switch(factor)
    {
      case GRIB_MISSING_LONG: return 1;
      case 0: return 1;
      case 1: return 0.1;
      case 2: return 0.01;
      case 3: return 0.001;
      case 4: return 0.0001;
      case 5: return 0.00001;
      case 6: return 0.000001;
      case 7: return 0.0000001;
      case 8: return 0.00000001;
      case 9: return 0.000000001;
      default: return 0;
    }
}

static
int calcLevel(int level_sf, long factor, long level)
{
  double result = 0;
  if(level != GRIB_MISSING_LONG) result = (double)level*grib2ScaleFactor(factor);
  if(level_sf) result *= level_sf;
  return (int)result;
}

static
void grib2GetLevel(grib_handle *gh, int *leveltype1, int *leveltype2, int *lbounds, int *level1,
                   int *level2, int *level_sf, int *level_unit)
{
  int status;
  long lpar;
  long factor;

  *leveltype1 = 0;
  *leveltype2 = -1;
  *lbounds    = 0;
  *level1     = 0;
  *level2     = 0;
  *level_sf   = 0;
  *level_unit = 0;

  status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar); //1 byte
  if ( status == 0 )
    {
      long llevel;

      *leveltype1 = (int) lpar;

      status = grib_get_long(gh, "typeOfSecondFixedSurface", &lpar); //1 byte
      /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
      if ( status == 0 ) *leveltype2 = (int)lpar;

      if ( *leveltype1 != 255 && *leveltype2 != 255 && *leveltype2 > 0 ) *lbounds = 1;
      switch(*leveltype1)
        {
          case GRIB2_LTYPE_REFERENCE:
            if(*leveltype2 == 1) *lbounds = 0;
            break;

          case GRIB2_LTYPE_LANDDEPTH:
            *level_sf = 1000;
            *level_unit = CDI_UNIT_M;
            break;

          case GRIB2_LTYPE_ISOBARIC:
            *level_sf = 1000;
            *level_unit = CDI_UNIT_PA;
            break;

          case GRIB2_LTYPE_SIGMA:
            *level_sf = 1000;
            *level_unit = 0;
            break;
        }

      GRIB_CHECK(grib_get_long(gh, "scaleFactorOfFirstFixedSurface", &factor), 0);      //1 byte
      GRIB_CHECK(grib_get_long(gh, "scaledValueOfFirstFixedSurface", &llevel), 0);      //4 byte
      *level1 = calcLevel(*level_sf, factor, llevel);

      if ( *lbounds )
        {
          GRIB_CHECK(grib_get_long(gh, "scaleFactorOfSecondFixedSurface", &factor), 0); //1 byte
          GRIB_CHECK(grib_get_long(gh, "scaledValueOfSecondFixedSurface", &llevel), 0); //4 byte
          *level2 = calcLevel(*level_sf, factor, llevel);
        }
    }
}

static
void gribGetLevel(grib_handle *gh, int* leveltype1, int* leveltype2, int* lbounds, int* level1, int* level2, int* level_sf, int* level_unit, var_tile_t* tiles)
{
  if ( gribEditionNumber(gh) <= 1 )
    {
      grib1GetLevel(gh, leveltype1, lbounds, level1, level2);
      *leveltype2 = -1;
      *level_sf = 0;
      *level_unit = 0;
    }
  else
    {
      grib2GetLevel(gh, leveltype1, leveltype2, lbounds, level1, level2, level_sf, level_unit);

      /* read in tiles attributes (if there are any) */
      tiles->tileindex = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEINDEX], -1);
      tiles->totalno_of_tileattr_pairs = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS], -1);
      tiles->tileClassification = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILE_CLASSIFICATION], -1);
      tiles->numberOfTiles = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_TILES], -1);
      tiles->numberOfAttributes = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_ATTR], -1);
      tiles->attribute = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEATTRIBUTE], -1);
    }
}

static
void gribapiGetString(grib_handle *gh, const char *key, char *string, size_t length)
{
  string[0] = 0;

  int ret = grib_get_string(gh, key, string, &length);
  if (ret != 0)
    {
      fprintf(stderr, "grib_get_string(gh, \"%s\", ...) failed!\n", key);
      GRIB_CHECK(ret, 0);
    }
  if      ( length == 8 && memcmp(string, "unknown", length) == 0 ) string[0] = 0;
  else if ( length == 2 && memcmp(string, "~", length)       == 0 ) string[0] = 0;
}

static
void gribapiAddRecord(stream_t * streamptr, int param, grib_handle *gh,
                      size_t recsize, off_t position, int datatype, int comptype, const char *varname,
                      int leveltype1, int leveltype2, int lbounds, int level1, int level2, int level_sf, int level_unit,
                      const var_tile_t *tiles, int lread_additional_keys)
{
  int levelID = 0;
  char stdname[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  long ens_index = 0, ens_count = 0, ens_forecast_type = 0;

  int vlistID = streamptr->vlistID;
  int tsID    = streamptr->curTsID;
  int recID   = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  int tsteptype = gribapiGetTsteptype(gh);
  // numavg  = ISEC1_AvgNum;
  int numavg = 0;

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype1);

  record->size      = recsize;
  record->position  = position;
  record->param     = param;
  record->ilevel    = level1;
  record->ilevel2   = level2;
  record->ltype     = leveltype1;
  record->tsteptype = (short)tsteptype;
  record->tiles = tiles ? *tiles : dummy_tiles;

  //FIXME: This may leave the variable name unterminated (which is the behavior that I found in the code).
  //       I don't know precisely how this field is used, so I did not change this behavior to avoid regressions,
  //       but I think that it would be better to at least add a line
  //
  //           record->varname[sizeof(record->varname) - 1] = 0;`
  //
  //       after the `strncpy()` call.
  //
  //       I would consider using strdup() (that requires POSIX-2008 compliance, though), or a similar homebrew approach.
  //       I. e. kick the fixed size array and allocate enough space, whatever that may be.
  strncpy(record->varname, varname, sizeof(record->varname));

  grid_t *grid = (grid_t *)Malloc(sizeof(*grid));
  gribapiGetGrid(gh, grid);

  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 0);
  int gridID = gridAdded.Id;
  if ( !gridAdded.isNew ) Free(grid);
  else if ( grid->projtype == CDI_PROJ_RLL )
    {
      double xpole = 0, ypole = 0, angle = 0;
      grib_get_double(gh, "latitudeOfSouthernPoleInDegrees",  &ypole);
      grib_get_double(gh, "longitudeOfSouthernPoleInDegrees", &xpole);
      grib_get_double(gh, "angleOfRotation", &angle);
      xpole -= 180;
      if ( fabs(ypole) > 0 ) ypole = -ypole; // change from south to north pole
      if ( fabs(angle) > 0 ) angle = -angle;

      gridDefParamRLL(gridID, xpole, ypole, angle);
    }
  else if ( grid->projtype == CDI_PROJ_LCC )
    {
      double a = 6367470., rf = 0;
      long earthIsOblate;
      grib_get_long(gh, "earthIsOblate", &earthIsOblate);
      if ( earthIsOblate ) { a = 6378160.; rf = 297.0; }
      double lon_0, lat_1, lat_2, xval_0, yval_0;
      long projflag = 0;
      grib_get_double(gh, "longitudeOfFirstGridPointInDegrees", &xval_0);
      grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &yval_0);
      grib_get_double(gh, "LoVInDegrees", &lon_0);
      grib_get_double(gh, "Latin1InDegrees", &lat_1);
      grib_get_double(gh, "Latin2InDegrees", &lat_2);
      grib_get_long(gh, "projectionCentreFlag", &projflag);
      bool lsouth = gribbyte_get_bit((int)projflag, 1);
      if ( lsouth ) { lat_1 = -lat_1; lat_2 = -lat_2; }

      double lat_0 = lat_2;
      double x_0 = grid_missval;
      double y_0 = grid_missval;

      if ( proj_lonlat_to_lcc_func )
        {
          x_0 = xval_0; y_0 = yval_0;
          proj_lonlat_to_lcc_func(grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, (size_t)1, &x_0, &y_0);
          if ( IS_NOT_EQUAL(x_0, grid_missval) && IS_NOT_EQUAL(y_0, grid_missval) )
            { x_0 = -x_0; y_0 = -y_0; }
        }
      gridDefParamLCC(gridID, grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0);
    }

  int zaxistype = gribapiGetZaxisType(gribEditionNumber(gh), leveltype1);

  switch (zaxistype)
    {
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        long lpar;
        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        /* FIXME: assert(lpar >= 0) */
        size_t vctsize = (size_t)lpar;
        if ( vctsize > 0 )
          {
            double *vctptr = (double *) Malloc(vctsize*sizeof(double));
            size_t dummy = vctsize;
            GRIB_CHECK(grib_get_double_array(gh, "pv", vctptr, &dummy), 0);
            varDefVCT(vctsize, vctptr);
            Free(vctptr);
          }
        break;
      }
    case ZAXIS_REFERENCE:
      {
        unsigned char uuid[CDI_UUID_SIZE];
        long lpar;
        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        if ( lpar != 6 ) fprintf(stderr, "Warning ...\n");
        GRIB_CHECK(grib_get_long(gh, "nlev", &lpar), 0);
        int nhlev = (int)lpar;
        GRIB_CHECK(grib_get_long(gh, "numberOfVGridUsed", &lpar), 0);
        int nvgrid = (int)lpar;
        size_t len = (size_t)CDI_UUID_SIZE;
        memset(uuid, 0, CDI_UUID_SIZE);
        GRIB_CHECK(grib_get_bytes(gh, "uuidOfVGrid", uuid, &len), 0);
        varDefZAxisReference(nhlev, nvgrid, uuid);
        break;
      }
    }

  // if ( datatype > 32 ) datatype = CDI_DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = CDI_DATATYPE_PACK;

  stdname[0] = 0;
  longname[0] = 0;
  units[0] = 0;

  if ( varname[0] != 0 )
    {
      size_t vlen = CDI_MAX_NAME;
      gribapiGetString(gh, "name", longname, vlen);
      vlen = CDI_MAX_NAME;
      gribapiGetString(gh, "units", units, vlen);
      vlen = CDI_MAX_NAME;
      int status = grib_get_string(gh, "cfName", stdname, &vlen);
      if ( status != 0 || vlen <= 1 || strncmp(stdname, "unknown", 7) == 0 )
        stdname[0] = 0;
    }
  // fprintf(stderr, "param %d name %s %s %s\n", param, name, longname, units);

  /* add the previously read record data to the (intermediate) list of records */
  int tile_index = 0, varID;
  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, level_sf, level_unit,
	       datatype, &varID, &levelID, tsteptype, numavg, leveltype1, leveltype2,
	       varname, stdname, longname, units, tiles, &tile_index);

  record->varID   = (short)varID;
  record->levelID = (short)levelID;

  varDefCompType(varID, comptype);

  /*
    Get the ensemble Info from the grib-2 Tables and update the intermediate datastructure.
    Further update to the "vlist" is handled in the same way as for GRIB-1 by "cdi_generate_vars"
  */
  if ( grib_get_long(gh, "typeOfEnsembleForecast", &ens_forecast_type) == 0 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfForecastsInEnsemble", &ens_count ), 0);
      GRIB_CHECK(grib_get_long(gh, "perturbationNumber", &ens_index ), 0);
    }

  if ( ens_index > 0 )
    varDefEnsembleInfo(varID, (int)ens_index, (int)ens_count, (int)ens_forecast_type);

  long typeOfGeneratingProcess = 0;
  if ( grib_get_long(gh, "typeOfGeneratingProcess", &typeOfGeneratingProcess) == 0 )
    varDefTypeOfGeneratingProcess(varID, (int) typeOfGeneratingProcess);

  long productDefinitionTemplate = 0;
  if ( grib_get_long(gh, "productDefinitionTemplateNumber", &productDefinitionTemplate) == 0 )
    varDefProductDefinitionTemplate(varID, (int) productDefinitionTemplate);

  long   lval;
  double dval;

  if (lread_additional_keys)
    for ( int i = 0; i < cdiNAdditionalGRIBKeys; i++ )
      {
        /* note: if the key is not defined, we do not throw an error! */
        if ( grib_get_long(gh, cdiAdditionalGRIBKeys[i], &lval) == 0 )
          varDefOptGribInt(varID, tile_index, lval, cdiAdditionalGRIBKeys[i]);
        if ( grib_get_double(gh, cdiAdditionalGRIBKeys[i], &dval) == 0 )
          varDefOptGribDbl(varID, tile_index, dval, cdiAdditionalGRIBKeys[i]);
      }

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      long center, subcenter;
      GRIB_CHECK(grib_get_long(gh, "centre", &center), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter), 0);
      int instID = institutInq((int)center, (int)subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef((int)center, (int)subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      long processID;
      if ( grib_get_long(gh, "generatingProcessIdentifier", &processID) == 0 )
	{
          /* FIXME: assert(processID >= INT_MIN && processID <= INT_MAX) */
	  int modelID = modelInq(varInqInst(varID), (int)processID, NULL);
	  if ( modelID == CDI_UNDEFID )
	    modelID = modelDef(varInqInst(varID), (int)processID, NULL);
	  varDefModel(varID, modelID);
	}
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);

      if ( pdis == 255 )
	{
	  int tabnum = pcat;
	  int tableID = tableInq(varInqModel(varID), tabnum, NULL);
	  if ( tableID == CDI_UNDEFID )
	    tableID = tableDef(varInqModel(varID), tabnum, NULL);
	  varDefTable(varID, tableID);
	}
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d  param = %d  zaxistype = %d  gridID = %d  levelID = %d",
	    varID, param, zaxistype, gridID, levelID);
}

static compvar2_t gribapiVarSet(int param, int level1, int level2, int leveltype,
                                int tsteptype, char *name, var_tile_t tiles_data)
{
  compvar2_t compVar;
  size_t maxlen = sizeof(compVar.name);
  size_t len = strlen(name);
  if ( len > maxlen ) len = maxlen;

  compVar.param     = param;
  compVar.level1    = level1;
  compVar.level2    = level2;
  compVar.ltype     = leveltype;
  compVar.tsteptype = tsteptype;
  memset(compVar.name, 0, maxlen);
  memcpy(compVar.name, name, len);
  compVar.tiles = tiles_data;

  return compVar;
}

static
int gribapiVarCompare(compvar2_t compVar, record_t record, int flag)
{
  compvar2_t compVar0;
  compVar0.param     = record.param;
  compVar0.level1    = record.ilevel;
  compVar0.level2    = record.ilevel2;
  compVar0.ltype     = record.ltype;
  compVar0.tsteptype = record.tsteptype;
  memcpy(compVar0.name, record.varname, sizeof(compVar.name));

  if ( flag == 0 )
    {
      if ( compVar0.tsteptype == TSTEP_INSTANT  && compVar.tsteptype == TSTEP_INSTANT3 ) compVar0.tsteptype = TSTEP_INSTANT3;
      if ( compVar0.tsteptype == TSTEP_INSTANT3 && compVar.tsteptype == TSTEP_INSTANT  ) compVar0.tsteptype = TSTEP_INSTANT;
    }

  compVar0.tiles = record.tiles;

  return memcmp(&compVar0, &compVar, sizeof(compvar2_t));
}

static
void ensureBufferSize(size_t requiredSize, size_t *curSize, void **buffer)
{
  if ( *curSize < requiredSize )
    {
      *curSize = requiredSize;
      *buffer = Realloc(*buffer, *curSize);
    }
}

static
grib_handle *gribapiGetDiskRepresentation(size_t recsize, size_t *buffersize, void **gribbuffer, int *outDatatype, int *outCompressionType, size_t *outUnzipsize)
{
  int gribversion = (int)((char*)*gribbuffer)[7];

  if ( gribversion <= 1 )
    {
      if ( gribGetZip(recsize, *gribbuffer, outUnzipsize) > 0 )
        {
          *outCompressionType = CDI_COMPRESS_SZIP;
          ensureBufferSize(*outUnzipsize + 100, buffersize, gribbuffer);
        }
      else
        {
          *outCompressionType = CDI_COMPRESS_NONE;
        }
    }

  grib_handle *gh = grib_handle_new_from_message(NULL, *gribbuffer, recsize);

  bool lieee = false;

  if ( gribversion > 1 )
    {
      size_t len = 256;
      char typeOfPacking[256];

      if ( grib_get_string(gh, "packingType", typeOfPacking, &len) == 0 )
        {
          // fprintf(stderr, "packingType %d %s\n", len, typeOfPacking);
          if      ( strncmp(typeOfPacking, "grid_jpeg", len) == 0 ) *outCompressionType = CDI_COMPRESS_JPEG;
          else if ( strncmp(typeOfPacking, "grid_ccsds", len) == 0 ) *outCompressionType = CDI_COMPRESS_SZIP;
          else if ( strncmp(typeOfPacking, "grid_ieee", len) == 0 ) lieee = true;
        }
    }

  if ( lieee )
    {
      *outDatatype = CDI_DATATYPE_FLT64;
      long precision;
      int status = grib_get_long(gh, "precision", &precision);
      if ( status == 0 && precision == 1 ) *outDatatype = CDI_DATATYPE_FLT32;
    }
  else
    {
      *outDatatype = CDI_DATATYPE_PACK;
      long bitsPerValue;
      if ( grib_get_long(gh, "bitsPerValue", &bitsPerValue) == 0 )
        {
          if ( bitsPerValue > 0 && bitsPerValue <= 32 ) *outDatatype = (int)bitsPerValue;
        }
    }

  return gh;
}

typedef enum { CHECKTIME_OK, CHECKTIME_SKIP, CHECKTIME_STOP, CHECKTIME_INCONSISTENT } checkTimeResult;
static checkTimeResult checkTime(stream_t* streamptr, compvar2_t compVar, const DateTime* verificationTime, const DateTime* expectedVTime) {
  // First determine whether the current record exists already.
  int recID = 0;
  for ( ; recID < streamptr->nrecs; recID++ )
    {
      if ( gribapiVarCompare(compVar, streamptr->tsteps[0].records[recID], 1) == 0 ) break;
    }
  int recordExists = recID < streamptr->nrecs;

  // Then we need to know whether the verification time is consistent.
  int consistentTime = !memcmp(verificationTime, expectedVTime, sizeof(*verificationTime));

  // Finally, we make a decision.
  if ( cdiInventoryMode == 1 )
    {
      if ( recordExists ) return CHECKTIME_STOP;
      if ( !consistentTime ) return CHECKTIME_INCONSISTENT;
    }
  else
    {
      if ( !consistentTime ) return CHECKTIME_STOP;
      if ( recordExists ) return CHECKTIME_SKIP;
    }

  return CHECKTIME_OK;
}

#define gribWarning(text, nrecs, timestep, varname, param, level1, level2) do \
  { \
    char paramstr[32]; \
    cdiParamToString(param, paramstr, sizeof(paramstr)); \
    Warning("Record %2d (name=%s id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, varname, paramstr, level1, level2, timestep, text); \
  } \
while(0)

int gribapiScanTimestep1(stream_t * streamptr)
{
  off_t recpos = 0;
  void *gribbuffer = NULL;
  size_t buffersize = 0;
  DateTime datetime0 = { .date = 10101, .time = 0 };
  int nrecs_scanned = 0;        //Only used for debug output.
  bool warn_time = true;
  // bool warn_numavg = true;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  grib_handle *gh = NULL;

  streamptr->curTsID = 0;

  int tsID  = tstepsNewEntry(streamptr);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  int fileID = streamptr->fileID;

  unsigned nrecs = 0;
  while ( true )
    {
      int level1 = 0, level2 = 0;
      size_t recsize = gribGetSize(fileID);
      recpos = fileGetPos(fileID);

      if ( recsize == 0 )
        {
          streamptr->ntsteps = 1;
          break;
        }
      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      int rstatus = gribRead(fileID, gribbuffer, &readsize); //Search for next 'GRIB', read the following record, and position file offset after it.
      if ( rstatus ) break;

      int datatype, comptype = 0;
      size_t unzipsize;
      gh = gribapiGetDiskRepresentation(recsize, &buffersize, &gribbuffer, &datatype, &comptype, &unzipsize);

      nrecs_scanned++;
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      int param = gribapiGetParam(gh);
      int leveltype1 = -1, leveltype2 = -1, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));

      int tsteptype = gribapiGetTsteptype(gh);

      int vdate = 0, vtime = 0;
      gribapiGetValidityDateTime(gh, &vdate, &vtime);
      DateTime datetime = { .date = vdate, .time = vtime };
      /*
      printf("%d %d %d\n", vdate, vtime, leveltype1);
      */

      if ( datetime0.date == 10101 && datetime0.time == 0 )
        {
          if( datetimeCmp(datetime, datetime0) || !nrecs )       //Do we really need this condition? I have included it in order not to change the number of times gribapiGetDataDateTime() etc. get called. But if those are sideeffect-free, this condition should be removed.
            {
              datetime0 = datetime;

              gribapiGetDataDateTime(gh, &rdate, &rtime);

              fcast = gribapiTimeIsFC(gh);
              if ( fcast ) tunit = gribapiGetTimeUnits(gh);
            }
        }

      if ( nrecs )
        {
          checkTimeResult result = checkTime(streamptr, gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles), &datetime, &datetime0);
          if ( result == CHECKTIME_STOP )
            {
              break;
            }
          else if ( result == CHECKTIME_SKIP )
            {
              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);
              continue;
            }
          else if ( result == CHECKTIME_INCONSISTENT && warn_time )
            {
              gribWarning("Inconsistent verification time!", nrecs_scanned, tsID+1, varname, param, level1, level2);
              warn_time = false;
            }
          assert(result == CHECKTIME_OK || result == CHECKTIME_INCONSISTENT);
        }
      /*
      if ( ISEC1_AvgNum )
        {
          if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
            {
              Message("Change numavg from %d to %d not allowed!",
                      taxis->numavg, ISEC1_AvgNum);
              warn_numavg = false;
            }
          else
            {
              taxis->numavg = ISEC1_AvgNum;
            }
        }
      */
      nrecs++;

      if ( CDI_Debug )
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4u %8d name=%s id=%s ltype=%d lev1=%d lev2=%d vdate=%d vtime=%d",
                nrecs, (int)recpos, varname, paramstr, leveltype1, level1, level2, vdate, vtime);
        }

      var_tile_t *ptiles = NULL;
      if ( memcmp(&tiles, &dummy_tiles, sizeof(var_tile_t)) != 0 ) ptiles = &tiles;
      gribapiAddRecord(streamptr, param, gh, recsize, recpos, datatype, comptype, varname,
                       leveltype1, leveltype2, lbounds, level1, level2, level_sf, level_unit, ptiles, 1);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if ( gh ) grib_handle_delete(gh);

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return CDI_EUFSTRUCT;

  cdi_generate_vars(streamptr);

  int taxisID = -1;
  if ( fcast )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
      taxis->type  = TAXIS_RELATIVE;
      taxis->rdate = rdate;
      taxis->rtime = rtime;
      taxis->unit  = tunit;
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      taxis->type  = TAXIS_ABSOLUTE;
    }

  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  int nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
        (record_t *) Realloc(streamptr->tsteps[0].records, (size_t)nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) Malloc((size_t)nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( int recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
        Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = true;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
        {
          streamptr->ntsteps = 0;
          for ( int varID = 0; varID < streamptr->nvars; varID++ )
            vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
        }
    }

  return 0;
}


int gribapiScanTimestep2(stream_t * streamptr)
{
  int rstatus = 0;
  off_t recpos = 0;
  DateTime datetime0 = { LONG_MIN, LONG_MIN };
  // int gridID;
  int recID;
  //  bool warn_numavg = true;
  grib_handle *gh = NULL;

  streamptr->curTsID = 1;

  int fileID  = streamptr->fileID;
  int vlistID = streamptr->vlistID;
  int taxisID = vlistInqTaxis(vlistID);

  void *gribbuffer = streamptr->record->buffer;
  size_t buffersize = streamptr->record->buffersize;

  int tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  int nrecords = streamptr->tsteps[tsID].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) Malloc((size_t)nrecords*sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      streamptr->tsteps[tsID].records[recID].position = streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     = streamptr->tsteps[0].records[recID].size;
    }

  int nrecs_scanned = nrecords; //Only used for debug output
  int rindex = 0;
  while ( true )
    {
      if ( rindex > nrecords ) break;

      size_t recsize = gribGetSize(fileID);
      recpos = fileGetPos(fileID);
      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      size_t unzipsize;
      if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
        ensureBufferSize(unzipsize + 100, &buffersize, &gribbuffer);

      nrecs_scanned++;
      gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      int param = gribapiGetParam(gh);
      int level1 = 0, level2 = 0, leveltype1, leveltype2, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));

      int vdate = 0, vtime = 0;
      gribapiGetValidityDateTime(gh, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;

              gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

	      taxis->unit  = gribapiGetTimeUnits(gh);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      int tsteptype = gribapiGetTsteptype(gh);
      /*
      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
		(taxis->numavg != ISEC1_AvgNum) )
	    {
	      warn_numavg = false;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}
      */
      DateTime datetime = {
        .date = vdate,
        .time = vtime
      };

      compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles);

      for ( recID = 0; recID < nrecords; recID++ )
        if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;

      if ( recID == nrecords )
	{
	  gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, varname, param, level1, level2);
	  return CDI_EUFSTRUCT;
	}

      if ( streamptr->tsteps[tsID].records[recID].used )
        {
          if ( cdiInventoryMode == 1 ) break;
          else
	    {
	      if ( datetimeCmp(datetime, datetime0) != 0 ) break;

              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);
	      continue;
	    }
	}

      streamptr->tsteps[tsID].records[recID].used = true;
      streamptr->tsteps[tsID].recIDs[rindex] = recID;

      if ( CDI_Debug )
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4d %8d name=%s id=%s ltype=%d lev1=%d lev2=%d vdate=%d vtime=%d",
                  nrecs_scanned, (int)recpos, varname, paramstr, leveltype1, level1, level2, vdate, vtime);
        }

      streamptr->tsteps[tsID].records[recID].size = recsize;

      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return CDI_EUFSTRUCT;
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      int varID = streamptr->tsteps[tsID].records[recID].varID;
      /*
      gridID = vlistInqVarGrid(vlistID, varID);
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}
      */
      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      grib_handle_delete(gh);
      gh = NULL;

      rindex++;
    }

  if ( gh ) grib_handle_delete(gh);

  int nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  int varID = streamptr->tsteps[tsID].records[recID].varID;
	  vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = true;
      streamptr->tsteps[tsID].position = recpos;
    }

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  return rstatus;
}


int gribapiScanTimestep(stream_t * streamptr)
{
  int vrecID, recID;
  //bool warn_numavg = true;
  int nrecs = 0;
  int vlistID = streamptr->vlistID;

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      void *gribbuffer = streamptr->record->buffer;
      size_t buffersize = streamptr->record->buffersize;

      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) Malloc((size_t)nrecs*sizeof(int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      int nrecs_scanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs*(tsID-1);    //Only used for debug output.
      int rindex = 0;
      off_t recpos = 0;
      DateTime datetime0 = { LONG_MIN, LONG_MIN };
      grib_handle *gh = NULL;
      char varname[256];
      while ( true )
	{
	  if ( rindex > nrecs ) break;

	  size_t recsize = gribGetSize(fileID);
	  recpos = fileGetPos(fileID);
	  if ( recsize == 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }

	  if ( rindex >= nrecs ) break;

          ensureBufferSize(recsize, &buffersize, &gribbuffer);

	  size_t readsize = recsize;
	  if (gribRead(fileID, gribbuffer, &readsize))
	    {
	      Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
		      streamptr->tsteps[tsID].recordSize);
	      break;
	    }

          size_t unzipsize;
	  if ( gribGetZip(recsize, gribbuffer, &unzipsize) > 0 )
            ensureBufferSize(unzipsize + 100, &buffersize, &gribbuffer);

          nrecs_scanned++;
	  gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
	  GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

          int param = gribapiGetParam(gh);
          int level1 = 0, level2 = 0, leveltype1, leveltype2 = -1, lbounds, level_sf, level_unit;
          var_tile_t tiles = dummy_tiles;
          gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

          varname[0] = 0;
	  gribapiGetString(gh, "shortName", varname, sizeof(varname));

          int vdate = 0, vtime = 0;
	  gribapiGetValidityDateTime(gh, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
              int taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;

                  gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

		  taxis->unit  = gribapiGetTimeUnits(gh);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }
	  /*
	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
		  warn_numavg = false;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }
	  */
          DateTime datetime = {
            .date  = vdate,
            .time  = vtime
          };

          int tsteptype = gribapiGetTsteptype(gh);

          compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles);

	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, varname, param, level1, level2);

	      if ( cdiInventoryMode == 1 )
		return CDI_EUFSTRUCT;
	      else
		continue;
	    }

	  if ( cdiInventoryMode != 1 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  if ( datetimeCmp(datetime, datetime0) != 0 ) break;

		  if ( CDI_Debug )
                    gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);

		  continue;
		}
	    }

          streamptr->tsteps[tsID].records[recID].used = true;
          streamptr->tsteps[tsID].recIDs[rindex] = recID;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

	  if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, level1);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex, (int)recpos, param, level1, vdate, vtime);

	  grib_handle_delete(gh);
	  gh = NULL;

	  rindex++;
	}

      if ( gh ) grib_handle_delete(gh);

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  gribWarning("Paramameter not found!", nrecs_scanned, tsID+1, varname, streamptr->tsteps[tsID].records[recID].param,
                      streamptr->tsteps[tsID].records[recID].ilevel, streamptr->tsteps[tsID].records[recID].ilevel2);
	  return CDI_EUFSTRUCT;
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = true;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;

      streamptr->record->buffer     = gribbuffer;
      streamptr->record->buffersize = buffersize;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return (int)streamptr->ntsteps;
}

#ifdef gribWarning
#undef gribWarning
#endif

int gribapiDecode(void *gribbuffer, int gribsize, double *data, long gridsize,
		  int unreduced, int *nmiss, double missval, int vlistID, int varID)
{
  int status = 0;
  long lpar;
  long numberOfPoints;
  size_t datasize;

  UNUSED(vlistID);
  UNUSED(varID);

  if ( unreduced )
    {
      static bool lwarn = true;

      if ( lwarn )
	{
	  lwarn = false;
	  Warning("Conversion of gaussian reduced grids unsupported!");
	}
    }

  size_t recsize = (size_t)gribsize;
  grib_handle *gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
  GRIB_CHECK(my_grib_set_double(gh, "missingValue", missval), 0);

  /* get the size of the values array*/
  GRIB_CHECK(grib_get_size(gh, "values", &datasize), 0);
  GRIB_CHECK(grib_get_long(gh, "numberOfPoints", &numberOfPoints), 0);

  // printf("values_size = %d  numberOfPoints = %ld\n", datasize, numberOfPoints);

  if ( gridsize != (long) datasize )
    Error("Internal problem: gridsize(%ld) != datasize(%zu)!", gridsize, datasize);
  size_t dummy = datasize;
  GRIB_CHECK(grib_get_double_array(gh, "values", data, &dummy), 0);

  GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
  int gridtype = (int) lpar;

  *nmiss = 0;
  if ( gridtype < 50 || gridtype > 53 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfMissing", &lpar), 0);
      *nmiss = (int) lpar;
      // printf("gridtype %d, nmiss %d\n", gridtype, nmiss);
    }

  grib_handle_delete(gh);

  return status;
}


static
void gribapiDefInstitut(grib_handle *gh, int vlistID, int varID)
{
  int instID;

  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID )
    instID = vlistInqInstitut(vlistID);
  else
    instID = vlistInqVarInstitut(vlistID, varID);

  if ( instID != CDI_UNDEFID )
    {
      long center, subcenter;
      long center0, subcenter0;

      center    = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);

      GRIB_CHECK(grib_get_long(gh, "centre", &center0), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter0), 0);

      if ( center != center0 )
	GRIB_CHECK(my_grib_set_long(gh, "centre", center), 0);
      if ( subcenter != subcenter0 )
	GRIB_CHECK(my_grib_set_long(gh, "subCentre", subcenter), 0);
    }
}

static
void gribapiDefModel(grib_handle *gh, int vlistID, int varID)
{
  int modelID;

  if ( vlistInqModel(vlistID) != CDI_UNDEFID )
    modelID = vlistInqModel(vlistID);
  else
    modelID = vlistInqVarModel(vlistID, varID);

  if ( modelID != CDI_UNDEFID )
    GRIB_CHECK(my_grib_set_long(gh, "generatingProcessIdentifier", modelInqGribID(modelID)), 0);
}

static
void gribapiDefParam(int editionNumber, grib_handle *gh, int param, const char *name, const char *stdname)
{
  bool ldefined = false;

  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if ( pnum < 0 )
    {
      size_t len;
      len = strlen(stdname);
      if ( len )
        {
          int status = my_grib_set_string(gh, "cfName", stdname, &len);
          if ( status == 0 ) ldefined = true;
          else Warning("grib_api: No match for cfName=%s", stdname);
        }

      if ( ldefined == false )
        {
          len = strlen(name);
          int status = my_grib_set_string(gh, "shortName", name, &len);
          if ( status == 0 ) ldefined = true;
          else Warning("grib_api: No match for shortName=%s", name);
        }
    }

  if ( ldefined == false )
    {
      if ( pnum < 0 ) pnum = -pnum;

      static bool lwarn_pnum = true;
      if ( pnum > 255 && lwarn_pnum )
        {
          Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum%256);
          lwarn_pnum = false;
          pnum = pnum%256;
        }

      if ( editionNumber <= 1 )
	{
          static bool lwarn_pdis = true;
	  if ( pdis != 255 && lwarn_pdis )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));
	      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
              lwarn_pdis = false;
	    }

	  GRIB_CHECK(my_grib_set_long(gh, "table2Version",        pcat), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "indicatorOfParameter", pnum), 0);
	}
      else
	{
	  GRIB_CHECK(my_grib_set_long(gh, "discipline",        pdis), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "parameterCategory", pcat), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "parameterNumber",   pnum), 0);
	}
    }

  // printf("param: %d.%d.%d %s\n", pnum, pcat, pdis, name);
}

static
int getTimeunitFactor(int timeunit)
{
  int factor = 1;

  switch (timeunit)
    {
    case TUNIT_SECOND:  factor =     1;  break;
    case TUNIT_MINUTE:  factor =    60;  break;
    case TUNIT_HOUR:    factor =  3600;  break;
    case TUNIT_3HOURS:  factor = 10800;  break;
    case TUNIT_6HOURS:  factor = 21600;  break;
    case TUNIT_12HOURS: factor = 43200;  break;
    case TUNIT_DAY:     factor = 86400;  break;
    default:            factor =  3600;  break;
    }

  return factor;
}

static
int grib2ProDefTempHasStatisticalDef(int proDefTempNum)
{
  int hasStatisticalDef = 0;

  switch (proDefTempNum)
    {
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
      case 13:
      case 14:
      case 34:
      case 42:
      case 43:
      case 46:
      case 47:
      case 61:
      case 91:
      case 1001:
      case 1101:
      case 40034:
               hasStatisticalDef = 1;  break;
      default: hasStatisticalDef = 0;  break;
    }

  return hasStatisticalDef;
}

static
void gribapiDefStepUnits(int editionNumber, grib_handle *gh, int timeunit, int proDefTempNum, int gcinit)
{
  long unitsOfTime;

  switch (timeunit)
    {
    case TUNIT_SECOND:  unitsOfTime = 13;  break;
    case TUNIT_MINUTE:  unitsOfTime =  0;  break;
    case TUNIT_HOUR:    unitsOfTime =  1;  break;
    case TUNIT_3HOURS:  unitsOfTime = 10;  break;
    case TUNIT_6HOURS:  unitsOfTime = 11;  break;
    case TUNIT_12HOURS: unitsOfTime = 12;  break;
    case TUNIT_DAY:     unitsOfTime =  2;  break;
    default:            unitsOfTime =  1;  break;
    }

  if ( !gcinit )
    {
      GRIB_CHECK(my_grib_set_long(gh, "stepUnits", unitsOfTime), 0);
      if ( editionNumber == 1 )
        {
          GRIB_CHECK(my_grib_set_long(gh, "unitOfTimeRange", unitsOfTime), 0);
        }
      else if ( grib2ProDefTempHasStatisticalDef(proDefTempNum) )
        {
          GRIB_CHECK(my_grib_set_long(gh, "indicatorOfUnitForTimeRange", unitsOfTime), 0);
          GRIB_CHECK(my_grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime), 0);
        }
      else
        {
	  // NOTE KNMI:  HIRLAM model files LAMH_D11 are in grib1 and do NOT have key indicatorOfUnitForTimeRange
	  // Watch out for compatibility issues.
          GRIB_CHECK(my_grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime), 0);
        }
    }
}

static
int gribapiDefSteptype(int editionNumber, grib_handle *gh, int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int gcinit)
{
  long proDefTempNum = 0;
  size_t len = 64;
  const char *stepType;

  if (tsteptype >= TSTEP_INSTANT && tsteptype <= TSTEP_RATIO)
    {
      stepType = cdiGribAPI_ts_str_map[tsteptype].sname;
      proDefTempNum = cdiGribAPI_ts_str_map[tsteptype].productionTemplate;
    }
  else
    {
      stepType = "instant";
      proDefTempNum = 0;
    }

  if ( typeOfGeneratingProcess == 4 )
    {
      if ( proDefTempNum == 8 ) proDefTempNum = 11;
      else                      proDefTempNum = 1;
    }

  if ( productDefinitionTemplate != -1 ) proDefTempNum = productDefinitionTemplate;

  if ( !gcinit )
    {
      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "productDefinitionTemplateNumber", proDefTempNum), 0);
      len = strlen(stepType);
      GRIB_CHECK(my_grib_set_string(gh, "stepType", stepType, &len), 0);
    }

  return (int)proDefTempNum;
}

static
void gribapiDefDateTimeAbs(int editionNumber, grib_handle *gh, int date, int time, int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int gcinit)
{
  (void ) gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);

  if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 0), 0);
  if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

  if ( date == 0 ) date = 10101;
  gribapiSetDataDateTime(gh, date, time);
}

static
int gribapiDefDateTimeRel(int editionNumber, grib_handle *gh, int rdate, int rtime, int vdate, int vtime,
                          int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int timeunit, int calendar, int gcinit)
{
  int status = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  long startStep = 0, endStep;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday1, &secofday1);

  if ( vdate == 0 && vtime == 0 ) { vdate = rdate; vtime = rtime; }

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  int factor = getTimeunitFactor(timeunit);

  if ( !(int)(fmod(days*86400.0 + secs, factor)))
    {
      int proDefTempNum = gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);

      gribapiDefStepUnits(editionNumber, gh, timeunit, proDefTempNum, gcinit);

      endStep = (int) ((days*86400.0 + secs)/factor);

      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 1), 0);
      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

      if ( rdate == 0 ) rdate = 10101;
      gribapiSetDataDateTime(gh, rdate, rtime);

      // printf(">>>>> tsteptype %d  startStep %ld  endStep %ld\n", tsteptype, startStep, endStep);

      // Product Definition Template Number: defined in GRIB_API file 4.0.table
      // point in time products:
      if ( (proDefTempNum >= 0 && proDefTempNum <=  7) || 
           proDefTempNum == 55 || proDefTempNum == 40055 ) // Tile
        startStep = endStep;

      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "forecastTime", startStep), 0);
      GRIB_CHECK(my_grib_set_long(gh, "endStep", endStep), 0);

      status = 0;
    }

  return status;
}

static
void gribapiDefTime(int editionNumber, int productDefinitionTemplate, int typeOfGeneratingProcess, grib_handle *gh,
                    int vdate, int vtime, int tsteptype, int numavg, int taxisID, int gcinit)
{
  int taxistype = -1;

  UNUSED(numavg);

  if ( taxisID != -1 ) taxistype = taxisInqType(taxisID);

  if ( typeOfGeneratingProcess == 196 )
    {
      vdate = 10101;
      vtime = 0;
      taxistype = TAXIS_ABSOLUTE;
    }
  /*
  else if ( typeOfGeneratingProcess == 9 )
    {
    }
  */

  if ( taxistype == TAXIS_RELATIVE )
    {
      int status;
      int calendar = taxisInqCalendar(taxisID);
      int rdate    = taxisInqRdate(taxisID);
      int rtime    = taxisInqRtime(taxisID);
      int timeunit = taxisInqTunit(taxisID);

      status = gribapiDefDateTimeRel(editionNumber, gh, rdate, rtime, vdate, vtime,
                                     productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, timeunit, calendar, gcinit);

      if ( status != 0 ) taxistype = TAXIS_ABSOLUTE;
    }

  if ( taxistype == TAXIS_ABSOLUTE )
    {
      gribapiDefDateTimeAbs(editionNumber, gh, vdate, vtime, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);
    }
}

struct gribApiMsg {
  size_t msgLen;
  const char *msg;
};

static struct gribApiMsg
getGribApiCompTypeMsg(int comptype, int gridsize)
{
  const char *mesg;
  size_t len;

  if ( comptype == CDI_COMPRESS_JPEG && gridsize > 1 )
    {
      static const char mesg_grid_jpeg[] = "grid_jpeg";
      len = sizeof (mesg_grid_jpeg) - 1;
      mesg = mesg_grid_jpeg;
    }
  else if ( comptype == CDI_COMPRESS_SZIP && gridsize > 1 )
    {
      static const char mesg_grid_ccsds[] = "grid_ccsds";
      len = sizeof (mesg_grid_ccsds) - 1;
      mesg = mesg_grid_ccsds;
    }
  else
    {
      static const char mesg_simple[] = "grid_simple";
      len = sizeof (mesg_simple) - 1;
      mesg = mesg_simple;
    }

  return (struct gribApiMsg){ .msgLen = len, .msg = mesg };
}


static
void gribapiDefGrid(int editionNumber, grib_handle *gh, int gridID, int comptype, bool lieee, int datatype, int nmiss, int gcinit)
{
  UNUSED(nmiss);
  bool lrotated = false;
  bool lcurvi = false;

  int gridtype = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);

  if ( editionNumber <= 1 )
    if ( gridtype == GRID_GME || gridtype == GRID_UNSTRUCTURED )
      gridtype = -1;

  if ( gridtype == GRID_GENERIC )
    {
      int xsize = gridInqXsize(gridID);
      int ysize = gridInqYsize(gridID);

      if ( (ysize ==  32 || ysize ==  48 || ysize ==  64 ||
	    ysize ==  96 || ysize == 160 || ysize == 192 ||
	    ysize == 240 || ysize == 320 || ysize == 384 ||
	    ysize == 480 || ysize == 768 ) &&
	   (xsize == 2*ysize || xsize == 1) )
	{
	  gridtype = GRID_GAUSSIAN;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridsize == 1 )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      int projID = gridInqProj(gridID);
      if ( projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION )
        {
          gridID = projID;
          gridtype = GRID_PROJECTION;
        }
      else
        {
          static bool lwarning = true;
          if ( lwarning && gridsize > 1 )
            {
              lwarning = false;
              Warning("Curvilinear grid is unsupported in GRIB format! Created wrong Grid Description Section!");
            }
          lcurvi = true;
          gridtype = GRID_LONLAT;
        }
    }

  if ( gridtype == GRID_PROJECTION )
    {
      if ( gridInqProjType(gridID) == CDI_PROJ_RLL )
        {
          gridtype = GRID_LONLAT;
          lrotated = true;
        }
      else if ( gridInqProjType(gridID) == CDI_PROJ_LCC )
        {
          gridtype = GRID_LCC;
        }
    }

  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
    {
      if ( editionNumber != 2 || lieee ) { comptype = 0; }

      if ( comptype )
        {
          struct gribApiMsg gaMsg = getGribApiCompTypeMsg(comptype, gridsize);
          size_t len = gaMsg.msgLen;
          const char *mesg = gaMsg.msg;
          GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
        }
    }

  if ( gcinit ) return;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
	double xfirst = 0, xlast = 0, xinc = 0;
	double yfirst = 0, ylast = 0, yinc = 0;

        const char *mesg;
        size_t len;
        if ( gridtype == GRID_GAUSSIAN )
          {
            static const char mesg_gaussian[] = "regular_gg";
            len = sizeof(mesg_gaussian) - 1;
            mesg = mesg_gaussian;
          }
        else if ( gridtype == GRID_GAUSSIAN_REDUCED )
          {
            static const char mesg_gaussian_reduced[] = "reduced_gg";
            len = sizeof(mesg_gaussian_reduced) - 1;
            mesg = mesg_gaussian_reduced;
          }
        else if ( lrotated )
          {
            static const char mesg_rot_lonlat[] = "rotated_ll";
            len = sizeof(mesg_rot_lonlat) - 1;
            mesg = mesg_rot_lonlat;
          }
        else
          {
            static const char mesg_regular_ll[] = "regular_ll";
            len = sizeof(mesg_regular_ll) - 1;
            mesg = mesg_regular_ll;
          }
        GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

	int nlon = gridInqXsize(gridID);
	int nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    nlon = 0;

	    int *rowlon = (int *) Malloc((size_t)nlat*sizeof(int));
	    long *pl    = (long *) Malloc((size_t)nlat*sizeof(long));
	    gridInqRowlon(gridID, rowlon);
	    for ( int i = 0; i < nlat; ++i ) pl[i] = rowlon[i];

            GRIB_CHECK(grib_set_long_array(gh, "pl", pl, (size_t)nlat), 0);

	    Free(pl);
	    Free(rowlon);

	    xfirst = 0;
	    xinc   =        360. * 0.5 / (double)nlat;
	    xlast  = 360. - 360. * 0.5 / (double)nlat;
	  }
	else
	  {
	    if ( nlon == 0 ) nlon = 1;
	    else
	      {
		xfirst = gridInqXval(gridID, 0);
                xlast  = gridInqXval(gridID, (lcurvi ? nlon*nlat : nlon) - 1);
		xinc   = fabs(gridInqXinc(gridID));
	      }
	  }

	if ( nlat == 0 ) nlat = 1;
	else
	  {
	    yfirst = gridInqYval(gridID, 0);
            ylast  = gridInqYval(gridID, (lcurvi ? nlon*nlat : nlat) - 1);
	    yinc   = fabs(gridInqYinc(gridID));
	  }

	if ( gridtype != GRID_GAUSSIAN_REDUCED ) GRIB_CHECK(my_grib_set_long(gh, "Ni", nlon), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", xfirst), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfLastGridPointInDegrees",  xlast), 0);
	GRIB_CHECK(my_grib_set_double(gh, "iDirectionIncrementInDegrees", xinc), 0);

	GRIB_CHECK(my_grib_set_long(gh, "Nj", nlat), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees",  yfirst), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfLastGridPointInDegrees",   ylast), 0);

        {
          long iscan = xfirst > xlast;
          GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", iscan), 0);
        }
        {
          long jscan = yfirst < ylast;
          GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", jscan), 0);
        }

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            GRIB_CHECK(my_grib_set_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", np), 0);
          }
	else
	  {
	    GRIB_CHECK(my_grib_set_double(gh, "jDirectionIncrementInDegrees", yinc), 0);
	  }

	if ( lrotated )
	  {
            double xpole = 0, ypole = 0, angle = 0;
            gridInqParamRLL(gridID, &xpole, &ypole, &angle);

            xpole += 180;
            if ( fabs(ypole) > 0 ) ypole = -ypole; // change from north to south pole
            if ( fabs(angle) > 0 ) angle = -angle;
            GRIB_CHECK(my_grib_set_double(gh, "latitudeOfSouthernPoleInDegrees",  ypole), 0);
            GRIB_CHECK(my_grib_set_double(gh, "longitudeOfSouthernPoleInDegrees", xpole), 0);
            GRIB_CHECK(my_grib_set_double(gh, "angleOfRotation", angle), 0);
          }

        if ( editionNumber != 2 ) { lieee = false; comptype = 0; }

        if ( lieee )
          {
            static const char mesg_grid_ieee[] = "grid_ieee";
            len = sizeof (mesg_grid_ieee) - 1;
            mesg = mesg_grid_ieee;
          }
        else
          {
            struct gribApiMsg gaMsg = getGribApiCompTypeMsg(comptype, gridsize);
            len = gaMsg.msgLen;
            mesg = gaMsg.msg;
          }
        GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
        if ( lieee )
          GRIB_CHECK(my_grib_set_long(gh, "precision", datatype == CDI_DATATYPE_FLT64 ? 2 : 1), 0);

        long uvRelativeToGrid = gridInqUvRelativeToGrid(gridID);
        if ( uvRelativeToGrid ) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);

	break;
      }
    case GRID_LCC:
      {
	int xsize = gridInqXsize(gridID);
	int ysize = gridInqYsize(gridID);

        double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;
	gridInqParamLCC(gridID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);
	gridVerifyGribParamLCC(grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);
        if ( xval_0 < 0 ) xval_0 += 360;
        bool lsouth = (lat_1 < 0);
        if ( lsouth ) { lat_1 = -lat_2; lat_2 = -lat_2; }
        int projflag = 0;
        if ( lsouth ) gribbyte_set_bit(&projflag, 1);

        double xinc = gridInqXinc(gridID);
        double yinc = gridInqYinc(gridID);

        static const char mesg[] = "lambert";
        size_t len = sizeof(mesg) -1;
        GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

	GRIB_CHECK(my_grib_set_long(gh, "Nx", xsize), 0);
	GRIB_CHECK(my_grib_set_long(gh, "Ny", ysize), 0);
	GRIB_CHECK(my_grib_set_long(gh, "DxInMetres", lround(xinc)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "DyInMetres", lround(yinc)), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", xval_0), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", yval_0), 0);
	GRIB_CHECK(my_grib_set_double(gh, "LoVInDegrees", lon_0), 0);
	GRIB_CHECK(my_grib_set_double(gh, "Latin1InDegrees", lat_1), 0);
	GRIB_CHECK(my_grib_set_double(gh, "Latin2InDegrees", lat_2), 0);
        GRIB_CHECK(my_grib_set_long(gh, "projectionCentreFlag", projflag), 0);

        long uvRelativeToGrid = gridInqUvRelativeToGrid(gridID);
        if ( uvRelativeToGrid ) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);
        long earthIsOblate = (IS_EQUAL(a, 6378160.) && IS_EQUAL(rf, 297.));
        if ( earthIsOblate ) GRIB_CHECK(my_grib_set_long(gh, "earthIsOblate", earthIsOblate), 0);

        int scanflag = 0;
        gribbyte_set_bit(&scanflag, 2);
        if ( editionNumber <= 1 )
          GRIB_CHECK(my_grib_set_long(gh, "scanningMode", (long)scanflag), 0);

	break;
      }
    case GRID_SPECTRAL:
      {
        {
          static const char mesg[] = "sh";
          size_t len = sizeof (mesg) -1;
          GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
        }
	{
          int trunc = gridInqTrunc(gridID);
          enum { numTruncAtt = 3 };
          static const char truncAttNames[numTruncAtt][2] = { "J", "K", "M" };
          for (size_t i = 0; i < numTruncAtt; ++i)
            GRIB_CHECK(my_grib_set_long(gh, truncAttNames[i], trunc), 0);
        }
	// GRIB_CHECK(my_grib_set_long(gh, "numberOfDataPoints", gridsize), 0);
        /*
        if ( lieee )
          {
            printf("spectral_ieee\n");
            if ( editionNumber == 2 ) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", gridsize, 0);
            static const char mesg[] = "spectral_ieee";
            size_t len = sizeof (mesg) -1;
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }
        else */ if ( gridInqComplexPacking(gridID) )
	  {
	    if ( editionNumber == 2 ) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", gridsize), 0);
            static const char mesg[] = "spectral_complex";
            size_t len = sizeof (mesg) -1;
	    GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
            enum { numTruncAtt = 3 };
            static const char truncAttNames[numTruncAtt][3]
              = { "JS", "KS", "MS" };
            for (size_t i = 0; i < numTruncAtt; ++i)
              GRIB_CHECK(my_grib_set_long(gh, truncAttNames[i], 20), 0);
	  }
	else
	  {
            static const char mesg[] = "spectral_simple";
            size_t len = sizeof (mesg) -1;
	    GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
	  }

	break;
      }
    case GRID_GME:
      {
	GRIB_CHECK(my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_GME), 0);

        int nd = 0, ni = 0, ni2 = 0, ni3 = 0;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
	GRIB_CHECK(my_grib_set_long(gh, "nd", nd), 0);
	GRIB_CHECK(my_grib_set_long(gh, "Ni", ni), 0);
	GRIB_CHECK(my_grib_set_long(gh, "n2", ni2), 0);
	GRIB_CHECK(my_grib_set_long(gh, "n3", ni3), 0);
	GRIB_CHECK(my_grib_set_long(gh, "latitudeOfThePolePoint", 90000000), 0);
	GRIB_CHECK(my_grib_set_long(gh, "longitudeOfThePolePoint", 0), 0);

	GRIB_CHECK(my_grib_set_long(gh, "numberOfDataPoints", gridsize), 0);
	GRIB_CHECK(my_grib_set_long(gh, "totalNumberOfGridPoints", gridsize), 0);

        if ( comptype == CDI_COMPRESS_SZIP )
          {
            static const char mesg[] = "grid_ccsds";
            size_t len = sizeof (mesg) -1;
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }

	break;
      }
    case GRID_UNSTRUCTURED:
      {
	static bool warning = true;

	int status = my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_UNSTRUCTURED);
	if ( status != 0 && warning )
	  {
	    warning = false;
	    Warning("Can't write reference grid!");
	    Warning("gridDefinitionTemplateNumber %d not found (grib2/template.3.%d.def)!",
		    GRIB2_GTYPE_UNSTRUCTURED, GRIB2_GTYPE_UNSTRUCTURED);
	  }
	else
	  {
            unsigned char uuid[CDI_UUID_SIZE];
            int position = gridInqPosition(gridID);
            int number = gridInqNumber(gridID);
            if ( position < 0 ) position = 0;
            if ( number < 0 ) number = 0;
	    GRIB_CHECK(my_grib_set_long(gh, "numberOfGridUsed", number), 0);
	    GRIB_CHECK(my_grib_set_long(gh, "numberOfGridInReference", position), 0);
            size_t len = CDI_UUID_SIZE;
            gridInqUUID(gridID, uuid);
	    if (grib_set_bytes(gh, "uuidOfHGrid", uuid, &len) != 0)
	      Warning("Can't write UUID!");
	  }

        if ( comptype == CDI_COMPRESS_SZIP )
          {
            static const char mesg[] = "grid_ccsds";
            size_t len = sizeof (mesg) -1;
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }

	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }
}

static
void getLevelFactor(double level, long *factor, long *out_scaled_value)
{
  double scaled_value  = level;
  long   iscaled_value = lround(scaled_value);
  long   i;

  const double eps = 1.e-8;
  for ( i=0; (fabs(scaled_value - (double) iscaled_value) >= eps) && i < 7; i++ )
    {
      scaled_value *= 10.;
      iscaled_value = lround(scaled_value);
    }

  (*factor)           = i;
  (*out_scaled_value) = iscaled_value;
}

static
void gribapiDefLevelType(grib_handle *gh, int gcinit, const char *keyname, long leveltype)
{
  bool lset = false;
  if ( (leveltype == GRIB1_LTYPE_ISOBARIC_PA || leveltype == 99 || leveltype == 100) && gribEditionNumber(gh) == 1 )
    {
      if ( gribGetLong(gh, "indicatorOfTypeOfLevel") != leveltype ) lset = true;
    }

  if ( !gcinit || lset ) GRIB_CHECK(my_grib_set_long(gh, keyname, leveltype), 0);
}

static
void grib1DefLevel(grib_handle *gh, int gcinit, long leveltype, bool lbounds, double level, double dlevel1, double dlevel2)
{
  gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", leveltype);

  if ( lbounds )
    {
      GRIB_CHECK(my_grib_set_long(gh, "topLevel", lround(dlevel1)), 0);
      GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", lround(dlevel2)), 0);
    }
  else
    {
      GRIB_CHECK(my_grib_set_long(gh, "level", lround(level)), 0);
    }
}

static
void grib2DefLevel(grib_handle *gh, int gcinit, long leveltype1, long leveltype2, bool lbounds, double level, double dlevel1, double dlevel2)
{
  gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", leveltype1);
  if ( lbounds ) gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", leveltype2);

  if ( !lbounds ) dlevel1 = level;

  long scaled_level, factor;
  getLevelFactor(dlevel1, &factor, &scaled_level);
  GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfFirstFixedSurface", factor), 0);
  GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfFirstFixedSurface", scaled_level), 0);

  if ( lbounds )
    {
      getLevelFactor(dlevel2, &factor, &scaled_level);
      GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfSecondFixedSurface", factor), 0);
      GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfSecondFixedSurface", scaled_level), 0);
    }
}

static
void gribapiDefLevel(int editionNumber, grib_handle *gh, int zaxisID, int levelID, int gcinit, int proddef_template_num)
{
  char units[CDI_MAX_NAME];
  bool lbounds = false;
  double dlevel1 = 0, dlevel2 = 0;

  int zaxistype = zaxisInqType(zaxisID);
  long ltype = zaxisInqLtype(zaxisID);
  long ltype2 = zaxisInqLtype2(zaxisID);
  double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) : levelID+1;

  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      lbounds = true;
      dlevel1 = zaxisInqLbound(zaxisID, levelID);
      dlevel2 = zaxisInqUbound(zaxisID, levelID);
    }
  else
    {
      dlevel1 = level;
      dlevel2 = 0;
    }

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s", zaxisNamePtr(zaxistype), zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  long grib_ltype = -1;
  if ( editionNumber <= 1 )
    grib_ltype = zaxisTypeToGrib1ltype(zaxistype);
  else
    grib_ltype = zaxisTypeToGrib2ltype(zaxistype);

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_HEIGHT:
    case ZAXIS_ALTITUDE:
    case ZAXIS_SIGMA:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
        if ( zaxistype == ZAXIS_HEIGHT )
          {
            double sf = 1;
            zaxisInqUnits(zaxisID, units);
            if ( units[1] == 'm' && !units[2] )
              {
                if      ( units[0] == 'c' ) sf = 0.01;
                else if ( units[0] == 'd' ) sf = 0.1;
                else if ( units[0] == 'k' ) sf = 1000;
              }
            if ( IS_NOT_EQUAL(sf, 1) )
              {
                level   *= sf;
                dlevel1 *= sf;
                dlevel2 *= sf;
              }
          }

        if ( editionNumber <= 1 )
          {
            grib1DefLevel(gh, gcinit, grib_ltype, lbounds, level, dlevel1, dlevel2);
          }
        else
          {
            /* PRODUCT DEFINITION TEMPLATE NUMBER 32:

               "Analysis or forecast at a horizontal level or in a horizontal layer at a point
                in time for simulate (synthetic) satellite data"

               The key/value pairs that are set in "grib2DefLevel" do not exist for this template.
            */
            if ( proddef_template_num != 32 )
              grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_LAKE_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM_TA:
    case ZAXIS_SEDIMENT_BOTTOM_TW:
    case ZAXIS_MIX_LAYER:
    case ZAXIS_ATMOSPHERE:
      {
        if ( editionNumber <= 1 )
          grib1DefLevel(gh, gcinit, grib_ltype, lbounds, level, dlevel1, dlevel2);
        else
          grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, lbounds, level, dlevel1, dlevel2);

        break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        if ( editionNumber <= 1 )
          {
            grib_ltype = lbounds ? GRIB1_LTYPE_HYBRID_LAYER : GRIB1_LTYPE_HYBRID;
            grib1DefLevel(gh, gcinit, grib_ltype, lbounds, level, dlevel1, dlevel2);
          }
        else
          {
            grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, lbounds, level, dlevel1, dlevel2);
          }

        if ( !gcinit )
          {
            int vctsize = zaxisInqVctSize(zaxisID);
            if ( vctsize > 0 )
              {
                GRIB_CHECK(my_grib_set_long(gh, "PVPresent", 1), 0);
                GRIB_CHECK(grib_set_double_array(gh, "pv", zaxisInqVctPtr(zaxisID), (size_t)vctsize), 0);
              }
          }

	break;
      }
    case ZAXIS_PRESSURE:
      {
	if ( level < 0 ) Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "Pa", 2) != 0 )
          {
            level   *= 100;
            dlevel1 *= 100;
            dlevel2 *= 100;
          }

        if ( editionNumber <= 1 )
          {
            double dum;
            if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
              grib_ltype = GRIB1_LTYPE_ISOBARIC_PA;
            else
              level /= 100;

            grib1DefLevel(gh, gcinit, grib_ltype, lbounds, level, dlevel1, dlevel2);
	  }
	else
	  {
            if ( ltype2 == -1 ) ltype2 = GRIB2_LTYPE_ISOBARIC;
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_ISOBARIC, ltype2, lbounds, level, dlevel1, dlevel2);
	  }

	break;
      }
    case ZAXIS_SNOW:
      {
        if ( editionNumber <= 1 )
          ; // not available
	else
          {
            grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	zaxisInqUnits(zaxisID, units);
        double sf; //scalefactor

	if ( editionNumber <= 1 )
	  {
	    if      ( memcmp(units, "mm", 2) == 0 ) sf =   0.1;
	    else if ( memcmp(units, "cm", 2) == 0 ) sf =   1; // cm
	    else if ( memcmp(units, "dm", 2) == 0 ) sf =  10;
	    else                                    sf = 100;

            grib1DefLevel(gh, gcinit, grib_ltype, lbounds, level*sf, dlevel1*sf, dlevel2*sf);
	  }
	else
	  {
	    if      ( memcmp(units, "mm", 2) == 0 ) sf = 0.001;
	    else if ( memcmp(units, "cm", 2) == 0 ) sf = 0.01;
	    else if ( memcmp(units, "dm", 2) == 0 ) sf = 0.1;
	    else                                    sf = 1; // meter

            grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, lbounds, level*sf, dlevel1*sf, dlevel2*sf);
	  }

	break;
      }
    case ZAXIS_REFERENCE:
      {
        if ( !gcinit ) GRIB_CHECK(my_grib_set_long(gh, "genVertHeightCoords", 1), 0);

        if ( editionNumber <= 1 )
          ; // not available
        else
          {
            if ( lbounds )
              {
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", grib_ltype);
                gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", grib_ltype);
                GRIB_CHECK(my_grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
            else
              {
                grib2DefLevel(gh, gcinit, GRIB2_LTYPE_REFERENCE, GRIB2_LTYPE_REFERENCE, lbounds, level, dlevel1, dlevel2);
              }

            int number = zaxisInqNumber(zaxisID);
            unsigned char uuid[CDI_UUID_SIZE];
            GRIB_CHECK(my_grib_set_long(gh, "NV", 6), 0);
            GRIB_CHECK(my_grib_set_long(gh, "nlev", zaxisInqNlevRef(zaxisID)), 0);
            GRIB_CHECK(my_grib_set_long(gh, "numberOfVGridUsed", number), 0);
            size_t len = CDI_UUID_SIZE;
            zaxisInqUUID(zaxisID, uuid);
            if ( grib_set_bytes(gh, "uuidOfVGrid", uuid, &len) != 0 ) Warning("Can't write UUID!");
          }

        break;
      }
    case ZAXIS_GENERIC:
      {
	if ( editionNumber <= 1 )
          {
            grib1DefLevel(gh, gcinit, ltype, lbounds, level, dlevel1, dlevel2);
          }
        else
          {
            grib2DefLevel(gh, gcinit, ltype, ltype, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}


int gribapiGetScanningMode(grib_handle *gh)
{
  long iScansNegatively;
  long jScansPositively;
  long jPointsAreConsecutive;

  GRIB_CHECK(grib_get_long(gh, "iScansNegatively", &iScansNegatively), 0);
  GRIB_CHECK(grib_get_long(gh, "jScansPositively", &jScansPositively), 0);
  GRIB_CHECK(grib_get_long(gh, "jPointsAreConsecutive", &jPointsAreConsecutive), 0);
  int scanningMode
    = 128*(bool)iScansNegatively
    + 64 *(bool)jScansPositively
    + 32 *(bool)jPointsAreConsecutive;
  if (cdiDebugExt>=30)
    printf("gribapiGetScanningMode(): Scanning mode = %02d (%1d%1d%1d)*32; \n",\
            scanningMode,(int)jPointsAreConsecutive,(int)jScansPositively,(int)iScansNegatively);

 return scanningMode;
}


void gribapiSetScanningMode(grib_handle *gh, int scanningMode)
{
   // 127: reserved for testing; generated test data will be in 64 scanning mode
  //if (scanningMode== 127)  scanningMode = 64;

  long iScansNegatively      = (scanningMode & 128)/128;
  long jScansPositively      = (scanningMode & 64)/64;
  long jPointsAreConsecutive = (scanningMode & 32)/32;

  if (cdiDebugExt>=30)
  {
    long paramId, levelTypeId, levelId, uvRelativeToGrid;
    GRIB_CHECK(grib_get_long(gh, "uvRelativeToGrid", &uvRelativeToGrid), 0);
    GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &paramId), 0);
    GRIB_CHECK(grib_get_long(gh, "indicatorOfTypeOfLevel", &levelTypeId), 0);
    GRIB_CHECK(grib_get_long(gh, "level", &levelId), 0);
    printf("gribapiSetScanningMode(): (param,ltype,level) = (%3d,%3d,%4d); Scanning mode = %02d (%1d%1d%1d)*32;  uvRelativeToGrid = %02d\n",\
            (int)paramId, (int)levelTypeId, (int)levelId,
            scanningMode,(int)jPointsAreConsecutive,(int)jScansPositively,(int)iScansNegatively,
            (int)uvRelativeToGrid);
  }

  GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", iScansNegatively), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", jScansPositively), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jPointsAreConsecutive", jPointsAreConsecutive), 0);
}


static void gribapiSetUvRelativeToGrid(grib_handle *gh, int mode)
{
  long uvRelativeToGridMode = mode;
  long uvRelativeToGridModeOld;

  GRIB_CHECK(grib_get_long(gh, "uvRelativeToGrid", &uvRelativeToGridModeOld), 0);

  if (cdiDebugExt>=30)
    printf("gribapiSetUvRelativeToGrid():  uvRelativeToGrid: %02d (old) => %02d (new); \n",(int)uvRelativeToGridModeOld,(int)uvRelativeToGridMode);

  GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGridMode), 0);
}


  /*
    TABLE 8. SCANNING MODE FLAG

    (GDS Octet 28)
    BIT     VALUE     MEANING
    1       0       Points scan in +i direction
            1       Points scan in -i direction
    2       0       Points scan in -j direction
            1       Points scan in +j direction
    3       0       Adjacent points in i direction are consecutive
                      (FORTRAN: (I,J))
            1       Adjacent points in j direction are consecutive
                    (FORTRAN: (J,I))

    => Scanning Mode     0 0 0 0 0 0 0 0  (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & North->South)
    => Scanning Mode     0 1 0 0 0 0 0 0  (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & South->North )
    => Scanning Mode     1 1 0 0 0 0 0 0  (96 dec)  +i, +j; j direction consecutive (column-major order South->North & West->East )

    NOTE:  South->North  - As if you would plot the data as image on the screen
                           where [0,0] of the data is the top-left pixel.

                           grib2ppm LAMH_D11_201302150000_00000_oro | display ppm:-
                           ImageMagick (display): [0,0] of an image belongs to the top-left pixel
    [DEFAULT] : 64 dec

    iScansNegatively = 0;
    jScansPositively = 1;
    jPointsAreConsecutive = 0;    => Scanning Mode 64

    cdo selindexbox,1,726,100,550 LAMH_D11_201302150000_00000_oro LAMH_D11_201302150000_00000_oro_cropped
    grib2ppm LAMH_D11_201302150000_00000_oro_cropped | /usr/bin/display ppm:- &
    # ^^^ this image will be missing the souther parts of data

    grib2ppm LAMH_D11_201302150000_00000_oro | /usr/bin/display ppm:- &
    # ^ full domain data
  */

#ifdef HIRLAM_EXTENSIONS
static void
verticallyFlipGridDefinitionWhenScanningModeChanged(grib_handle *gh, double yfirst, double ylast, double yinc )
{
  /*
  Nj = 550;
  latitudeOfFirstGridPointInDegrees = -30.8;
  latitudeOfLastGridPointInDegrees = 24.1;
  iScansNegatively = 0;
  jScansPositively = 0;
  jPointsAreConsecutive = 0;
  jDirectionIncrementInDegrees = 0.1;

  When switching from scanning mode 0 <=> 64
  yfirst = -30.8 + (550-1)*0.1

  yfirst = yfirst + (ysize-1) * yinc
  yinc   = -1.0*yinc

  */


  //long jDim=0;
  //GRIB_CHECK(grib_get_long(gh, "Nj", &jDim), 0);

  double latitudeOfFirstGridPointInDegrees;
  double latitudeOfLastGridPointInDegrees;
  double jDirectionIncrementInDegrees;

  //GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &latitudeOfFirstGridPointInDegrees), 0);  // yfirst
  //GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees", &latitudeOfLastGridPointInDegrees), 0);    // ylast
  //GRIB_CHECK(grib_get_double(gh, "jDirectionIncrementInDegrees", &jDirectionIncrementInDegrees), 0);  // yinc

  if (cdiDebugExt>=10)
  {
      Message(" BEFORE: yfirst = %f; ylast = %f; yinc = %f; ", yfirst,ylast, yinc);
  }

  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", ylast), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfLastGridPointInDegrees", yfirst), 0);
  //yinc *= -1.0; // don't set yinc here ...
  //GRIB_CHECK(my_grib_set_double(gh, "jDirectionIncrementInDegrees", yinc), 0);

  if (cdiDebugExt>=10)
  {
    GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &latitudeOfFirstGridPointInDegrees), 0);  // yfirst
    GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees", &latitudeOfLastGridPointInDegrees), 0);    // ylast
    GRIB_CHECK(grib_get_double(gh, "jDirectionIncrementInDegrees", &jDirectionIncrementInDegrees), 0);  // yinc
    Message("CHANGED INTO:  yfirst = %f, ylast = %f, yinc = %f",latitudeOfFirstGridPointInDegrees,latitudeOfLastGridPointInDegrees, jDirectionIncrementInDegrees);
  }
}

static void
convertDataScanningMode(int scanModeIN, int scanModeOUT, double *data,
                        int gridsize, int iDim, int jDim)
{
  int i,j;
  int idxIN, idxOUT;

   // 127: reserved for testing; it will generate test data in 64 scanning mode
  if (scanModeOUT== 127)  // fill with testdata ...
  {
      scanModeOUT = 64;
      if (cdiDebugExt>=30) printf("convertDataScanningMode(): Generating test data in 64 scanning mode..\n");
      for (j=0; j<jDim; j++)
      {
        int jXiDim = j*iDim;
        for (i=0; i<iDim; i++)
        {
          idxIN = i + jXiDim;
          data[idxIN] = (double) (100.0*j +i);
        }
      }
  }

  if ( (iDim*jDim)!= gridsize)
  {
    if (cdiDebugExt>=30) printf("convertDataScanningMode(): ERROR: (iDim*jDim)!= gridsize;  (%d * %d) != %d\n", iDim,jDim, gridsize);
    return;
  }
  if (cdiDebugExt>=30) printf("convertDataScanningMode(): scanModeIN=%02d => scanModeOUT=%02d ; where: (iDim * jDim == gridsize)  (%d*%d == %d)\n",scanModeIN, scanModeOUT, iDim,jDim, gridsize);

  if (cdiDebugExt>=100)
  {
      printf("convertDataScanningMode(): data IN:\n");
      for (j=0; j<jDim; j++)
      {
        int jXiDim = j*iDim;
        for (i=0; i<iDim; i++)
        {
          idxIN = i + jXiDim;
          printf("%03.0f, ",data[idxIN]);
        }
        printf("\n");
      }
  }

  if (scanModeIN==scanModeOUT)
  {
    if (cdiDebugExt>=30) printf("convertDataScanningMode(): INFO: Nothing to do;  scanModeIN==scanModeOUT..\n");
    return;
  }

  if (0)
  {
      return;
      if (scanModeOUT==00)
      {
          if (cdiDebugExt>0) printf("convertDataScanningMode(): Leave data unchaged BUT set scanModeOUT=00.\n");
          // CHECK:  Looks like that GRIB-API provide (no matter what) data in the scannning mode 00, even it is store in the gribfile as 64 !!
          return;
      }
  }
  double *dataCopy = NULL;
  dataCopy = (double *) malloc(gridsize*sizeof(double));

  memcpy((void*)dataCopy,(void*) data, gridsize*sizeof(double));

  if (scanModeIN==64)           // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
  {                             // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
                                // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
      if (scanModeOUT==00)
      // CHECK:  Looks like that GRIB-API provide (no matter what) data in the scannning mode 00, even it is store in the gribfile as 64 !!
#define VERTICAL_FLIP
#ifdef VERTICAL_FLIP
      { // flip the data vertically ..
        idxIN= 0; idxOUT= (jDim-1)*iDim;
        if (cdiDebugExt>=30) printf("convertDataScanningMode():  copying rows nr. (%04d : %04d)\n",0,jDim-1);
        for (j=0; j<jDim; j++)
        {
          memcpy((void*)&data[idxOUT], (void*)&dataCopy[idxIN], iDim*sizeof(double));
          idxIN  += iDim; idxOUT -= iDim;
        }
      } // end if (scanModeOUT==00)*/
#endif
#ifdef HORIZONTAL_FLIP
      { // flip data horizontally ...
        if (1)
        {
            if (cdiDebugExt>=30) printf("convertDataScanningMode():  copying columns nr. (%04d : %04d);\n", 0, iDim-1);
            for (i=0; i<iDim; i++)
            {
              for (j=0; j<jDim; j++)
              {
                int jXiDim = j*iDim;
                idxIN  = i           + jXiDim;
                //data[idxIN] = (double) (100.0*j +i);  // just some testdata ..
                idxOUT = iDim - i -1 + jXiDim;
                //printf("[%03d=>%03d] = %f;",idxIN,idxOUT,dataCopy[idxIN]);
                data[idxOUT] =  dataCopy[idxIN];
              }
            }
        }
      } // end if (scanModeOUT==00)
#endif

      if (scanModeOUT==96)
      { // transpose the data
        if (cdiDebugExt>=30) printf("convertDataScanningMode():  transpose data rows=>columns nr. (%04d : %04d) => (%04d : %04d);\n", 0, iDim-1, 0, jDim-1);
        for (j=0; j<jDim; j++)
        {
          int jXiDim = j*iDim;
          for (i=0; i<iDim; i++)
          {
            idxIN  = i + jXiDim;
            idxOUT = j + i*jDim;
            //printf("[%03d=>%03d] = %f;",idxIN,idxOUT,dataCopy[idxIN]);
            data[idxOUT] =  dataCopy[idxIN];
          }
          //printf(".\n");
        }
      } // end if (scanModeOUT==96)
  } // end if (scanModeIN==64)

  if (scanModeIN==00)           // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
  {                             // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
                               // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
    if (scanModeOUT==64)
      { // flip the data vertically ..
        idxIN= 0; idxOUT= (jDim-1)*iDim;
        for (j=0; j<jDim; j++)
        {
          if (cdiDebugExt>=25) printf("convertDataScanningMode():  copying row nr. %04d; [idxIN=%08d] => [idxOUT=%08d]\n",j, idxIN, idxOUT);
          memcpy((void*)&data[idxOUT], (void*)&dataCopy[idxIN], iDim*sizeof(double));
          idxIN  += iDim; idxOUT -= iDim;
        }
      } // end if (scanModeOUT==64)

      if (scanModeOUT==96)
      { // transpose the data
        int jInv;
        for (j=0; j<jDim; j++)
        {
          if (cdiDebugExt>=30) printf("convertDataScanningMode():  processing row nr. %04d;\n", j);
          jInv = (jDim-1) -j;
          for (i=0; i<iDim; i++)
            data[j + i*jDim] =  dataCopy[i + jInv*iDim];  // source data has -j
        }
      } // end if (scanModeOUT==96)
  } // end if (scanModeIN==00)

  if (scanModeIN==96)           // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
  {                             // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
                                // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
    if (scanModeOUT==64)
      { // transpose the data
        for (j=0; j<jDim; j++)
        {
          if (cdiDebugExt>=30) printf("convertDataScanningMode():  processing row nr. %04d;\n", j);
          int jXiDim = j*iDim;
          for (i=0; i<iDim; i++)
            //data[j + i*jDim] =  dataCopy[i + j*iDim];
            data[i + jXiDim] =  dataCopy[j + i*jDim];
        }
      } // end if (scanModeOUT==64)

      if (scanModeOUT==00)
      { // transpose the data
        idxIN= 0; idxOUT= 0;
        int jInv;
        for (j=0; j<jDim; j++)
        {
          if (cdiDebugExt>=30) printf("convertDataScanningMode():  processing row nr. %04d;\n", j);
          jInv = (jDim-1) -j;
          int jXiDim = j*iDim;
          for (i=0; i<iDim; i++)
            //data[jInv + iXjDim] =  dataCopy[i + jXiDim];  // target data has -j
            data[i + jXiDim] =  dataCopy[jInv + i*jDim];  // target data has -j
        }
      } // end if (scanModeOUT==00)
  } // end if (scanModeIN==96)

  if (cdiDebugExt>=100)
  {
      printf("convertDataScanningMode(): data OUT (new scanning mode):\n");
      for (j=0; j<jDim; j++)
      {
        int jXiDim = j*iDim;
        for (i=0; i<iDim; i++)
        {
          idxIN = i + jXiDim;
          printf("%03.0f, ",data[idxIN]);
        }
        printf("\n");
      }
  }

  free(dataCopy); return;
}
#endif //HIRLAM_EXTENSIONS

static
void gribapiSetExtMode(grib_handle *gh, int gridID, long datasize, const double *data)
{
  /*
  Nj = 550;
  latitudeOfFirstGridPointInDegrees = -30.8;
  latitudeOfLastGridPointInDegrees = 24.1;
  iScansNegatively = 0;
  jScansPositively = 0;
  jPointsAreConsecutive = 0;
  jDirectionIncrementInDegrees = 0.1; */
#ifndef HIRLAM_EXTENSIONS
  (void)data;
  (void)datasize;
#endif
  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_GENERIC || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GAUSSIAN_REDUCED || gridtype == GRID_PROJECTION )
    {
#ifdef HIRLAM_EXTENSIONS
      int scanModeIN = gridInqScanningMode(gridID);

      if (cdiDebugExt>=100)
        {
          int gridsize = gridInqSize(gridID);
          Message("(scanModeIN=%d; gridsize=%d", scanModeIN, gridsize);
        }

      if ( cdiGribDataScanningMode.active )   // allowed modes: <0, 64, 96>; Default is 64
        {
          int iDim = gridInqXsize(gridID);
          int jDim = gridInqYsize(gridID);

          double yfirst = gridInqYval(gridID,      0);
          double ylast  = gridInqYval(gridID, jDim-1);
          double yinc   = gridInqYinc(gridID);

          int scanModeOUT = cdiGribDataScanningMode.value;
          convertDataScanningMode(scanModeIN, scanModeOUT, (double*)data, datasize, iDim, jDim);
          // This will overrule the old scanning mode of the given grid
          if (cdiDebugExt>=10) Message("Set GribDataScanningMode (%d) => (%d)", scanModeIN, cdiGribDataScanningMode.value);
          gribapiSetScanningMode(gh, cdiGribDataScanningMode.value);

          if (((scanModeIN==00) && (cdiGribDataScanningMode.value==64)) ||
              ((scanModeIN==64) && (cdiGribDataScanningMode.value==00)) )
            verticallyFlipGridDefinitionWhenScanningModeChanged(gh, yfirst, ylast, yinc);
        }
      else
        {
          if (cdiDebugExt>=100) Message("Set GribDataScanningMode => (%d) based on used grid", scanModeIN);
          gribapiSetScanningMode(gh, scanModeIN);
        }
#endif

      if ( cdiGribChangeModeUvRelativeToGrid.active )
        {
          // this will overrule/change the UvRelativeToGrid flag;
          // typically when the wind is rotated with respect to north pole
          if (cdiDebugExt>=100) Message("Set ModeUvRelativeToGrid =>%d ( note grid has: %d)", cdiGribChangeModeUvRelativeToGrid.mode, gridInqUvRelativeToGrid(gridID));
          GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", (long) cdiGribChangeModeUvRelativeToGrid.mode), 0);
        }
      else
        {
          if (cdiDebugExt>=100) Message("Set ModeUvRelativeToGrid =>%d based on used grid", gridInqUvRelativeToGrid(gridID));
          gribapiSetUvRelativeToGrid(gh, gridInqUvRelativeToGrid(gridID));
        }
    }
}

/* #define GRIBAPIENCODETEST 1 */

size_t gribapiEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg,
		     long datasize, const double *data, int nmiss, void **gribbuffer, size_t *gribbuffersize,
		     int comptype, void *gribContainer)
{
  size_t recsize = 0;
  void *dummy = NULL;
  bool lieee = false;
  /*  int ensID, ensCount, forecast_type; *//* Ensemble Data */
  int typeOfGeneratingProcess;
  int productDefinitionTemplate;
  long bitsPerValue;
  long editionNumber = 2;
  char name[256];
  char stdname[256];

  // extern unsigned char _grib_template_GRIB2[];

  int param    = vlistInqVarParam(vlistID, varID);
  int datatype = vlistInqVarDatatype(vlistID, varID);
  typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID);
  productDefinitionTemplate = vlistInqVarProductDefinitionTemplate(vlistID, varID);

  vlistInqVarName(vlistID, varID, name);
  vlistInqVarStdname(vlistID, varID, stdname);

#if defined(GRIBAPIENCODETEST)
  grib_handle *gh = (grib_handle *) gribHandleNew(editionNumber);
#else
  gribContainer_t *gc = (gribContainer_t *) gribContainer;
  assert(gc != NULL);
  grib_handle *gh = (struct grib_handle *)gc->gribHandle;
#endif

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber == 2 )
    {
      if ( typeOfGeneratingProcess == -1 ) typeOfGeneratingProcess = 0;
      if ( ! gc->init ) GRIB_CHECK(my_grib_set_long(gh, "typeOfGeneratingProcess", typeOfGeneratingProcess), 0);
    }

  /*
  if( vlistInqVarEnsemble( vlistID,  varID, &ensID, &ensCount, &forecast_type ) )
    {
      GRIB_CHECK(my_grib_set_long(gh, "typeOfEnsembleForecast", forecast_type ), 0);
      GRIB_CHECK(my_grib_set_long(gh, "numberOfForecastsInEnsemble", ensCount ), 0);
      GRIB_CHECK(my_grib_set_long(gh, "perturbationNumber", ensID ), 0);
    }
  */

  gribapiDefTime((int)editionNumber, productDefinitionTemplate, typeOfGeneratingProcess,
                 gh, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID), gc->init);

  if ( ! gc->init ) gribapiDefInstitut(gh, vlistID, varID);
  if ( ! gc->init ) gribapiDefModel(gh, vlistID, varID);

  if ( ! gc->init ) gribapiDefParam((int)editionNumber, gh, param, name, stdname);

  if ( editionNumber == 2 && (datatype == CDI_DATATYPE_FLT32 || datatype == CDI_DATATYPE_FLT64) ) lieee = true;

  /* bitsPerValue have to be defined before call to DefGrid (complex packing) */
  //  if ( lieee == false )
    {
      bitsPerValue = grbBitsPerValue(datatype);
      GRIB_CHECK(my_grib_set_long(gh, "bitsPerValue", bitsPerValue), 0);
    }

  gribapiDefGrid((int)editionNumber, gh, gridID, comptype, lieee, datatype, nmiss, gc->init);

  gribapiDefLevel((int)editionNumber, gh, zaxisID, levelID, gc->init, productDefinitionTemplate);

  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  //if (!gc->init)
  {
    int ret = 0;

    /* NOTE: Optional key/value pairs: Note that we do not distinguish
     *       between tiles here! */

    for ( int i=0; i<vlistptr->vars[varID].opt_grib_nentries; i++ )
      {
        if ( vlistptr->vars[varID].opt_grib_kvpair[i].update )
          {
            //DR: Fix for multi-level fields (otherwise only the 1st level is correct)
            if ( zaxisInqSize(zaxisID)==(levelID+1) )
              vlistptr->vars[varID].opt_grib_kvpair[i].update = false;

            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_double)
              {
                if ( CDI_Debug )
                  Message("key \"%s\"  :   double value = %g\n",
                          vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                my_grib_set_double(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                   vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                GRIB_CHECK(ret, 0);
                }
            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_int)
              {
                if ( CDI_Debug )
                  Message("key \"%s\"  :   integer value = %d\n",
                          vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                my_grib_set_long(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                 (long) vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                GRIB_CHECK(ret, 0);
              }
          }
      }
  }

  if ( nmiss > 0 )
    {
      GRIB_CHECK(my_grib_set_long(gh, "bitmapPresent", 1), 0);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", vlistInqVarMissval(vlistID, varID)), 0);
    }

  gribapiSetExtMode(gh, gridID, datasize, data);

  GRIB_CHECK(grib_set_double_array(gh, "values", data, (size_t)datasize), 0);

  /* get the size of coded message  */
  GRIB_CHECK(grib_get_message(gh, (const void **)&dummy, &recsize), 0);
  recsize += 512; /* add some space for possible filling */
  *gribbuffersize = recsize;
  *gribbuffer = Malloc(*gribbuffersize);

  /* get a copy of the coded message */
  GRIB_CHECK(grib_get_message_copy(gh, *gribbuffer, &recsize), 0);

#if defined(GRIBAPIENCODETEST)
  gribHandleDelete(gh);
#endif

  gc->init = true;

  return recsize;
}


void gribapiChangeParameterIdentification(grib_handle *gh, int code, int ltype, int lev)
{
  long  indicatorOfParameter,  indicatorOfTypeOfLevel,  level; //  timeRangeIndicator: could be included later
  indicatorOfParameter = code;
  indicatorOfTypeOfLevel = ltype;
  level = lev;

  if (indicatorOfParameter!=-1) GRIB_CHECK(my_grib_set_long(gh, "indicatorOfParameter", indicatorOfParameter), 0);
  if (indicatorOfTypeOfLevel!=-1) GRIB_CHECK(my_grib_set_long(gh, "indicatorOfTypeOfLevel", indicatorOfTypeOfLevel), 0);
  if (level!=-1) GRIB_CHECK(my_grib_set_long(gh, "level", level), 0);
}

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
