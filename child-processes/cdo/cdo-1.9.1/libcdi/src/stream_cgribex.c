#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"
#include "stream_cgribex.h"

#if  defined  (HAVE_LIBCGRIBEX)
#  include "cgribex.h"
#endif

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
} compvar_t;


#if  defined  (HAVE_LIBCGRIBEX)
static
int cgribexGetGridType(int *isec2)
{
  int gridtype = GRID_GENERIC;

  switch (ISEC2_GridType)
    {
    case  GRIB1_GTYPE_LATLON:     { gridtype = GRID_LONLAT;     break; }
    case  GRIB1_GTYPE_LATLON_ROT: { gridtype = GRID_PROJECTION; break; }
    case  GRIB1_GTYPE_LCC:        { gridtype = GRID_LCC;        break; }
    case  GRIB1_GTYPE_GAUSSIAN:   { gridtype = ISEC2_Reduced ? GRID_GAUSSIAN_REDUCED : GRID_GAUSSIAN; break; }
    case  GRIB1_GTYPE_SPECTRAL:   { gridtype = GRID_SPECTRAL;   break; }
    case  GRIB1_GTYPE_GME:        { gridtype = GRID_GME;        break; }
    }

  return gridtype;
}

static
bool cgribexGetIsRotated(int *isec2)
{
  return (ISEC2_GridType == GRIB1_GTYPE_LATLON_ROT);
}

static
int cgribexGetZaxisHasBounds(int grb_ltype)
{
  int lbounds = 0;

  switch (grb_ltype)
    {
    case GRIB1_LTYPE_SIGMA_LAYER:
    case GRIB1_LTYPE_HYBRID_LAYER:
    case GRIB1_LTYPE_LANDDEPTH_LAYER:
      {
	lbounds = 1;
	break;
      }
    }

  return lbounds;
}

static
int cgribexGetTimeUnit(int *isec1)
{
  int timeunit = TUNIT_HOUR;
  static bool lprint = true;

  switch ( ISEC1_TimeUnit )
    {
    case ISEC1_TABLE4_MINUTE:    timeunit = TUNIT_MINUTE;    break;
    case ISEC1_TABLE4_QUARTER:   timeunit = TUNIT_QUARTER;   break;
    case ISEC1_TABLE4_30MINUTES: timeunit = TUNIT_30MINUTES; break;
    case ISEC1_TABLE4_HOUR:      timeunit = TUNIT_HOUR;      break;
    case ISEC1_TABLE4_3HOURS:    timeunit = TUNIT_3HOURS;    break;
    case ISEC1_TABLE4_6HOURS:    timeunit = TUNIT_6HOURS;    break;
    case ISEC1_TABLE4_12HOURS:   timeunit = TUNIT_12HOURS;   break;
    case ISEC1_TABLE4_DAY:       timeunit = TUNIT_DAY;       break;
    default:
      if ( lprint )
	{
	  Message("GRIB time unit %d unsupported!", ISEC1_TimeUnit);
	  lprint = false;
	}
      break;
    }

  return timeunit;
}

static
bool cgribexTimeIsFC(int *isec1)
{
  bool isFC = (ISEC1_TimeRange == 10 && ISEC1_TimePeriod1 == 0 && ISEC1_TimePeriod2 == 0) ? false : true;

  return isFC;
}

static
int cgribexGetTsteptype(int timerange)
{
  static bool lprint = true;

  int tsteptype = TSTEP_INSTANT;
  switch ( timerange )
    {
    case  0:  tsteptype = TSTEP_INSTANT;  break;
    case  1:  tsteptype = TSTEP_INSTANT2; break;
    case  2:  tsteptype = TSTEP_RANGE;    break;
    case  3:  tsteptype = TSTEP_AVG;      break;
    case  4:  tsteptype = TSTEP_ACCUM;    break;
    case  5:  tsteptype = TSTEP_DIFF;     break;
    case 10:  tsteptype = TSTEP_INSTANT3; break;
    default:
      if ( lprint )
	{
	  Message("Time range indicator %d unsupported, set to 0!", timerange);
	  lprint = false;
	}
      break;
    }

  return tsteptype;
}

static
void cgribexGetGrid(stream_t *streamptr, int *isec2, int *isec4, grid_t *grid, int iret)
{
  bool compyinc = true;
  int gridtype = cgribexGetGridType(isec2);
  int projtype = (gridtype == GRID_PROJECTION && cgribexGetIsRotated(isec2)) ? CDI_PROJ_RLL : CDI_UNDEFID;
  if ( gridtype == GRID_LCC )
    {
      gridtype = GRID_PROJECTION;
      projtype = CDI_PROJ_LCC;
    }

  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED && iret != -801 )
    {
      int nlon = 0;
      for ( int ilat = 0; ilat < ISEC2_NumLat; ++ilat )
        if ( ISEC2_RowLon(ilat) > nlon ) nlon = ISEC2_RowLon(ilat);
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = nlon;
      ISEC4_NumValues = nlon*ISEC2_NumLat;
      compyinc = false;
    }

  grid_init(grid);
  cdiGridTypeInit(grid, gridtype, 0);

  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_RLL )
    {
      bool ijDirectionIncrementGiven = gribbyte_get_bit(ISEC2_ResFlag, 1);
      bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);

      if ( uvRelativeToGrid ) grid->uvRelativeToGrid = 1;
      if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
        Error("numberOfPoints (%d) and gridSize (%d) differ!", ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);

      grid->size   = ISEC4_NumValues;
      grid->x.size = ISEC2_NumLon;
      grid->y.size = ISEC2_NumLat;
      if ( gridtype == GRID_GAUSSIAN ) grid->np = ISEC2_NumPar;
      grid->x.inc  = 0;
      grid->y.inc  = 0;
      grid->x.flag = 0;
      /* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
      {
        if ( grid->x.size > 1 )
          {
            bool recompinc = true;

            if ( ISEC2_LastLon < ISEC2_FirstLon && ISEC2_LastLon < 0 ) ISEC2_LastLon += 360000;

            if ( ijDirectionIncrementGiven && ISEC2_LonIncr > 0 )
              {
                if ( abs(ISEC2_LastLon - (ISEC2_FirstLon+ISEC2_LonIncr*(grid->x.size-1))) <= 2 )
                  {
                    recompinc = false;
                    grid->x.inc = ISEC2_LonIncr * 0.001;
                  }
              }

            /* recompute xinc if necessary */
            if ( recompinc ) grid->x.inc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->x.size-1);

            /* correct xinc if necessary */
            if ( ISEC2_FirstLon == 0 && ISEC2_LastLon > 354000 && ISEC2_LastLon < 360000 )
              {
                double xinc = 360. / grid->x.size;
                if ( fabs(grid->x.inc-xinc) > 0.0 )
                  {
                    grid->x.inc = xinc;
                    if ( CDI_Debug ) Message("set xinc to %g", grid->x.inc);
                  }
              }
          }
        grid->x.first = ISEC2_FirstLon * 0.001;
        grid->x.last  = ISEC2_LastLon  * 0.001;
        grid->x.flag  = 2;
      }
      grid->y.flag = 0;
      /* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
      {
        if ( grid->y.size > 1 && compyinc )
          {
            bool recompinc = true;
            if ( ijDirectionIncrementGiven && ISEC2_LatIncr > 0 )
              {
                if ( abs(ISEC2_LastLat - (ISEC2_FirstLat+ISEC2_LatIncr*(grid->y.size-1))) <= 2 )
                  {
                    recompinc = false;
                    grid->y.inc = ISEC2_LatIncr * 0.001;
                  }
              }

            /* recompute yinc if necessary */
            if ( recompinc ) grid->y.inc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->y.size - 1);
          }
        grid->y.first = ISEC2_FirstLat * 0.001;
        grid->y.last  = ISEC2_LastLat  * 0.001;
        grid->y.flag  = 2;
      }
    }
  else if ( gridtype == GRID_GAUSSIAN_REDUCED )
    {
      bool ijDirectionIncrementGiven = gribbyte_get_bit(ISEC2_ResFlag, 1);
      bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);
      if ( uvRelativeToGrid ) grid->uvRelativeToGrid = 1;
      grid->np      = ISEC2_NumPar;
      grid->size    = ISEC4_NumValues;
      grid->rowlon  = ISEC2_RowLonPtr;
      grid->nrowlon = ISEC2_NumLat;
      grid->y.size  = ISEC2_NumLat;
      grid->x.inc   = 0;
      grid->y.inc   = 0;
      grid->x.flag  = 0;
      /* if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 ) */
      {
        if ( grid->x.size > 1 )
          {
            if ( ISEC2_LastLon < ISEC2_FirstLon && ISEC2_LastLon < 0 ) ISEC2_LastLon += 360000;

            if ( ijDirectionIncrementGiven && ISEC2_LonIncr > 0 )
              grid->x.inc = ISEC2_LonIncr * 0.001;
            else
              grid->x.inc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->x.size - 1);
          }
        grid->x.first = ISEC2_FirstLon * 0.001;
        grid->x.last  = ISEC2_LastLon  * 0.001;
        grid->x.flag  = 2;
      }
      grid->y.flag = 0;
      /* if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 ) */
      {
        if ( grid->y.size > 1 )
          {
            if ( ijDirectionIncrementGiven && ISEC2_LatIncr > 0 )
              grid->y.inc = ISEC2_LatIncr * 0.001;
            else
              grid->y.inc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->y.size - 1);
          }
        grid->y.first = ISEC2_FirstLat * 0.001;
        grid->y.last  = ISEC2_LastLat  * 0.001;
        grid->y.flag  = 2;
      }
    }
  else if ( projtype == CDI_PROJ_LCC )
    {
      bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);
      if ( uvRelativeToGrid ) grid->uvRelativeToGrid = 1;

      if ( ISEC4_NumValues != ISEC2_NumLon*ISEC2_NumLat )
        Error("numberOfPoints (%d) and gridSize (%d) differ!", ISEC4_NumValues, ISEC2_NumLon*ISEC2_NumLat);

      grid->size   = ISEC4_NumValues;
      grid->x.size = ISEC2_NumLon;
      grid->y.size = ISEC2_NumLat;

      grid->x.first = 0;
      grid->x.last  = 0;
      grid->x.inc   = ISEC2_Lambert_dx;
      grid->y.first = 0;
      grid->y.last  = 0;
      grid->y.inc   = ISEC2_Lambert_dy;
      grid->x.flag  = 2;
      grid->y.flag  = 2;
    }
  else if ( gridtype == GRID_SPECTRAL )
    {
      grid->size  = ISEC4_NumValues;
      grid->trunc = ISEC2_PentaJ;
      if ( ISEC2_RepMode == 2 )
        grid->lcomplex = 1;
      else
        grid->lcomplex = 0;
    }
  else if ( gridtype == GRID_GME )
    {
      grid->size  = ISEC4_NumValues;
      grid->gme.nd  = ISEC2_GME_ND;
      grid->gme.ni  = ISEC2_GME_NI;
      grid->gme.ni2 = ISEC2_GME_NI2;
      grid->gme.ni3 = ISEC2_GME_NI3;
    }
  else if ( gridtype == GRID_GENERIC )
    {
      grid->size  = ISEC4_NumValues;
      grid->x.size = 0;
      grid->y.size = 0;
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }

  grid->type = gridtype;
  grid->projtype = projtype;
}

static
void cgribexGetLevel(int *isec1, int *leveltype, int *level1, int *level2)
{
  *leveltype = ISEC1_LevelType;
  *level1 = ISEC1_Level1;
  *level2 = ISEC1_Level2;
  if ( *leveltype == GRIB1_LTYPE_ISOBARIC ) *level1 *= 100;
  else if ( *leveltype == GRIB1_LTYPE_99 || *leveltype == GRIB1_LTYPE_ISOBARIC_PA ) *leveltype = GRIB1_LTYPE_ISOBARIC;
}

static
void cgribexAddRecord(stream_t * streamptr, int param, int *isec1, int *isec2, double *fsec2, double *fsec3,
		      int *isec4, size_t recsize, off_t position, int datatype, int comptype, int lmv, int iret)
{
  int varID;
  int levelID = 0;

  int vlistID = streamptr->vlistID;
  int tsID    = streamptr->curTsID;
  int recID   = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  int tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);
  int numavg    = ISEC1_AvgNum;

  int leveltype, level1, level2;
  cgribexGetLevel(isec1, &leveltype, &level1, &level2);

  /* fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype); */

  record->size      = recsize;
  record->position  = position;
  record->param     = param;
  record->ilevel    = level1;
  record->ilevel2   = level2;
  record->ltype     = leveltype;
  record->tsteptype = (short)tsteptype;

  grid_t *gridptr = (grid_t*) Malloc(sizeof(*gridptr));
  cgribexGetGrid(streamptr, isec2, isec4, gridptr, iret);

  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, gridptr, 0);
  int gridID = gridAdded.Id;
  if ( gridAdded.isNew )
    {
      if ( gridptr->nrowlon )
        {
          size_t nrowlon = (size_t) gridptr->nrowlon;
          int *rowlon = gridptr->rowlon;
          gridptr->rowlon = (int*) Malloc(nrowlon * sizeof(int));
          memcpy(gridptr->rowlon, rowlon, nrowlon * sizeof(int));
        }
      else if ( gridptr->projtype == CDI_PROJ_RLL )
        {
          double xpole =   ISEC2_LonSP*0.001 - 180;
          double ypole = - ISEC2_LatSP*0.001;
          double angle = - FSEC2_RotAngle;
          gridDefParamRLL(gridID, xpole, ypole, angle);
        }
      else if ( gridptr->projtype == CDI_PROJ_LCC )
        {
          double a = 6367470., rf = 0;
          bool earthIsOblate = gribbyte_get_bit(ISEC2_ResFlag, 2);
          if ( earthIsOblate ) { a = 6378160.; rf = 297.0; }
          double xval_0 = ISEC2_FirstLon * 0.001;
          double yval_0 = ISEC2_FirstLat * 0.001;
          double lon_0  = ISEC2_Lambert_Lov * 0.001;
          double lat_1  = ISEC2_Lambert_LatS1 * 0.001;
          double lat_2  = ISEC2_Lambert_LatS2 * 0.001;
          bool lsouth = gribbyte_get_bit(ISEC2_Lambert_ProjFlag, 1);
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
    }
  else
    Free(gridptr);

  int zaxistype = grib1ltypeToZaxisType(leveltype);

  if ( zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
    {
      size_t vctsize = (size_t)ISEC2_NumVCP;
      double *vctptr = &fsec2[10];

      varDefVCT(vctsize, vctptr);
    }

  int lbounds = cgribexGetZaxisHasBounds(leveltype);

  if ( datatype > 32 ) datatype = CDI_DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = CDI_DATATYPE_PACK;

  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, 0, 0,
	       datatype, &varID, &levelID, tsteptype, numavg, leveltype, -1,
               NULL, NULL, NULL, NULL, NULL, NULL);

  record->varID   = (short)varID;
  record->levelID = (short)levelID;

  varDefCompType(varID, comptype);

  if ( ISEC1_LocalFLag )
    {
      if      ( ISEC1_CenterID == 78  && isec1[36] == 253 ) // DWD local extension
        varDefEnsembleInfo(varID, isec1[54], isec1[53], isec1[52]);
      else if ( ISEC1_CenterID == 252 && isec1[36] ==   1 ) // MPIM local extension
        varDefEnsembleInfo(varID, isec1[38], isec1[39], isec1[37]);
    }

  if ( lmv ) varDefMissval(varID, FSEC3_MissVal);

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      int center    = ISEC1_CenterID;
      int subcenter = ISEC1_SubCenterID;
      int instID    = institutInq(center, subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef(center, subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      int modelID = modelInq(varInqInst(varID), ISEC1_ModelID, NULL);
      if ( modelID == CDI_UNDEFID )
	modelID = modelDef(varInqInst(varID), ISEC1_ModelID, NULL);
      varDefModel(varID, modelID);
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int tableID = tableInq(varInqModel(varID), ISEC1_CodeTable, NULL);
      if ( tableID == CDI_UNDEFID )
	tableID = tableDef(varInqModel(varID), ISEC1_CodeTable, NULL);
      varDefTable(varID, tableID);
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
}

static
void MCH_get_undef(int *isec1, double *undef_pds, double *undef_eps)
{
  /* 2010-01-13: Oliver Fuhrer */
  if ( ISEC1_CenterID == 215 ) {
    if (isec1[34] != 0 && isec1[34] != 255) {
      if (isec1[34] & 2) {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,-isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,-isec1[35]);
        }
        *undef_eps = pow(10.0,-isec1[35]-1);
      } else {
        if (isec1[34] & 1) {
          *undef_pds = -0.99*pow(10.0,+isec1[35]);
        } else {
          *undef_pds = +0.99*pow(10.0,+isec1[35]);
        }
        *undef_eps = pow(10.0,isec1[35]-1);
      }
    }
  }
}

static
void cgribexDecodeHeader(int *isec0, int *isec1, int *isec2, double *fsec2,
			 int *isec3, double *fsec3, int *isec4, double *fsec4,
			 int *gribbuffer, int recsize, int *lmv, int *iret)
{
  int ipunp = 0, iword = 0;

  memset(isec1, 0, 256*sizeof(int));
  memset(isec2, 0, 32*sizeof(int));

  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
	   ipunp, (int *) gribbuffer, recsize, &iword, "J", iret);

  if ( !(ISEC1_Sec2Or3Flag & 128) ) isec2[0] = -1; // default generic grid

  *lmv = 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;
      MCH_get_undef(isec1, &undef_pds, &undef_eps);
      FSEC3_MissVal = undef_pds;
      *lmv = 1;
    }
}

static
compvar_t cgribexVarSet(int param, int level1, int level2, int leveltype, int trange)
{
  int tsteptype = cgribexGetTsteptype(trange);

  compvar_t compVar;
  compVar.param     = param;
  compVar.level1    = level1;
  compVar.level2    = level2;
  compVar.ltype     = leveltype;
  compVar.tsteptype = tsteptype;

  return compVar;
}

static inline int
cgribexVarCompare(compvar_t compVar, record_t record, int flag)
{
  bool vinst = compVar.tsteptype == TSTEP_INSTANT || compVar.tsteptype == TSTEP_INSTANT2 || compVar.tsteptype == TSTEP_INSTANT3;
  bool rinst = record.tsteptype == TSTEP_INSTANT || record.tsteptype == TSTEP_INSTANT2 || record.tsteptype == TSTEP_INSTANT3;
  int tstepDiff = (!((flag == 0) & (vinst && rinst)))
                & (compVar.tsteptype != record.tsteptype);
  int rstatus = (compVar.param != record.param)
    |           (compVar.level1 != record.ilevel)
    |           (compVar.level2 != record.ilevel2)
    |           (compVar.ltype != record.ltype)
    |           tstepDiff;
  return rstatus;
}
#endif

#define gribWarning(text, nrecs, timestep, paramstr, level1, level2) \
            Warning("Record %2d (id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, paramstr, level1, level2, timestep, text)

#if  defined  (HAVE_LIBCGRIBEX)

static inline void
cgribexScanTsFixNtsteps(stream_t *streamptr, off_t recpos)
{
  if ( streamptr->ntsteps == -1 )
    {
      int tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = true;
      streamptr->tsteps[tsID].position = recpos;
    }
}

static inline void
cgribexScanTsConstAdjust(stream_t *streamptr, taxis_t *taxis)
{
  int vlistID = streamptr->vlistID;
  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
	{
	  streamptr->ntsteps = 0;
	  for ( int varID = 0; varID < streamptr->nvars; varID++ )
            vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
    }
}


int cgribexScanTimestep1(stream_t *streamptr)
{
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  void *gribbuffer = NULL;
  size_t buffersize = 0;
  int rstatus;
  int param = 0;
  int leveltype = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0 = { LONG_MIN, LONG_MIN };
  size_t readsize;
  unsigned nrecords, recID;
  int nrecs_scanned = 0;
  int datatype;
  size_t recsize = 0;
  bool warn_time = true;
  bool warn_numavg = true;
  int taxisID = -1;
  int rdate = 0, rtime = 0, tunit = 0;
  bool fcast = false;
  int vlistID;
  int comptype;
  size_t unzipsize;
  char paramstr[32];
  int nskip = cdiSkipRecords;

  streamptr->curTsID = 0;

  int *isec0 = streamptr->record->sec0;
  int *isec1 = streamptr->record->sec1;
  int *isec2 = streamptr->record->sec2;
  int *isec3 = streamptr->record->sec3;
  int *isec4 = streamptr->record->sec4;

  int tsID  = tstepsNewEntry(streamptr);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  int fileID = streamptr->fileID;

  while ( nskip-- > 0 )
    {
      recsize = gribGetSize(fileID);
      if ( recsize == 0 )
	Error("Skipping of %d records failed!", cdiSkipRecords);

      recpos  = fileGetPos(fileID);
      fileSetPos(fileID, (off_t)recsize, SEEK_CUR);
    }

  unsigned nrecs = 0;
  while ( true )
    {
      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);

      if ( recsize == 0 )
	{
	  if ( nrecs == 0 )
	    Error("No GRIB records found!");

	  streamptr->ntsteps = 1;
	  break;
	}
      if ( recsize > buffersize )
	{
	  buffersize = recsize;
	  gribbuffer = Realloc(gribbuffer, buffersize);
	}

      readsize = recsize;
      rstatus = gribRead(fileID, (unsigned char *)gribbuffer, &readsize);
      if ( rstatus ) break;

      comptype = CDI_COMPRESS_NONE;
      if ( gribGetZip(recsize, (unsigned char *)gribbuffer, &unzipsize) > 0 )
	{
	  comptype = CDI_COMPRESS_SZIP;
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( buffersize < unzipsize )
	    {
	      buffersize = unzipsize;
	      gribbuffer = Realloc(gribbuffer, buffersize);
	    }
	}

      nrecs_scanned++;
      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, (int)recsize, &lmv, &iret);

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      cgribexGetLevel(isec1, &leveltype, &level1, &level2);

      gribDateTime(isec1, &vdate, &vtime);

      datatype = (ISEC4_NumBits > 0 && ISEC4_NumBits <= 32) ? ISEC4_NumBits : CDI_DATATYPE_PACK;

      if ( nrecs == 0 )
	{
	  datetime0.date = vdate;
	  datetime0.time = vtime;
	  rdate = gribRefDate(isec1);
	  rtime = gribRefTime(isec1);
	  tunit = cgribexGetTimeUnit(isec1);
	  fcast = cgribexTimeIsFC(isec1);
	}
      else
	{
	  datetime.date = vdate;
	  datetime.time = vtime;
	  compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      if ( cgribexVarCompare(compVar, streamptr->tsteps[0].records[recID], 0) == 0 ) break;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      if ( recID < nrecs ) break;
	      if ( warn_time )
		if ( datetimeCmp(datetime, datetime0) != 0 )
		  {
                    gribWarning("Inconsistent verification time!", nrecs_scanned, tsID+1, paramstr, level1, level2);
		    warn_time = false;
		  }
	    }
	  else
	    {
	      if ( datetimeCmp(datetime, datetime0) != 0 ) break;

	      if ( recID < nrecs )
		{
		  gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr, level1, level2);
		  continue;
		}
	    }
	}

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
	    {
	      Warning("Changing numavg from %d to %d not supported!", taxis->numavg, ISEC1_AvgNum);
	      warn_numavg = false;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      nrecs++;

      if ( CDI_Debug )
	Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

      cgribexAddRecord(streamptr, param, isec1, isec2, fsec2, fsec3,
		       isec4, recsize, recpos, datatype, comptype, lmv, iret);
    }

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return CDI_EUFSTRUCT;

  cdi_generate_vars(streamptr);

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
      taxis->unit  = tunit;
    }

  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  nrecords = (unsigned)streamptr->tsteps[0].nallrecs;
  if ( nrecords < (unsigned)streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = (int)nrecords;
      streamptr->tsteps[0].records =
      (record_t *) Realloc(streamptr->tsteps[0].records, nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) Malloc(nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = (int)nrecords;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = (int)recID;

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  cgribexScanTsFixNtsteps(streamptr, recpos);
  cgribexScanTsConstAdjust(streamptr, taxis);

  return 0;
}


int cgribexScanTimestep2(stream_t * streamptr)
{
  int rstatus = 0;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  int param = 0;
  int leveltype = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0 = { LONG_MIN, LONG_MIN };
  int varID, gridID;
  size_t readsize;
  int nrecs, recID;
  size_t recsize = 0;
  bool warn_numavg = true;
  int tsteptype;
  size_t unzipsize;
  char paramstr[32];

  streamptr->curTsID = 1;

  int *isec0 = streamptr->record->sec0;
  int *isec1 = streamptr->record->sec1;
  int *isec2 = streamptr->record->sec2;
  int *isec3 = streamptr->record->sec3;
  int *isec4 = streamptr->record->sec4;

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
  if ( nrecords ) streamptr->tsteps[1].recIDs = (int *) Malloc((size_t)nrecords * sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      varID = streamptr->tsteps[0].records[recID].varID;
      streamptr->tsteps[tsID].records[recID].position =	streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     =	streamptr->tsteps[0].records[recID].size;
    }

  int nrecs_scanned = nrecords;
  int rindex = 0;
  while ( true )
    {
      if ( rindex > nrecords ) break;

      recsize = gribGetSize(fileID);
      recpos  = fileGetPos(fileID);
      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      if ( recsize > buffersize )
	{
	  buffersize = recsize;
	  gribbuffer = Realloc(gribbuffer, buffersize);
	}

      readsize = recsize;
      rstatus = gribRead(fileID, (unsigned char *)gribbuffer, &readsize);
      if ( rstatus ) break;

      if ( gribGetZip(recsize, (unsigned char *)gribbuffer, &unzipsize) > 0 )
	{
	  unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	  if ( buffersize < unzipsize )
	    {
	      buffersize = unzipsize;
	      gribbuffer = Realloc(gribbuffer, buffersize);
	    }
	}

      cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			  (int *) gribbuffer, (int)recsize, &lmv, &iret);

      nrecs_scanned++;

      param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      cgribexGetLevel(isec1, &leveltype, &level1, &level2);

      gribDateTime(isec1, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;
	      taxis->rdate = gribRefDate(isec1);
	      taxis->rtime = gribRefTime(isec1);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->unit  = cgribexGetTimeUnit(isec1);
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);

      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
        	(taxis->numavg != ISEC1_AvgNum) )
	    {
	  /*
	      Warning("Changing numavg from %d to %d not supported!", taxis->numavg, ISEC1_AvgNum);
	  */
	      warn_numavg = false;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}

      datetime.date  = vdate;
      datetime.time  = vtime;

      compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	}

      if ( recID == nrecords )
	{
	  gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, paramstr, level1, level2);
	  return CDI_EUFSTRUCT;
	}

      if ( cdiInventoryMode == 1 )
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      break;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = true;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	}
      else
	{
	  if ( streamptr->tsteps[tsID].records[recID].used )
	    {
	      if ( datetimeCmp(datetime, datetime0) != 0 ) break;

              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr, level1, level2);
	      continue;
	    }
	  else
	    {
	      streamptr->tsteps[tsID].records[recID].used = true;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	}

      if ( CDI_Debug )
	Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

      streamptr->tsteps[tsID].records[recID].size = recsize;

      if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return CDI_EUFSTRUCT;
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      varID = streamptr->tsteps[tsID].records[recID].varID;
      gridID = vlistInqVarGrid(vlistID, varID);
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}

      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      rindex++;
    }

  nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  varID = streamptr->tsteps[tsID].records[recID].varID;
          vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  cgribexScanTsFixNtsteps(streamptr, recpos);

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  return rstatus;
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
int cgribexScanTimestep(stream_t * streamptr)
{
  int rstatus = 0;
  double fsec2[512], fsec3[2], *fsec4 = NULL;
  int lmv = 0, iret = 0;
  size_t recsize = 0;
  off_t recpos = 0;
  void *gribbuffer;
  size_t buffersize = 0;
  int fileID;
  int param = 0;
  int leveltype = 0, level1 = 0, level2 = 0, vdate = 0, vtime = 0;
  DateTime datetime, datetime0 = { LONG_MIN, LONG_MIN };
  int vrecID, recID;
  bool warn_numavg = true;
  size_t readsize;
  int taxisID = -1;
  int rindex, nrecs = 0;
  int nrecs_scanned;
  size_t unzipsize;
  char paramstr[32];

  /*
  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }
  */
  int *isec0 = streamptr->record->sec0;
  int *isec1 = streamptr->record->sec1;
  int *isec2 = streamptr->record->sec2;
  int *isec3 = streamptr->record->sec3;
  int *isec4 = streamptr->record->sec4;

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      gribbuffer = streamptr->record->buffer;
      buffersize = streamptr->record->buffersize;

      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) Malloc((size_t)nrecs * sizeof (int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      nrecs_scanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs*(tsID-1);
      rindex = 0;
      while ( true )
	{
	  if ( rindex > nrecs ) break;

	  recsize = gribGetSize(fileID);
	  recpos  = fileGetPos(fileID);
	  if ( recsize == 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }
	  if ( recsize > 0 && recsize > buffersize )
	    {
	      buffersize = recsize;
	      gribbuffer = Realloc(gribbuffer, buffersize);
	    }

	  if ( rindex >= nrecs ) break;

	  readsize = recsize;
	  rstatus = gribRead(fileID, (unsigned char *)gribbuffer, &readsize);
	  if ( rstatus )
	    {
	      Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
                      streamptr->tsteps[tsID].recordSize);
	      break;
	    }

	  if ( gribGetZip(recsize, (unsigned char *)gribbuffer, &unzipsize) > 0 )
	    {
	      unzipsize += 100; /* need 0 to 1 bytes for rounding of bds */
	      if ( buffersize < unzipsize )
		{
		  buffersize = unzipsize;
		  gribbuffer = Realloc(gribbuffer, buffersize);
		}
	    }

	  cgribexDecodeHeader(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4,
			      (int *) gribbuffer, (int)recsize, &lmv, &iret);

          nrecs_scanned++;

	  param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
          cdiParamToString(param, paramstr, sizeof(paramstr));

          cgribexGetLevel(isec1, &leveltype, &level1, &level2);

	  gribDateTime(isec1, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
              int vlistID = streamptr->vlistID;
	      taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;
		  taxis->rdate = gribRefDate(isec1);
		  taxis->rtime = gribRefTime(isec1);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->unit  = cgribexGetTimeUnit(isec1);
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }

	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
	      /*
	          Warning("Changing numavg from %d to %d not supported!", streamptr->tsteps[tsID].taxis.numavg, ISEC1_AvgNum);
	      */
		  warn_numavg = false;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }

	  datetime.date  = vdate;
	  datetime.time  = vtime;

	  compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);

	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, paramstr, level1, level2);

	      if ( cdiInventoryMode == 1 )
		return CDI_EUFSTRUCT;
	      else
		continue;
	    }

	  if ( cdiInventoryMode == 1 )
	    {
	      streamptr->tsteps[tsID].records[recID].used = true;
	      streamptr->tsteps[tsID].recIDs[rindex] = recID;
	    }
	  else
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  char paramstr_[32];
		  cdiParamToString(param, paramstr_, sizeof(paramstr_));

		  if ( datetimeCmp(datetime, datetime0) != 0 ) break;

		  if ( CDI_Debug )
                    gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, paramstr_, level1, level2);

		  continue;
		}
	      else
		{
		  streamptr->tsteps[tsID].records[recID].used = true;
		  streamptr->tsteps[tsID].recIDs[rindex] = recID;
		}
	    }

	  if ( CDI_Debug )
            Message("Read record %2d (id=%s lev1=%d lev2=%d) %8d %6d", nrecs_scanned, paramstr, level1, level2, vdate, vtime);

	  if ( cgribexVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, level1);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  rindex++;
	}

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  cdiParamToString(streamptr->tsteps[tsID].records[recID].param, paramstr, sizeof(paramstr));
	  gribWarning("Paramameter not found!", nrecs_scanned, tsID+1, paramstr,
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

  rstatus = (int)streamptr->ntsteps;

  return rstatus;
}
#endif

#ifdef gribWarning
#undef gribWarning
#endif

#if  defined  (HAVE_LIBCGRIBEX)
int cgribexDecode(int memtype, void *gribbuffer, int gribsize, void *data, long datasize,
		  int unreduced, int *nmiss, double missval)
{
  int status = 0;
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  float fsec2f[512], fsec3f[2];
  double fsec2[512], fsec3[2];
  char hoper[2];

  if ( unreduced ) strcpy(hoper, "R");
  else             strcpy(hoper, "D");

  FSEC3_MissVal = missval;

  if ( memtype == MEMTYPE_FLOAT )
    gribExSP(isec0, isec1, isec2, fsec2f, isec3, fsec3f, isec4, (float*) data,
             (int) datasize, (int*) gribbuffer, gribsize, &iword, hoper, &iret);
  else
    gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double*) data,
             (int) datasize, (int*) gribbuffer, gribsize, &iword, hoper, &iret);

  *nmiss = (ISEC1_Sec2Or3Flag & 64) ? ISEC4_NumValues - ISEC4_NumNonMissValues : 0;

  if ( ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255) )
    {
      double undef_pds, undef_eps;
      MCH_get_undef(isec1, &undef_pds, &undef_eps);

      *nmiss = 0;
      if ( memtype == MEMTYPE_FLOAT )
        {
          float *restrict dataf = (float*) data;
          for ( long i = 0; i < datasize; i++ )
            if ( (fabs(dataf[i]-undef_pds) < undef_eps) || IS_EQUAL(dataf[i],FSEC3_MissVal) ) {
              dataf[i] = (float)missval;
              (*nmiss)++;
            }
        }
      else
        {
          double *restrict datad = (double*) data;
          for ( long i = 0; i < datasize; i++ )
            if ( (fabs(datad[i]-undef_pds) < undef_eps) || IS_EQUAL(datad[i],FSEC3_MissVal) ) {
              datad[i] = missval;
              (*nmiss)++;
            }
        }
    }

  return status;
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
static
void cgribexDefInstitut(int *isec1, int vlistID, int varID)
{
  int instID = (vlistInqInstitut(vlistID) != CDI_UNDEFID) ? vlistInqInstitut(vlistID) : vlistInqVarInstitut(vlistID, varID);
  if ( instID != CDI_UNDEFID )
    {
      ISEC1_CenterID    = institutInqCenter(instID);
      ISEC1_SubCenterID = institutInqSubcenter(instID);
    }
}

static
void cgribexDefModel(int *isec1, int vlistID, int varID)
{
  int modelID = (vlistInqModel(vlistID) != CDI_UNDEFID) ? vlistInqModel(vlistID) : vlistInqVarModel(vlistID, varID);
  if ( modelID != CDI_UNDEFID )
    ISEC1_ModelID = modelInqGribID(modelID);
}

static
void cgribexDefParam(int *isec1, int param)
{
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  if ( pnum < 0 ) pnum = -pnum;

  static bool lwarn_pdis = true;
  if ( pdis != 255 && lwarn_pdis )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
      lwarn_pdis = false;
    }

  static bool lwarn_pnum = true;
  if ( pnum > 255 && lwarn_pnum )
    {
      Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum%256);
      lwarn_pnum = false;
      pnum = pnum%256;
    }

  ISEC1_CodeTable = pcat;
  ISEC1_Parameter = pnum;
}

static
int cgribexDefTimerange(int tsteptype, int factor, int calendar,
			int rdate, int rtime, int vdate, int vtime, int *pip1, int *pip2)
{
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday1, &secofday1);

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  int timerange = -1;
  int ip1 = 0, ip2 = 0;
  if ( !(int)(fmod(days*86400.0 + secs, factor)) )
    {
      int ip = (int) ((days*86400.0 + secs)/factor);

      if ( (ip > 255) && (tsteptype == TSTEP_INSTANT) ) tsteptype = TSTEP_INSTANT3;

      switch ( tsteptype )
	{
	case TSTEP_INSTANT:  timerange =  0; ip1 = ip; ip2 = 0;  break;
	case TSTEP_INSTANT2: timerange =  1; ip1 = 0;  ip2 = 0;  break;
	case TSTEP_RANGE:    timerange =  2; ip1 = 0;  ip2 = ip; break;
	case TSTEP_AVG:      timerange =  3; ip1 = 0;  ip2 = ip; break;
	case TSTEP_ACCUM:    timerange =  4; ip1 = 0;  ip2 = ip; break;
	case TSTEP_DIFF:     timerange =  5; ip1 = 0;  ip2 = ip; break;
	case TSTEP_INSTANT3:
	default:             timerange = 10; ip1 = ip/256; ip2 = ip%256; break;
	}
    }

  *pip1 = ip1;
  *pip2 = ip2;

  return timerange;
}

static
int cgribexDefDateTime(int *isec1, int timeunit, int date, int time)
{
  int year, month, day, hour, minute, second;
  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  int century =  year / 100;

  ISEC1_Year = year - century*100;

  if ( year < 0 )
    {
      century = -century;
      ISEC1_Year = -ISEC1_Year;
    }

  if ( ISEC1_Year == 0 )
    {
      century -= 1;
      ISEC1_Year = 100;
    }

  century += 1;
  if ( year < 0 ) century = -century;

  ISEC1_Month  = month;
  ISEC1_Day    = day;
  ISEC1_Hour   = hour;
  ISEC1_Minute = minute;

  ISEC1_Century = century;

  int factor = 1;
  switch (timeunit)
    {
    case TUNIT_MINUTE:    factor =    60; ISEC1_TimeUnit = ISEC1_TABLE4_MINUTE;    break;
    case TUNIT_QUARTER:   factor =   900; ISEC1_TimeUnit = ISEC1_TABLE4_QUARTER;   break;
    case TUNIT_30MINUTES: factor =  1800; ISEC1_TimeUnit = ISEC1_TABLE4_30MINUTES; break;
    case TUNIT_HOUR:      factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    case TUNIT_3HOURS:    factor = 10800; ISEC1_TimeUnit = ISEC1_TABLE4_3HOURS;    break;
    case TUNIT_6HOURS:    factor = 21600; ISEC1_TimeUnit = ISEC1_TABLE4_6HOURS;    break;
    case TUNIT_12HOURS:   factor = 43200; ISEC1_TimeUnit = ISEC1_TABLE4_12HOURS;   break;
    case TUNIT_DAY:       factor = 86400; ISEC1_TimeUnit = ISEC1_TABLE4_DAY;       break;
    default:              factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    }

  return factor;
}

static
void cgribexDefTime(int *isec1, int vdate, int vtime, int tsteptype, int numavg, int taxisID)
{
  int timetype = TAXIS_ABSOLUTE;
  int timerange = 0;
  int timeunit = TUNIT_HOUR;

  if ( taxisID != -1 )
    {
      timetype = taxisInqType(taxisID);
      timeunit = taxisInqTunit(taxisID);
    }

  if ( timetype == TAXIS_RELATIVE )
    {
      int calendar = taxisInqCalendar(taxisID);
      int rdate    = taxisInqRdate(taxisID);
      int rtime    = taxisInqRtime(taxisID);

      int factor = cgribexDefDateTime(isec1, timeunit, rdate, rtime);
      int ip1 = 0, ip2 = 0;
      timerange = cgribexDefTimerange(tsteptype, factor, calendar,
				      rdate, rtime, vdate, vtime, &ip1, &ip2);

      if ( ip2 > 0xFF && timeunit < TUNIT_YEAR )
        {
          timeunit++;
          factor = cgribexDefDateTime(isec1, timeunit, rdate, rtime);
          timerange = cgribexDefTimerange(tsteptype, factor, calendar,
                                          rdate, rtime, vdate, vtime, &ip1, &ip2);
        }

      if ( timerange == -1 || timerange == 3 )
	{
	  timetype = TAXIS_ABSOLUTE;
	}
      /*
      else if ( timerange == 10 )
	{
	  if ( ip1 < 0 || ip1 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFFFF ) timetype = TAXIS_ABSOLUTE;
	}
      */
      else
	{
	  if ( ip1 < 0 || ip1 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	  if ( ip2 < 0 || ip2 > 0xFF   ) timetype = TAXIS_ABSOLUTE;
	}

      ISEC1_TimeRange   = timerange;
      ISEC1_TimePeriod1 = ip1;
      ISEC1_TimePeriod2 = ip2;
    }

  if ( timetype == TAXIS_ABSOLUTE )
    {
      (void) cgribexDefDateTime(isec1, timeunit, vdate, vtime);

      /*
      if ( numavg > 0 )
	ISEC1_TimeRange = 0;
      else
      */
      if ( ISEC1_TimeRange != 3 )
	ISEC1_TimeRange   = 10;

      ISEC1_TimePeriod1 = 0;
      ISEC1_TimePeriod2 = 0;
    }

  ISEC1_AvgNum         = numavg;
  ISEC1_AvgMiss        = 0;
  ISEC1_DecScaleFactor = 0;
}

static
void cgribexDefGrid(int *isec1, int *isec2, double *fsec2, int *isec4, int gridID)
{
  bool lrotated = false;
  bool lcurvi = false;

  memset(isec2, 0, 16*sizeof(int));

  ISEC1_Sec2Or3Flag = 128;

  int gridtype = gridInqType(gridID);

  ISEC1_GridDefinition = 255;

  if ( gridtype == GRID_GENERIC )
    {
      int gridsize = gridInqSize(gridID);
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
          if ( lwarning && gridInqSize(gridID) > 1 )
            {
              lwarning = false;
              Warning("Curvilinear grid is unsupported in GRIB1! Created wrong Grid Description Section!");
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

  ISEC2_Reduced  = false;
  ISEC2_ScanFlag = 0;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
	double xfirst = 0, xlast = 0, xinc = 0;
	double yfirst = 0, ylast = 0, yinc = 0;

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          ISEC2_GridType = GRIB1_GTYPE_GAUSSIAN;
        else if ( gridtype == GRID_LONLAT && lrotated )
	  ISEC2_GridType = GRIB1_GTYPE_LATLON_ROT;
	else
	  ISEC2_GridType = GRIB1_GTYPE_LATLON;

	int nlon = gridInqXsize(gridID);
	int nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_Reduced = true;
	    nlon = 0;
	    gridInqRowlon(gridID, ISEC2_RowLonPtr);
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

	ISEC2_NumLon   = nlon;
	ISEC2_NumLat   = nlat;
	ISEC2_FirstLat = (int)lround(yfirst*1000);
	ISEC2_LastLat  = (int)lround(ylast*1000);
	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    ISEC2_FirstLon = 0;
	    ISEC2_LastLon  = (int)lround(1000*(360.-360./(nlat*2)));
	    ISEC2_LonIncr  = (int)lround(1000*360./(nlat*2));
	  }
	else
	  {
	    ISEC2_FirstLon = (int)lround(xfirst*1000);
	    ISEC2_LastLon  = (int)lround(xlast*1000);
	    ISEC2_LonIncr  = (int)lround(xinc*1000);
	  }

	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            ISEC2_NumPar = np;
          }
	else
	  {
	    ISEC2_LatIncr = (int)lround(yinc*1000);
	  }

	if ( ISEC2_NumLon > 1 && ISEC2_NumLat == 1 )
	  if ( ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0 ) ISEC2_LatIncr = ISEC2_LonIncr;

	if ( ISEC2_NumLon == 1 && ISEC2_NumLat > 1 )
	  if ( ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0 ) ISEC2_LonIncr = ISEC2_LatIncr;

        ISEC2_ResFlag = 0;
        if ( ISEC2_LatIncr && ISEC2_LonIncr )   gribbyte_set_bit(&ISEC2_ResFlag, 1);
        if ( gridInqUvRelativeToGrid(gridID) )  gribbyte_set_bit(&ISEC2_ResFlag, 5);

	if ( lrotated )
          {
            double xpole = 0, ypole = 0, angle = 0;
            gridInqParamRLL(gridID, &xpole, &ypole, &angle);

	    ISEC2_LatSP = - (int)lround(ypole * 1000);
	    ISEC2_LonSP =   (int)lround((xpole + 180) * 1000);
            if ( fabs(angle) > 0 ) angle = -angle;
            FSEC2_RotAngle = angle;
          }

        ISEC2_ScanFlag = 0;
	if ( ISEC2_LastLon < ISEC2_FirstLon ) gribbyte_set_bit(&ISEC2_ScanFlag, 1); // East -> West
	if ( ISEC2_LastLat > ISEC2_FirstLat ) gribbyte_set_bit(&ISEC2_ScanFlag, 2); // South -> North

	break;
      }
    case GRID_LCC:
      {
	int xsize = gridInqXsize(gridID);
        int ysize = gridInqYsize(gridID);

        double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;
	gridInqParamLCC(gridID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);
	gridVerifyGribParamLCC(grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);
        bool lsouth = (lat_1 < 0);
        if ( lsouth ) { lat_1 = -lat_2; lat_2 = -lat_2; }

        double xinc = gridInqXinc(gridID);
        double yinc = gridInqYinc(gridID);

	ISEC2_GridType = GRIB1_GTYPE_LCC;
	ISEC2_NumLon   = xsize;
	ISEC2_NumLat   = ysize;
	ISEC2_FirstLon       = (int)lround(xval_0 * 1000);
	ISEC2_FirstLat       = (int)lround(yval_0 * 1000);
	ISEC2_Lambert_Lov    = (int)lround(lon_0 * 1000);
	ISEC2_Lambert_LatS1  = (int)lround(lat_1 * 1000);
	ISEC2_Lambert_LatS2  = (int)lround(lat_2 * 1000);
	ISEC2_Lambert_dx     = (int)lround(xinc);
	ISEC2_Lambert_dy     = (int)lround(yinc);
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_LonSP  = 0;
	ISEC2_Lambert_ProjFlag = 0;
        if ( lsouth ) gribbyte_set_bit(&ISEC2_Lambert_ProjFlag, 1);

        bool earthIsOblate = (IS_EQUAL(a, 6378160.) && IS_EQUAL(rf, 297.));
        ISEC2_ResFlag = 0;
        if ( ISEC2_Lambert_dx && ISEC2_Lambert_dy ) gribbyte_set_bit(&ISEC2_ResFlag, 1);
        if ( earthIsOblate )                        gribbyte_set_bit(&ISEC2_ResFlag, 2);
        if ( gridInqUvRelativeToGrid(gridID) )      gribbyte_set_bit(&ISEC2_ResFlag, 5);

	ISEC2_ScanFlag = 0;
        gribbyte_set_bit(&ISEC2_ScanFlag, 2); // South -> North

	break;
      }
    case GRID_SPECTRAL:
      {
	ISEC2_GridType = GRIB1_GTYPE_SPECTRAL;
	ISEC2_PentaJ   = gridInqTrunc(gridID);
	ISEC2_PentaK   = ISEC2_PentaJ;
	ISEC2_PentaM   = ISEC2_PentaJ;
	ISEC2_RepType  = 1;
	isec4[2]       = 128;
	if ( gridInqComplexPacking(gridID) && ISEC2_PentaJ >= 21 )
	  {
	    ISEC2_RepMode  = 2;
	    isec4[3]       = 64;
	    isec4[16]      = 0;
	    isec4[17]      = 20;
	    isec4[18]      = 20;
	    isec4[19]      = 20;
	  }
	else
	  {
	    ISEC2_RepMode  = 1;
	    isec4[3]       = 0;
	  }
	break;
      }
    case GRID_GME:
      {
	ISEC2_GridType   = GRIB1_GTYPE_GME;
        int nd = 0, ni = 0, ni2 = 0, ni3 = 0;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
	ISEC2_GME_ND     = nd;
	ISEC2_GME_NI     = ni;
	ISEC2_GME_NI2    = ni2;
	ISEC2_GME_NI3    = ni3;
	ISEC2_GME_AFlag  = 0;
	ISEC2_GME_LatPP  = 90000;
	ISEC2_GME_LonPP  = 0;
	ISEC2_GME_LonMPL = 0;
	ISEC2_GME_BFlag  = 0;
	break;
      }
    case GRID_GENERIC:
      {
        ISEC1_Sec2Or3Flag = 0;
	break;
      }
    default:
      {
        ISEC1_Sec2Or3Flag = 0;
	Warning("CGRIBEX library doesn't support %s grids, grid information will be lost!", gridNamePtr(gridtype));
	break;
      }
    }


  if ( cdiGribChangeModeUvRelativeToGrid.active )
    {
      // this will overrule/change the UvRelativeToGrid flag;
      // typically when the wind is rotated with respect to north pole
      bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);
      if      ( uvRelativeToGrid && !cdiGribChangeModeUvRelativeToGrid.mode )
        gribbyte_clear_bit(&ISEC2_ResFlag, 5);
      else if ( !uvRelativeToGrid && cdiGribChangeModeUvRelativeToGrid.mode )
        gribbyte_set_bit(&ISEC2_ResFlag, 5);
    }
}

static
void isec1DefLevel(int *isec1, int leveltype, int level1, int level2)
{
  ISEC1_LevelType = leveltype;
  ISEC1_Level1    = level1;
  ISEC1_Level2    = level2;
}

static
void cgribexDefLevel(int *isec1, int *isec2, double *fsec2, int zaxisID, int levelID)
{
  char units[CDI_MAX_NAME];
  static bool lwarning_vct = true;

  int zaxistype = zaxisInqType(zaxisID);
  int ltype = zaxisInqLtype(zaxisID);

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s",
	      zaxisNamePtr(zaxistype), zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  ISEC2_NumVCP = 0;

  int grib_ltype = zaxisTypeToGrib1ltype(zaxistype);

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_ALTITUDE:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
        isec1DefLevel(isec1, grib_ltype, (int)(zaxisInqLevel(zaxisID, levelID)), 0);
	break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_ATMOSPHERE:
      {
        isec1DefLevel(isec1, grib_ltype, 0, 0);
	break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
          isec1DefLevel(isec1, GRIB1_LTYPE_HYBRID_LAYER, (int)(zaxisInqLbound(zaxisID, levelID)),
                        (int)(zaxisInqUbound(zaxisID, levelID)));
	else
          isec1DefLevel(isec1, GRIB1_LTYPE_HYBRID, (int)(zaxisInqLevel(zaxisID, levelID)), 0);

	int vctsize = zaxisInqVctSize(zaxisID);
	if ( vctsize > 255 )
	  {
	    ISEC2_NumVCP = 0;
	    if ( lwarning_vct )
	      {
		Warning("VCT size of %d is too large (maximum is 255). Set to 0!", vctsize);
		lwarning_vct = false;
	      }
	  }
	else
	  {
	    ISEC2_NumVCP = vctsize;
	    zaxisInqVct(zaxisID, &fsec2[10]);
	  }
	break;
      }
    case ZAXIS_PRESSURE:
      {
	double level = zaxisInqLevel(zaxisID, levelID);
	if ( level < 0 ) Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( (units[0] != 'P') | (units[1] != 'a') ) level *= 100;

	double dum;
	if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
          grib_ltype = GRIB1_LTYPE_ISOBARIC_PA;
	else
          level = level/100;

        isec1DefLevel(isec1, grib_ltype, (int) level, 0);
	break;
      }
    case ZAXIS_HEIGHT:
      {
	double level = zaxisInqLevel(zaxisID, levelID);

	zaxisInqUnits(zaxisID, units);
        if ( units[1] == 'm' && !units[2] )
          {
            if      ( units[0] == 'c' ) level *= 0.01;
            else if ( units[0] == 'd' ) level *= 0.1;
            else if ( units[0] == 'k' ) level *= 1000;
          }

        isec1DefLevel(isec1, grib_ltype, (int) level, 0);
	break;
      }
    case ZAXIS_SIGMA:
      {
	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
          isec1DefLevel(isec1, GRIB1_LTYPE_SIGMA_LAYER, (int)(zaxisInqLbound(zaxisID, levelID)),
                        (int)(zaxisInqUbound(zaxisID, levelID)));
	else
          isec1DefLevel(isec1, GRIB1_LTYPE_SIGMA, (int)(zaxisInqLevel(zaxisID, levelID)), 0);

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	zaxisInqUnits(zaxisID, units);

	double factor = 100; // default: meter
        if      ( units[0] == 'm' && units[1] == 'm' ) factor =   0.1;
        else if ( units[0] == 'c' && units[1] == 'm' ) factor =   1;
        else if ( units[0] == 'd' && units[1] == 'm' ) factor =  10;

	if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
          isec1DefLevel(isec1, GRIB1_LTYPE_LANDDEPTH_LAYER, (int) (factor*zaxisInqLbound(zaxisID, levelID)),
                        (int) (factor*zaxisInqUbound(zaxisID, levelID)));
	else
          isec1DefLevel(isec1, GRIB1_LTYPE_LANDDEPTH, (int) (factor*zaxisInqLevel(zaxisID, levelID)), 0);

	break;
      }
    case ZAXIS_GENERIC:
      {
        isec1DefLevel(isec1, ltype, (int)(zaxisInqLevel(zaxisID, levelID)), 0);
	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}

static
void cgribexDefaultSec0(int *isec0)
{
  ISEC0_GRIB_Len     = 0;
  ISEC0_GRIB_Version = 0;
}

static
void cgribexDefaultSec1(int *isec1)
{
  ISEC1_CenterID    = 0;
  ISEC1_SubCenterID = 0;
  ISEC1_LocalFLag   = 0;
}

static
void cgribexDefaultSec4(int *isec4)
{
  for ( int i = 2; i <= 10; ++i ) isec4[i] = 0;
}

static
void cgribexDefEnsembleVar(int *isec1, int vlistID, int varID)
{
  int ensID, ensCount, forecast_type;

  /* For Ensemble info  */

  //Put1Byte(isec1[36]);        /* MPIM local GRIB use definition identifier  */
                                /*    (extension identifier)                  */
  //Put1Byte(isec1[37]);        /* type of ensemble forecast                  */
  //Put2Byte(isec1[38]);        /* individual ensemble member                 */
  //Put2Byte(isec1[39]);        /* number of forecasts in ensemble            */

  if ( vlistInqVarEnsemble(vlistID, varID, &ensID, &ensCount, &forecast_type) )
    {
      if ( ISEC1_CenterID == 252 )
        {
          ISEC1_LocalFLag = 1;
          isec1[36] = 1;

          isec1[37] =  forecast_type;
          isec1[38] =  ensID;
          isec1[39] =  ensCount;
        }
    }
}
#endif


#if  defined  (HAVE_LIBCGRIBEX)
size_t cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg,
		     long datasize, const void *data, int nmiss, void *gribbuffer, size_t gribbuffersize)
{
  int iret = 0, iword = 0;
  int isec0[2], isec1[4096], isec2[4096], isec3[2], isec4[512];
  float fsec2f[512], fsec3f[2];
  double fsec2[512], fsec3[2];

  memset(isec1, 0, 256*sizeof(int));
  fsec2[0] = 0; fsec2[1] = 0;
  fsec2f[0] = 0; fsec2f[1] = 0;

  int gribsize = (int)(gribbuffersize / sizeof(int));
  int param    = vlistInqVarParam(vlistID, varID);

  cgribexDefaultSec0(isec0);
  cgribexDefaultSec1(isec1);
  cgribexDefaultSec4(isec4);

  cgribexDefInstitut(isec1, vlistID, varID);
  cgribexDefModel(isec1, vlistID, varID);

  int datatype = vlistInqVarDatatype(vlistID, varID);

  cgribexDefParam(isec1, param);
  cgribexDefTime(isec1, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID));
  cgribexDefGrid(isec1, isec2, fsec2, isec4, gridID);
  cgribexDefLevel(isec1, isec2, fsec2, zaxisID, levelID);

  cgribexDefEnsembleVar(isec1, vlistID, varID);

  ISEC4_NumValues = gridInqSize(gridID);
  ISEC4_NumBits   = grbBitsPerValue(datatype);

  if ( nmiss > 0 )
    {
      FSEC3_MissVal = vlistInqVarMissval(vlistID, varID);
      ISEC1_Sec2Or3Flag |= 64;
    }

  if ( isec4[2] == 128 && isec4[3] == 64 )
    {
      if ( memtype == MEMTYPE_FLOAT )
        isec4[16] = (int) (1000*calculate_pfactor_float((const float*) data, ISEC2_PentaJ, isec4[17]));
      else
        isec4[16] = (int) (1000*calculate_pfactor_double((const double*) data, ISEC2_PentaJ, isec4[17]));
      if ( isec4[16] < -10000 ) isec4[16] = -10000;
      if ( isec4[16] >  10000 ) isec4[16] =  10000;
    }
  //printf("isec4[16] %d\n", isec4[16]);

  if ( memtype == MEMTYPE_FLOAT )
    {
      size_t numVCP = ISEC2_NumVCP > 0 ? (size_t)ISEC2_NumVCP : (size_t)0;
      for ( size_t i = 0; i < numVCP; ++i ) fsec2f[10+i] = (float)fsec2[10+i];
      fsec3f[ 1] = (float)fsec3[ 1];
    }

  if ( memtype == MEMTYPE_FLOAT )
    gribExSP(isec0, isec1, isec2, fsec2f, isec3, fsec3f, isec4, (float*) data,
             (int) datasize, (int*) gribbuffer, gribsize, &iword, "C", &iret);
  else
    gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double*) data,
             (int) datasize, (int*) gribbuffer, gribsize, &iword, "C", &iret);

  if ( iret ) Error("Problem during GRIB encode (errno = %d)!", iret);

  size_t nbytes = (size_t)iword * sizeof(int);
  return nbytes;
}


typedef struct
{
  void *gribbuffer;
  size_t gribbuffersize;
  unsigned char *pds;
  unsigned char *gds;
  unsigned char *bms;
  unsigned char *bds;
} cgribex_handle;


int grib1Sections(unsigned char *gribbuffer, long gribbufsize, unsigned char **pdsp,
		  unsigned char **gdsp, unsigned char **bmsp, unsigned char **bdsp, long *gribrecsize);

void *cgribex_handle_new_from_meassage(void *gribbuffer, size_t gribbuffersize)
{
  cgribex_handle *gh = (cgribex_handle*) Malloc(sizeof(cgribex_handle));
  gh->gribbuffer = NULL;
  gh->gribbuffersize = 0;
  gh->pds = NULL;

  if ( gribbuffersize && gribbuffer )
    {
      unsigned char *pds = NULL, *gds = NULL, *bms = NULL, *bds = NULL;
      long gribrecsize;
      int status = grib1Sections((unsigned char *)gribbuffer, (long)gribbuffersize, &pds, &gds, &bms, &bds, &gribrecsize);
      if ( status >= 0 )
        {
          gh->gribbuffer = gribbuffer;
          gh->gribbuffersize = gribbuffersize;
          gh->pds = pds;
          gh->gds = gds;
          gh->bms = bms;
          gh->bds = bds;
        }
    }

  return (void*)gh;
}


void cgribex_handle_delete(void *gh)
{
  if ( gh ) Free(gh);
}


void cgribexChangeParameterIdentification(void *gh, int code, int ltype, int lev)
{
  if ( !gh ) return;

  unsigned char *pds = ((cgribex_handle*)gh)->pds;
  if ( !pds ) return;

  pds[8]  = (unsigned char) code;
  pds[9]  = (unsigned char) ltype;
  pds[10] = (unsigned char) lev;
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
