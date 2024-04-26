#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "stream_cgribex.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "gribapi.h"
#include "file.h"
#include "cgribex.h"  /* gribZip gribGetZip gribGinfo */

int cdiDebugExt                        =  0;      //  Debug level for the KNMI extensions
#ifdef HIRLAM_EXTENSIONS
// *** RELATED to GRIB only ***
int cdiGribUseTimeRangeIndicator        = 0;       // normaly cdo looks in grib for attribute called "stepType"
                                                   // but NWP models such as Harmonie 37h1.2, use "timeRangeIndicator"
                                                   // where:  0: for instanteneous fields; 4: for accumulated fields
#endif // HIRLAM_EXTENSIONS


// Regarding operation to change parameter identification:
// change if cdiGribChangeParameterID.active
struct cdiGribParamChange cdiGribChangeParameterID;

// Used only for CDO module Selmulti
void streamGrbChangeParameterIdentification(int code, int ltype, int lev)
{
  // NOTE this is a "PROXY" function for gribapiChangeParameterIdentification();
  // This just sets the globals. There are probably better solutions to this.
  // The parameter change is done by function  gribapiChangeParameterIdentification() in stream_gribapi.c
  // Setting this control variable to true will cause calling fnc. gribapiChangeParameterIdentification later.
  // After grib attributes have been changed this variable goes to false.
  cdiGribChangeParameterID.active = true;
  cdiGribChangeParameterID.code = code;
  cdiGribChangeParameterID.ltype = ltype;
  cdiGribChangeParameterID.lev = lev;
}

struct cdiGribModeChange cdiGribChangeModeUvRelativeToGrid;

// Used only for CDO module WindTrans
void streamGrbChangeModeUvRelativeToGrid(int mode)
{
  cdiGribChangeModeUvRelativeToGrid.active = true;
  cdiGribChangeModeUvRelativeToGrid.mode = (mode > 0);
}

struct cdiGribScanModeChange cdiGribDataScanningMode;

void streamGrbDefDataScanningMode(int scanmode)
{
  cdiGribDataScanningMode.active = true;
  cdiGribDataScanningMode.value = scanmode;
}

int  streamGrbInqDataScanningMode(void)
{
  return cdiGribDataScanningMode.value;
}


int grib1ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB1_LTYPE_SURFACE:            zaxistype = ZAXIS_SURFACE;                break;
    case GRIB1_LTYPE_CLOUD_BASE:         zaxistype = ZAXIS_CLOUD_BASE;             break;
    case GRIB1_LTYPE_CLOUD_TOP:          zaxistype = ZAXIS_CLOUD_TOP;              break;
    case GRIB1_LTYPE_ISOTHERM0:          zaxistype = ZAXIS_ISOTHERM_ZERO;          break;
    case GRIB1_LTYPE_TOA:                zaxistype = ZAXIS_TOA;                    break;
    case GRIB1_LTYPE_SEA_BOTTOM:         zaxistype = ZAXIS_SEA_BOTTOM;             break;
    case GRIB1_LTYPE_ATMOSPHERE:         zaxistype = ZAXIS_ATMOSPHERE;             break;
    case GRIB1_LTYPE_MEANSEA:            zaxistype = ZAXIS_MEANSEA;                break;
    case GRIB1_LTYPE_99:
    case GRIB1_LTYPE_ISOBARIC_PA:
    case GRIB1_LTYPE_ISOBARIC:           zaxistype = ZAXIS_PRESSURE;               break;
    case GRIB1_LTYPE_HEIGHT:             zaxistype = ZAXIS_HEIGHT;                 break;
    case GRIB1_LTYPE_ALTITUDE:           zaxistype = ZAXIS_ALTITUDE;	           break;
    case GRIB1_LTYPE_SIGMA:
    case GRIB1_LTYPE_SIGMA_LAYER:        zaxistype = ZAXIS_SIGMA;	           break;
    case GRIB1_LTYPE_HYBRID:
    case GRIB1_LTYPE_HYBRID_LAYER:       zaxistype = ZAXIS_HYBRID;	           break;
    case GRIB1_LTYPE_LANDDEPTH:
    case GRIB1_LTYPE_LANDDEPTH_LAYER:    zaxistype = ZAXIS_DEPTH_BELOW_LAND;       break;
    case GRIB1_LTYPE_ISENTROPIC:         zaxistype = ZAXIS_ISENTROPIC;             break;
    case GRIB1_LTYPE_SEADEPTH:           zaxistype = ZAXIS_DEPTH_BELOW_SEA;        break;
    case GRIB1_LTYPE_LAKE_BOTTOM:        zaxistype = ZAXIS_LAKE_BOTTOM;            break;
    case GRIB1_LTYPE_SEDIMENT_BOTTOM:    zaxistype = ZAXIS_SEDIMENT_BOTTOM;        break;
    case GRIB1_LTYPE_SEDIMENT_BOTTOM_TA: zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;     break;
    case GRIB1_LTYPE_SEDIMENT_BOTTOM_TW: zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;     break;
    case GRIB1_LTYPE_MIX_LAYER:          zaxistype = ZAXIS_MIX_LAYER;              break;
    }

  return zaxistype;
}


int grib2ltypeToZaxisType(int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  switch ( grib_ltype )
    {
    case GRIB2_LTYPE_SURFACE:            zaxistype = ZAXIS_SURFACE;                break;
    case GRIB2_LTYPE_CLOUD_BASE:         zaxistype = ZAXIS_CLOUD_BASE;             break;
    case GRIB2_LTYPE_CLOUD_TOP:          zaxistype = ZAXIS_CLOUD_TOP;              break;
    case GRIB2_LTYPE_ISOTHERM0:          zaxistype = ZAXIS_ISOTHERM_ZERO;          break;
    case GRIB2_LTYPE_TOA:                zaxistype = ZAXIS_TOA;                    break;
    case GRIB2_LTYPE_SEA_BOTTOM:         zaxistype = ZAXIS_SEA_BOTTOM;             break;
    case GRIB2_LTYPE_ATMOSPHERE:         zaxistype = ZAXIS_ATMOSPHERE;             break;
    case GRIB2_LTYPE_MEANSEA:            zaxistype = ZAXIS_MEANSEA;                break;
    case GRIB2_LTYPE_ISOBARIC:           zaxistype = ZAXIS_PRESSURE;               break;
    case GRIB2_LTYPE_HEIGHT:             zaxistype = ZAXIS_HEIGHT;                 break;
    case GRIB2_LTYPE_ALTITUDE:           zaxistype = ZAXIS_ALTITUDE;               break;
    case GRIB2_LTYPE_SIGMA:              zaxistype = ZAXIS_SIGMA;                  break;
    case GRIB2_LTYPE_HYBRID:
 /* case GRIB2_LTYPE_HYBRID_LAYER: */    zaxistype = ZAXIS_HYBRID;                 break;
    case GRIB2_LTYPE_LANDDEPTH:
 /* case GRIB2_LTYPE_LANDDEPTH_LAYER: */ zaxistype = ZAXIS_DEPTH_BELOW_LAND;       break;
    case GRIB2_LTYPE_ISENTROPIC:         zaxistype = ZAXIS_ISENTROPIC;             break;
    case GRIB2_LTYPE_SNOW:               zaxistype = ZAXIS_SNOW;                   break;
    case GRIB2_LTYPE_SEADEPTH:           zaxistype = ZAXIS_DEPTH_BELOW_SEA;        break;
    case GRIB2_LTYPE_LAKE_BOTTOM:        zaxistype = ZAXIS_LAKE_BOTTOM;            break;
    case GRIB2_LTYPE_SEDIMENT_BOTTOM:    zaxistype = ZAXIS_SEDIMENT_BOTTOM;        break;
    case GRIB2_LTYPE_SEDIMENT_BOTTOM_TA: zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;     break;
    case GRIB2_LTYPE_SEDIMENT_BOTTOM_TW: zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;     break;
    case GRIB2_LTYPE_MIX_LAYER:          zaxistype = ZAXIS_MIX_LAYER;              break;
    case GRIB2_LTYPE_REFERENCE:          zaxistype = ZAXIS_REFERENCE;              break;
    }

  return zaxistype;
}


int zaxisTypeToGrib1ltype(int zaxistype)
{
  int grib_ltype = -1;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:               grib_ltype = GRIB1_LTYPE_SURFACE;            break;
    case ZAXIS_GENERIC:               grib_ltype = -1;                             break;
    case ZAXIS_HYBRID:                grib_ltype = -1;                             break;
    case ZAXIS_HYBRID_HALF:           grib_ltype = -1;                             break;
    case ZAXIS_PRESSURE:              grib_ltype = GRIB1_LTYPE_ISOBARIC;           break;
    case ZAXIS_HEIGHT:                grib_ltype = GRIB1_LTYPE_HEIGHT;             break;
    case ZAXIS_DEPTH_BELOW_SEA:       grib_ltype = GRIB1_LTYPE_SEADEPTH;           break;
    case ZAXIS_DEPTH_BELOW_LAND:      grib_ltype = GRIB1_LTYPE_LANDDEPTH;          break;
    case ZAXIS_ISENTROPIC:            grib_ltype = GRIB1_LTYPE_ISENTROPIC;         break;
    case ZAXIS_TRAJECTORY:            grib_ltype = -1;                             break;
    case ZAXIS_ALTITUDE:              grib_ltype = GRIB1_LTYPE_ALTITUDE;           break;
    case ZAXIS_SIGMA:                 grib_ltype = GRIB1_LTYPE_SIGMA;              break;
    case ZAXIS_MEANSEA:               grib_ltype = GRIB1_LTYPE_MEANSEA;            break;
    case ZAXIS_TOA:                   grib_ltype = GRIB1_LTYPE_TOA;                break;
    case ZAXIS_SEA_BOTTOM:            grib_ltype = GRIB1_LTYPE_SEA_BOTTOM;         break;
    case ZAXIS_ATMOSPHERE:            grib_ltype = GRIB1_LTYPE_ATMOSPHERE;         break;
    case ZAXIS_CLOUD_BASE:            grib_ltype = GRIB1_LTYPE_CLOUD_BASE;         break;
    case ZAXIS_CLOUD_TOP:             grib_ltype = GRIB1_LTYPE_CLOUD_TOP;          break;
    case ZAXIS_ISOTHERM_ZERO:         grib_ltype = GRIB1_LTYPE_ISOTHERM0;          break;
    case ZAXIS_SNOW:                  grib_ltype = -1;                             break;
    case ZAXIS_LAKE_BOTTOM:           grib_ltype = GRIB1_LTYPE_LAKE_BOTTOM;        break;
    case ZAXIS_SEDIMENT_BOTTOM:       grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM;    break;
    case ZAXIS_SEDIMENT_BOTTOM_TA:    grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM_TA; break;
    case ZAXIS_SEDIMENT_BOTTOM_TW:    grib_ltype = GRIB1_LTYPE_SEDIMENT_BOTTOM_TW; break;
    case ZAXIS_MIX_LAYER:             grib_ltype = GRIB1_LTYPE_MIX_LAYER;          break;
    case ZAXIS_REFERENCE:             grib_ltype = -1;                             break;
    }

  return grib_ltype;
}


int zaxisTypeToGrib2ltype(int zaxistype)
{
  int grib_ltype = -1;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:               grib_ltype = GRIB2_LTYPE_SURFACE;            break;
    case ZAXIS_GENERIC:               grib_ltype = -1;                             break;
    case ZAXIS_HYBRID:                grib_ltype = GRIB2_LTYPE_HYBRID;             break;
    case ZAXIS_HYBRID_HALF:           grib_ltype = GRIB2_LTYPE_HYBRID;             break;
    case ZAXIS_PRESSURE:              grib_ltype = GRIB2_LTYPE_ISOBARIC;           break;
    case ZAXIS_HEIGHT:                grib_ltype = GRIB2_LTYPE_HEIGHT;             break;
    case ZAXIS_DEPTH_BELOW_SEA:       grib_ltype = GRIB2_LTYPE_SEADEPTH;           break;
    case ZAXIS_DEPTH_BELOW_LAND:      grib_ltype = GRIB2_LTYPE_LANDDEPTH;          break;
    case ZAXIS_ISENTROPIC:            grib_ltype = GRIB2_LTYPE_ISENTROPIC;         break;
    case ZAXIS_TRAJECTORY:            grib_ltype = -1;                             break;
    case ZAXIS_ALTITUDE:              grib_ltype = GRIB2_LTYPE_ALTITUDE;           break;
    case ZAXIS_SIGMA:                 grib_ltype = GRIB2_LTYPE_SIGMA;              break;
    case ZAXIS_MEANSEA:               grib_ltype = GRIB2_LTYPE_MEANSEA;            break;
    case ZAXIS_TOA:                   grib_ltype = GRIB2_LTYPE_TOA;                break;
    case ZAXIS_SEA_BOTTOM:            grib_ltype = GRIB2_LTYPE_SEA_BOTTOM;         break;
    case ZAXIS_ATMOSPHERE:            grib_ltype = GRIB2_LTYPE_ATMOSPHERE;         break;
    case ZAXIS_CLOUD_BASE:            grib_ltype = GRIB2_LTYPE_CLOUD_BASE;         break;
    case ZAXIS_CLOUD_TOP:             grib_ltype = GRIB2_LTYPE_CLOUD_TOP;          break;
    case ZAXIS_ISOTHERM_ZERO:         grib_ltype = GRIB2_LTYPE_ISOTHERM0;          break;
    case ZAXIS_SNOW:                  grib_ltype = GRIB2_LTYPE_SNOW;               break;
    case ZAXIS_LAKE_BOTTOM:           grib_ltype = GRIB2_LTYPE_LAKE_BOTTOM;        break;
    case ZAXIS_SEDIMENT_BOTTOM:       grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM;    break;
    case ZAXIS_SEDIMENT_BOTTOM_TA:    grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM_TA; break;
    case ZAXIS_SEDIMENT_BOTTOM_TW:    grib_ltype = GRIB2_LTYPE_SEDIMENT_BOTTOM_TW; break;
    case ZAXIS_MIX_LAYER:             grib_ltype = GRIB2_LTYPE_MIX_LAYER;          break;
    case ZAXIS_REFERENCE:             grib_ltype = GRIB2_LTYPE_REFERENCE;          break;
    }

  return grib_ltype;
}


int grbBitsPerValue(int datatype)
{
  int bitsPerValue = 16;

  if ( datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64 )
    Error("CDI/GRIB library does not support complex numbers!");

  if ( datatype != CDI_UNDEFID )
    {
      if ( datatype > 0 && datatype <= 32 )
	bitsPerValue = datatype;
      else if ( datatype == CDI_DATATYPE_FLT64 )
	bitsPerValue = 24;
      else
	bitsPerValue = 16;
    }

  return bitsPerValue;
}


/*
int grbInqRecord(stream_t * streamptr, int *varID, int *levelID)
{
  int status;

  status = cgribexInqRecord(streamptr, varID, levelID);

  return (status);
}
*/

void grbDefRecord(stream_t * streamptr)
{
  UNUSED(streamptr);
}

static
int grbScanTimestep1(stream_t * streamptr)
{
  int status = CDI_EUFTYPE;

#if  defined  (HAVE_LIBCGRIBEX)
  int filetype  = streamptr->filetype;

  if ( filetype == CDI_FILETYPE_GRB )
    status = cgribexScanTimestep1(streamptr);
#endif
#if defined(HAVE_LIBCGRIBEX) && defined (HAVE_LIBGRIB_API)
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep1(streamptr);
#endif

  return status;
}

static
int grbScanTimestep2(stream_t * streamptr)
{
  int status = CDI_EUFTYPE;

#if  defined  (HAVE_LIBCGRIBEX)
  int filetype = streamptr->filetype;

  if ( filetype == CDI_FILETYPE_GRB )
    {
      status = cgribexScanTimestep2(streamptr);
    }
#endif
#if defined(HAVE_LIBCGRIBEX) && defined (HAVE_LIBGRIB_API)
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep2(streamptr);
#endif

  return status;
}

static
int grbScanTimestep(stream_t * streamptr)
{
  int status = CDI_EUFTYPE;
  int filetype  = streamptr->filetype;

#if  defined  (HAVE_LIBCGRIBEX)
  if ( filetype == CDI_FILETYPE_GRB )
    status = cgribexScanTimestep(streamptr);
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep(streamptr);
#else
    Error("Sufficient GRIB support unavailable!");
#endif

  return status;
}


#if  defined  (HAVE_LIBGRIB)
int grbInqContents(stream_t * streamptr)
{
  streamptr->curTsID = 0;

  int status = grbScanTimestep1(streamptr);
  if ( status == 0 && streamptr->ntsteps == -1 ) status = grbScanTimestep2(streamptr);

  int fileID = streamptr->fileID;
  fileSetPos(fileID, 0, SEEK_SET);

  return status;
}
#endif

int grbInqTimestep(stream_t * streamptr, int tsID)
{
  if ( tsID == 0 && streamptr->rtsteps == 0 )
    Error("Call to cdiInqContents missing!");

  if ( CDI_Debug )
    Message("tsid = %d rtsteps = %d", tsID, streamptr->rtsteps);

  int ntsteps = CDI_UNDEFID;
  while ( (tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID )
    {
      ntsteps = grbScanTimestep(streamptr);
      if ( ntsteps == CDI_EUFSTRUCT )
	{
	  streamptr->ntsteps = streamptr->rtsteps;
	  break;
	}
    }

  int nrecs;

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}


void streamInqGRIBinfo(int streamID, int *intnum, float *fltnum, off_t *bignum)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  int filetype = streamptr->filetype;

  if ( filetype == CDI_FILETYPE_GRB )
    {
      int tsID     = streamptr->curTsID;
      int vrecID   = streamptr->tsteps[tsID].curRecID;
      int recID    = streamptr->tsteps[tsID].recIDs[vrecID];
      off_t recpos = streamptr->tsteps[tsID].records[recID].position;
      int zip      = streamptr->tsteps[tsID].records[recID].zip;

      void *gribbuffer = streamptr->record->buffer;
      size_t gribbuffersize = streamptr->record->buffersize;

      if ( zip > 0 )
	Error("Compressed GRIB records unsupported!");
      else
        grib_info_for_grads(recpos, (long)gribbuffersize, (unsigned char *) gribbuffer, intnum, fltnum, bignum);
    }
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
