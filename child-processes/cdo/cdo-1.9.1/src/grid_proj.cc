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

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(HAVE_LIBPROJ)
#include "proj_api.h"
#endif

#include <stdio.h>
#include <stdarg.h> /* va_list */

#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"
#include "grid_proj.h"



static
char *gen_param(const char *fmt, ...)
{
  va_list args;
  char str[256];

  va_start(args, fmt);

  int len = vsprintf(str, fmt, args);

  va_end(args);

  len++;
  char *rstr = (char*) Malloc(len*sizeof(char));
  memcpy(rstr, str, len*sizeof(char));

  return rstr;
}

static
void verify_lcc_parameter(double lon_0, double lat_0, double lat_1, double lat_2, double a, double rf, double x_0, double y_0)
{
  const char *projection = "lambert_conformal_conic";

  if ( IS_NOT_EQUAL(a, grid_missval)  && a > 1.e10 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "earth_radius");
  if ( IS_NOT_EQUAL(rf, grid_missval) && rf > 0 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "inverse_flattening");
  if ( lon_0 < -360 || lon_0 > 360 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "longitude_of_central_meridian");
  if ( lat_0 < -90 || lat_0 > 90 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "latitude_of_central_meridian");
  if ( lat_1 < -90 || lat_1 > 90 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "standard_parallel");
  if ( lat_2 < -90 || lat_2 > 90 ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "standard_parallel");
  if ( IS_NOT_EQUAL(x_0, grid_missval) && (x_0 < -1.e20 || x_0 > 1.e20) ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "false_easting");
  if ( IS_NOT_EQUAL(y_0, grid_missval) && (y_0 < -1.e20 || y_0 > 1.e20) ) cdoWarning("%s mapping parameter %s out of bounds!", projection, "false_northing");
}


int proj_lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, size_t nvals, double *xvals, double *yvals)
{
  int status = 0;
#if defined(HAVE_LIBPROJ)
  char *params[20];

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( IS_NOT_EQUAL(a, missval) && a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  if ( IS_NOT_EQUAL(rf, missval) && rf > 0 ) params[nbpar++] = gen_param("rf=%g", rf);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);
  params[nbpar++] = gen_param("units=m");
  //  params[nbpar++] = gen_param("no_defs");

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) status = -1;

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  if ( status == 0 )
    {
      projUV p;
      for ( size_t i = 0; i < nvals; i++ )
        {
          p.u = xvals[i]*DEG_TO_RAD;
          p.v = yvals[i]*DEG_TO_RAD;
          p = pj_fwd(p, proj);
          xvals[i] = p.u;
          yvals[i] = p.v;
        }

      pj_free(proj);
    }
#else
  status = -1;
#endif

  if ( status == -1 )
    for ( size_t i = 0; i < nvals; i++ )
      {
        xvals[i] = missval;
        yvals[i] = missval;
      }

  return status;
}

static
void lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                   double a, double rf, size_t nvals, double *xvals, double *yvals)
{
  int status = proj_lonlat_to_lcc(missval, lon_0, lat_0, lat_1, lat_2, a, rf, nvals, xvals, yvals);
#if defined(HAVE_LIBPROJ)
  if ( status == -1 ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));
#else
  if ( status == -1 ) cdoAbort("proj4 support not compiled in!");
#endif
}


int cdo_lonlat_to_lcc(int gridID, size_t nvals, double *xvals, double *yvals)
{
  double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;
  gridInqParamLCC(gridID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);

  lonlat_to_lcc(grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, nvals, xvals, yvals);

  return 0;
}


int proj_lcc_to_lonlat(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, double x_0, double y_0, size_t nvals, double *xvals, double *yvals)
{
  int status = 0;
#if defined(HAVE_LIBPROJ)
  char *params[20];

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( IS_NOT_EQUAL(a, grid_missval)  && a  > 0 ) params[nbpar++] = gen_param("a=%g", a);
  if ( IS_NOT_EQUAL(rf, grid_missval) && rf > 0 ) params[nbpar++] = gen_param("rf=%g", rf);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);
  if ( IS_NOT_EQUAL(x_0, grid_missval) ) params[nbpar++] = gen_param("x_0=%g", x_0);
  if ( IS_NOT_EQUAL(y_0, grid_missval) ) params[nbpar++] = gen_param("y_0=%g", y_0);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) status = -1;

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */
  
  if ( status == 0 )
    {
      projUV p;
      for ( size_t i = 0; i < nvals; i++ )
        {
          p.u = xvals[i];
          p.v = yvals[i];
          p = pj_inv(p, proj);
          xvals[i] = p.u*RAD_TO_DEG;
          yvals[i] = p.v*RAD_TO_DEG;
        }

      pj_free(proj);
    }
#else
  status = -1;
#endif

  if ( status == -1 )
    for ( size_t i = 0; i < nvals; i++ )
      {
        xvals[i] = missval;
        yvals[i] = missval;
      }

  return status;
}

static
void lcc_to_lonlat(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                   double a, double rf, double x_0, double y_0, size_t nvals, double *xvals, double *yvals)
{
  int status = proj_lcc_to_lonlat(missval, lon_0, lat_0, lat_1, lat_2, a, rf, x_0, y_0, nvals, xvals, yvals);
#if defined(HAVE_LIBPROJ)
  if ( status == -1 ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));
#else
  if ( status == -1 ) cdoAbort("proj4 support not compiled in!");
#endif
}


int cdo_lcc_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals)
{
  const char *projection = "lambert_conformal_conic";

  double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;
  gridInqParamLCC(gridID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);

  int status = 0;
  if ( !status && IS_EQUAL(lon_0, grid_missval) ) { status = 1; cdoWarning("%s mapping parameter %s missing!", projection, "longitude_of_central_meridian"); }
  if ( !status && IS_EQUAL(lat_0, grid_missval) ) { status = 1; cdoWarning("%s mapping parameter %s missing!", projection, "latitude_of_central_meridian"); }
  if ( !status && IS_EQUAL(lat_1, grid_missval) ) { status = 1; cdoWarning("%s mapping parameter %s missing!", projection, "standard_parallel"); }
  if ( !status && IS_EQUAL(x_0, grid_missval) && IS_EQUAL(y_0, grid_missval) && IS_NOT_EQUAL(xval_0, grid_missval) && IS_NOT_EQUAL(yval_0, grid_missval) )
    {
#if defined(HAVE_LIBPROJ)
      x_0 = xval_0; y_0 = yval_0;
      lonlat_to_lcc(grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, 1, &x_0, &y_0);
      x_0 = -x_0; y_0 = -y_0;
#else
      status = 1; cdoWarning("%s mapping parameter %s missing!", projection, "false_easting and false_northing");
#endif
    }

  if ( status ) cdoAbort("%s mapping parameter missing!", projection);

  verify_lcc_parameter(lon_0, lat_0, lat_1, lat_2, a, rf, x_0, y_0);

  lcc_to_lonlat(grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, x_0, y_0, nvals, xvals, yvals);

  return 0;
}


void grid_def_param_sinu(int gridID)
{
  const char *projection = "sinusoidal";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, (int)strlen(projection)+1, projection);
  const char *mapvarname = "Sinusoidal";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(mapvarname)+1, mapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)strlen(projection), projection);
}


void grid_def_param_laea(int gridID, double a, double lon_0, double lat_0)
{
  const char *projection = "lambert_azimuthal_equal_area";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, (int)strlen(projection)+1, projection);
  const char *mapvarname = "Lambert_AEA";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(mapvarname)+1, mapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)strlen(projection), projection);
  
  cdiDefAttFlt(gridID, CDI_GLOBAL, "earth_radius", CDI_DATATYPE_FLT64, 1, &a);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "longitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lat_0);
}


void cdo_sinu_to_lonlat(size_t nvals, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=sinu");
  params[nbpar++] = gen_param("ellps=WGS84");

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, params);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  projUV p;
  for ( size_t i = 0; i < nvals; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
bool cdiInqAttConvertedToFloat(int gridID, int atttype, const char *attname, int attlen, double *attflt)
{
  bool status = true;

  if ( atttype == CDI_DATATYPE_INT32 )
    {
      std::vector<int> attint(attlen);
      cdiInqAttInt(gridID, CDI_GLOBAL, attname, attlen, attint.data());
      for ( int i = 0; i < attlen; ++i ) attflt[i] = (double)attint[i];
    }
  else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
    {
      cdiInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, attflt);
    }
  else
    {
      status = false;
    }

  return status;
}

static
void grid_inq_param_laea(int gridID, double *a, double *lon_0, double *lat_0, double *x_0, double *y_0)
{
  *a = 0; *lon_0 = 0; *lat_0 = 0, *x_0 = 0, *y_0 = 0;

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_PROJECTION )
    {
      const char *projection = "lambert_azimuthal_equal_area";
      char mapping[CDI_MAX_NAME]; mapping[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
      if ( mapping[0] && strcmp(mapping, projection) == 0 )
        {
          int atttype, attlen;
          char attname[CDI_MAX_NAME+1];

          int natts;
          cdiInqNatts(gridID, CDI_GLOBAL, &natts);

          for ( int iatt = 0; iatt < natts; ++iatt )
            {
              cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
              if ( attlen != 1 ) continue;

              double attflt;
              if ( cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, &attflt) )
                {
                  if      ( strcmp(attname, "earth_radius") == 0 )                    *a     = attflt;
                  else if ( strcmp(attname, "longitude_of_projection_origin") == 0 )  *lon_0 = attflt;
                  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )   *lat_0 = attflt;
                  else if ( strcmp(attname, "false_easting")  == 0 )                  *x_0   = attflt;
                  else if ( strcmp(attname, "false_northing") == 0 )                  *y_0   = attflt;
                }
            }
        }
      else
        cdoWarning("%s mapping parameter missing!", projection);
    }
}


void cdo_laea_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];
  
  double a, lon_0, lat_0, x_0, y_0;
  grid_inq_param_laea(gridID, &a, &lon_0, &lat_0, &x_0, &y_0);

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=laea");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  if ( IS_NOT_EQUAL(x_0,0) ) params[nbpar++] = gen_param("x_0=%g", x_0);
  if ( IS_NOT_EQUAL(y_0,0) ) params[nbpar++] = gen_param("y_0=%g", y_0);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  projUV p;
  for ( size_t i = 0; i < nvals; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}


void cdo_proj_to_lonlat(char *proj4param, size_t nvals, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[99];

  int nbpar;
  for ( nbpar = 0; nbpar < 99; ++nbpar )
    {
      while ( *proj4param == ' ' || *proj4param == '+' ) proj4param++;
      if ( *proj4param == 0 ) break;
      char *cstart = proj4param;
      while ( *proj4param != ' ' && *proj4param != 0 ) proj4param++;
      char *cend = proj4param;
      size_t len = cend - cstart;
      if ( len <= 0 ) break;
      bool lend = *cend == 0;
      if ( !lend ) *cend = 0;
      params[nbpar] = strdup(cstart);
      if ( !lend ) *cend = ' ';
    }

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  projUV p;
  for ( size_t i = 0; i < nvals; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}
