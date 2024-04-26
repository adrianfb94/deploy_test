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

#include <cdi.h>
#include "cdo_int.h"


#define MAX_LINE_LEN 65536


typedef struct {
  double *vals;
  double *lbounds;
  double *ubounds;
  double *vct;
  size_t  vctsize;
  int     type;
  int     datatype;
  size_t  size;
  bool    scalar;
  char    name[CDI_MAX_NAME];
  char    longname[CDI_MAX_NAME];
  char    units[CDI_MAX_NAME];
}
zaxis_t;


void zaxisInit(zaxis_t *zaxis)
{
  zaxis->vals        = NULL;
  zaxis->lbounds     = NULL;
  zaxis->ubounds     = NULL;
  zaxis->vct         = NULL;
  zaxis->type        = CDI_UNDEFID;
  zaxis->datatype    = CDI_UNDEFID;
  zaxis->vctsize     = 0;
  zaxis->size        = 0;
  zaxis->scalar      = false;
  zaxis->name[0]     = 0;
  zaxis->longname[0] = 0;
  zaxis->units[0]    = 0;
}

static
int getoptname(char *optname, const char *optstring, int nopt)
{
  int nerr = 0;

  const char *pname = optstring;
  const char *pend  = optstring;

  for ( int i = 0; i < nopt; i++ )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL ) break;
      pname = pend + 1;
    }

  if ( pend )
    {
      pend = strchr(pname, ',');
      size_t namelen = (pend == NULL) ? strlen(pname) : (size_t)(pend-pname);
      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return nerr;
}


int zaxisDefine(zaxis_t zaxis)
{
  if ( zaxis.type == CDI_UNDEFID ) cdoAbort("zaxistype undefined!");
  if ( zaxis.size ==  0 ) cdoAbort("zaxis size undefined!");

  int zaxisID = zaxisCreate(zaxis.type, (int)zaxis.size);

  if ( zaxis.size == 1 && zaxis.scalar ) zaxisDefScalar(zaxisID);

  if ( zaxis.datatype != CDI_UNDEFID ) zaxisDefDatatype(zaxisID, zaxis.datatype);

  if ( zaxis.vals )
    {
      zaxisDefLevels(zaxisID, zaxis.vals);
      Free(zaxis.vals);
    }
  if ( zaxis.lbounds )
    {
      zaxisDefLbounds(zaxisID, zaxis.lbounds);
      Free(zaxis.lbounds);
    }
  if ( zaxis.ubounds )
    {
      zaxisDefUbounds(zaxisID, zaxis.ubounds);
      Free(zaxis.ubounds);
    }

  if ( zaxis.name[0] )     zaxisDefName(zaxisID, zaxis.name);
  if ( zaxis.longname[0] ) zaxisDefLongname(zaxisID, zaxis.longname);
  if ( zaxis.units[0] )    zaxisDefUnits(zaxisID, zaxis.units);

  if ( zaxis.type == ZAXIS_HYBRID || zaxis.type == ZAXIS_HYBRID_HALF )
    {
      if ( zaxis.vctsize && zaxis.vct )
	zaxisDefVct(zaxisID, (int)zaxis.vctsize, zaxis.vct);
      else
	cdoWarning("vct undefined!");	    
    }

  return zaxisID;
}


typedef struct {
  keyValues_t *kv;
  bool isValid;
} kvmap_t;

static
void zaxis_read_data(size_t nkv, kvmap_t *kvmap, zaxis_t *zaxis, size_t *iatt, const char *dname)
{
  // char uuidStr[256];

  for ( size_t ik = 0; ik < nkv; ++ik )
    {
      if ( !kvmap[ik].isValid ) continue;

      keyValues_t *kv = kvmap[ik].kv;
      const char *key = kv->key;
      // size_t nvalues = kv->nvalues;
      const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;
      // bool lv1 = (kv->nvalues == 1);

      // printf("%s = ", key); if ( value  ) printf("%s", kv->value); printf("\n");

      if ( STR_IS_EQ(key, "zaxistype") )
        {
          const char *zaxistype = parameter2word(value);

          if      ( STR_IS_EQ(zaxistype, "pressure") )          zaxis->type = ZAXIS_PRESSURE;
          else if ( STR_IS_EQ(zaxistype, "hybrid_half") )       zaxis->type = ZAXIS_HYBRID_HALF;
          else if ( STR_IS_EQ(zaxistype, "hybrid") )            zaxis->type = ZAXIS_HYBRID;
          else if ( STR_IS_EQ(zaxistype, "height") )            zaxis->type = ZAXIS_HEIGHT;
          else if ( STR_IS_EQ(zaxistype, "depth_below_sea") )   zaxis->type = ZAXIS_DEPTH_BELOW_SEA;
          else if ( STR_IS_EQ(zaxistype, "depth_below_land") )  zaxis->type = ZAXIS_DEPTH_BELOW_LAND;
          else if ( STR_IS_EQ(zaxistype, "isentropic") )        zaxis->type = ZAXIS_ISENTROPIC;
          else if ( STR_IS_EQ(zaxistype, "surface") )           zaxis->type = ZAXIS_SURFACE;
          else if ( STR_IS_EQ(zaxistype, "generic") )           zaxis->type = ZAXIS_GENERIC;
	  else cdoAbort("Invalid zaxisname : %s (zaxis description file: %s)", zaxistype, dname);
        }
      else if ( STR_IS_EQ(key, "datatype") )
        {
          const char *datatype = parameter2word(value);

          if      ( STR_IS_EQ(datatype, "double") )  zaxis->datatype = CDI_DATATYPE_FLT64;
          else if ( STR_IS_EQ(datatype, "float") )   zaxis->datatype = CDI_DATATYPE_FLT32;
	  else cdoAbort("Invalid datatype : %s (zaxis description file: %s)", datatype, dname);
        }
      else if ( STR_IS_EQ(key, "size") ) zaxis->size = parameter2int(value);
      else if ( STR_IS_EQ(key, "scalar") ) zaxis->scalar = parameter2bool(value);
      else if ( STR_IS_EQ(key, "vctsize") ) zaxis->vctsize = parameter2int(value);
      else if ( STR_IS_EQ(key, "name") ) strcpy(zaxis->name, parameter2word(value));
      else if ( STR_IS_EQ(key, "units") ) strcpy(zaxis->units, parameter2word(value));
      else if ( STR_IS_EQ(key, "longname") ) strcpy(zaxis->longname, value);
      else if ( STR_IS_EQ(key, "levels") )
        {
          if ( zaxis->size == 0 ) cdoAbort("size undefined (zaxis description file: %s)!", dname);
          zaxis->vals = (double*) Malloc(zaxis->size*sizeof(double));
          for ( size_t i = 0; i < zaxis->size; ++i ) zaxis->vals[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "lbounds") )
        {
          if ( zaxis->size == 0 ) cdoAbort("size undefined (zaxis description file: %s)!", dname);
          zaxis->lbounds = (double*) Malloc(zaxis->size*sizeof(double));
          for ( size_t i = 0; i < zaxis->size; ++i ) zaxis->lbounds[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "ubounds") )
        {
          if ( zaxis->size == 0 ) cdoAbort("size undefined (zaxis description file: %s)!", dname);
          zaxis->ubounds = (double*) Malloc(zaxis->size*sizeof(double));
          for ( size_t i = 0; i < zaxis->size; ++i ) zaxis->ubounds[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "vct") )
        {
          if ( zaxis->vctsize == 0 ) cdoAbort("vctsize undefined (zaxis description file: %s)!", dname);
          zaxis->vct = (double*) Malloc(zaxis->vctsize*sizeof(double));
          for ( size_t i = 0; i < zaxis->vctsize; ++i ) zaxis->vct[i] = parameter2double(kv->values[i]);
        }
      else
        {
          *iatt = ik; break;
        }
    }
}

static
void zaxis_read_attributes(size_t iatt, size_t nkv, kvmap_t *kvmap, int zaxisID)
{
  const char *reserved_keys[] = {"zaxistype", "size", "scalar", "vctsize", "name", "units", "longname", "levels", "lbounds", "ubounds", "vct"};
  int num_rkeys = sizeof(reserved_keys) / sizeof(char*);
  const char *attkey0 = NULL;

  for ( size_t ik = iatt; ik < nkv; ++ik )
    {
      if ( !kvmap[ik].isValid ) continue;

      keyValues_t *kv = kvmap[ik].kv;
      const char *key = kv->key;
      size_t nvalues = kv->nvalues;
      const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;

      if ( ik == iatt ) attkey0 = key;
      else
        {
          for ( int n = 0; n < num_rkeys; ++n )
            if ( STR_IS_EQ(key, reserved_keys[n]) )
              cdoAbort("Found reserved key word >%s< in attribute names! Check name or position of >%s<.", key, attkey0);
        }
      
      int dtype = literals_find_datatype(nvalues, kv->values);

      if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
        {
          int *ivals = (int*) Malloc(nvalues*sizeof(int));
          for ( size_t i = 0; i < nvalues; ++i ) ivals[i] = literal_to_int(kv->values[i]);
          cdiDefAttInt(zaxisID, CDI_GLOBAL, key, dtype, nvalues, ivals);
          Free(ivals);
        }
      else if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
        {
          double *dvals = (double*) Malloc(nvalues*sizeof(double));
          for ( size_t i = 0; i < nvalues; ++i ) dvals[i] = literal_to_double(kv->values[i]);
          cdiDefAttFlt(zaxisID, CDI_GLOBAL, key, dtype, nvalues, dvals);
          Free(dvals);
        }
      else
        {
          int len = (value && *value) ? (int) strlen(value) : 0;
          cdiDefAttTxt(zaxisID, CDI_GLOBAL, key, len, value);
        }
    }
}


int zaxisFromFile(FILE *gfp, const char *dname)
{
  list_t *pmlist = namelist_to_pmlist(gfp, dname);
  if ( pmlist == NULL ) return -1;
  list_t *kvlist = *(list_t **)pmlist->head->data;
  if ( kvlist == NULL ) return -1;

  size_t nkv = list_size(kvlist);
  if ( nkv == 0 ) return -1;
  kvmap_t *kvmap = (kvmap_t*) Malloc(nkv*sizeof(kvmap_t));
  for ( size_t i = 0; i < nkv; ++i ) kvmap[i].isValid = false;

  size_t ik = 0;
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      if ( ik == 0 && !STR_IS_EQ(kv->key, "zaxistype") )
        cdoAbort("First zaxis description key word must be >zaxistype< (found: %s)!", kv->key);

      if ( kv->nvalues == 0 )
        {
          cdoWarning("Z-axis description key word %s has no values, skipped!", kv->key);
        }
      else
        {
          kvmap[ik].isValid = true;
          kvmap[ik].kv = kv;
        }
      ik++;
    }

  zaxis_t zaxis;
  zaxisInit(&zaxis);

  size_t iatt = 0;
  zaxis_read_data(nkv, kvmap, &zaxis, &iatt, dname);

  int zaxisID = (zaxis.type == CDI_UNDEFID) ? CDI_UNDEFID : zaxisDefine(zaxis);
  if ( zaxisID != CDI_UNDEFID && iatt > 0 ) zaxis_read_attributes(iatt, nkv, kvmap, zaxisID);

  list_destroy(pmlist);

  Free(kvmap);

  return zaxisID;
}

static
void gen_zaxis_height(zaxis_t *zaxis, const char *pline)
{
  int zaxistype = ZAXIS_HEIGHT;

  if ( CDO_CMOR_Mode ) zaxis->scalar = true;

  if ( *pline != 0 )
    {
      if ( *pline == '_' ) pline++;
      else return;

      if ( *pline == 0 ) return;

      if ( ! isdigit((int) *pline) && !ispunct((int) *pline) ) return;

      char *endptr = (char *) pline;
      double value = strtod(pline, &endptr);
      if ( *endptr != 0 )
        {
          pline = endptr;
          if ( *pline == '_' ) pline++;
          
          if ( *pline == 0 ) return;
          const char *units = pline;

          zaxis->type = zaxistype;
          zaxis->size = 1;
          // zaxis->scalar = true;
          double *levels = (double*) Malloc(sizeof(double));
          *levels = value;
          zaxis->vals = levels;
          strcpy(zaxis->units, units);

          size_t len = strlen(units);
          if ( len > 2 && units[len-2] == '_' && units[len-1] == 's' )
            {
              zaxis->units[len-2] = 0;
              zaxis->scalar = true;
            }
        }
    }
}


int zaxisFromName(const char *zaxisnameptr)
{
  int zaxisID = CDI_UNDEFID;
  size_t len;

  char *zaxisname = strdup(zaxisnameptr);
  strtolower(zaxisname);

  zaxis_t zaxis;
  zaxisInit(&zaxis);

  const char *pline = zaxisname;
  if ( cmpstr(pline, "surface") == 0 ) /* surface */
    {
      zaxis.type = ZAXIS_SURFACE;
      zaxis.size = 1;
      zaxis.vals = (double*) Malloc(zaxis.size*sizeof(double));
      zaxis.vals[0] = 0;
    }
  else if ( cmpstrlen(zaxisname, "height", len) == 0 )
    {
      pline = &zaxisname[len];
      gen_zaxis_height(&zaxis, pline);
    }

  if ( zaxis.type != CDI_UNDEFID ) zaxisID = zaxisDefine(zaxis);

  free(zaxisname);

  return zaxisID;
}


int cdoDefineZaxis(const char *zaxisfile)
{
  int zaxisID = CDI_UNDEFID;

  FILE *zfp = fopen(zaxisfile, "r");
  if ( zfp == NULL )
    {
      zaxisID = zaxisFromName(zaxisfile);

      if ( zaxisID == CDI_UNDEFID ) cdoAbort("Open failed on %s!", zaxisfile);
    }
  else
    {
      zaxisID = zaxisFromFile(zfp, zaxisfile);
      fclose(zfp);
    }

  if ( zaxisID == CDI_UNDEFID ) cdoAbort("Invalid zaxis description file %s!", zaxisfile);

  return zaxisID;
}


void defineZaxis(const char *zaxisarg)
{
  char zaxisfile[4096];
  int nfile = 0;

  while ( getoptname(zaxisfile, zaxisarg, nfile++) == 0 )
    {      
      (void) cdoDefineZaxis(zaxisfile);
    }
}

static
int ztype2ltype(int zaxistype)
{
  int ltype = CDI_UNDEFID;

  // clang-format off
  if      ( zaxistype == ZAXIS_SURFACE           )  ltype =   1;
  else if ( zaxistype == ZAXIS_PRESSURE          )  ltype = 100;
  else if ( zaxistype == ZAXIS_ALTITUDE          )  ltype = 103;
  else if ( zaxistype == ZAXIS_HEIGHT            )  ltype = 105;
  else if ( zaxistype == ZAXIS_SIGMA             )  ltype = 107;
  else if ( zaxistype == ZAXIS_HYBRID            )  ltype = 109;
  else if ( zaxistype == ZAXIS_HYBRID_HALF       )  ltype = 109;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_LAND  )  ltype = 111;
  else if ( zaxistype == ZAXIS_ISENTROPIC        )  ltype = 113;
  else if ( zaxistype == ZAXIS_DEPTH_BELOW_SEA   )  ltype = 160;
  // clang-format on

  return ltype;
}


int zaxis2ltype(int zaxisID)
{
  int ltype = zaxisInqLtype(zaxisID);
  if ( ltype <= 0 )
    {
      int zaxistype = zaxisInqType(zaxisID);
      ltype = ztype2ltype(zaxistype);
    }

  return ltype;
}
