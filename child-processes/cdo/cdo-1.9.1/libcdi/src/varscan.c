#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif


#include "cdi.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "dmemory.h"
#include "resource_handle.h"
#include "varscan.h"
#include "vlist.h"
#include "grid.h"
#include "zaxis.h"
#include "subtype.h"


static size_t Vctsize = 0;
static double *Vct = NULL;

static int numberOfVerticalLevels = 0;
static int numberOfVerticalGrid = 0;
static unsigned char uuidVGrid[CDI_UUID_SIZE];


typedef struct
{
  int      level1;
  int      level2;
  int      recID;
  int      lindex;
}
leveltable_t;


typedef struct
{
  int           subtypeIndex; /* corresponding tile in subtype_t structure (subtype->self) */

  unsigned      nlevels;
  int           levelTableSize;
  leveltable_t* levelTable;
} subtypetable_t;


typedef struct
{
  int            varID;
  int            param;
  int            prec;
  int            tsteptype;
  int            timaccu;
  int            gridID;
  int            zaxistype;
  int            ltype1;     /* GRIB first level type */
  int            ltype2;     /* GRIB second level type */
  int            lbounds;
  int            level_sf;
  int            level_unit;
  int            zaxisID;

  int            nsubtypes_alloc;
  int            nsubtypes;
  subtypetable_t *recordTable;    /* ~ two-dimensional record list [nsubtypes_alloc][levelTableSize] */

  int            instID;
  int            modelID;
  int            tableID;
  int            comptype;       // compression type
  int            complevel;      // compression level
  short          timave;
  bool           lmissval;
  double         missval;
  char          *name;
  char          *stdname;
  char          *longname;
  char          *units;
  ensinfo_t     *ensdata;
  int            typeOfGeneratingProcess;
  int            productDefinitionTemplate;

  /* meta-data for specification of tiles (currently only GRIB-API: */
  subtype_t     *tiles;

  int                 opt_grib_nentries;                       /* current no. key-value pairs            */
  int                 opt_grib_kvpair_size;                    /* current allocated size                 */
  opt_key_val_pair_t *opt_grib_kvpair;                         /* (optional) list of keyword/value pairs */
}
vartable_t;


static vartable_t *vartable;
static unsigned varTablesize = 0;
static unsigned varTableUsed = 0;

static
void paramInitEntry(unsigned varID, int param)
{
  vartable[varID].varID          = (int)varID;
  vartable[varID].param          = param;
  vartable[varID].prec           = 0;
  vartable[varID].tsteptype      = TSTEP_INSTANT;
  vartable[varID].timave         = 0;
  vartable[varID].timaccu        = 0;
  vartable[varID].gridID         = CDI_UNDEFID;
  vartable[varID].zaxistype      = 0;
  vartable[varID].ltype1         = 0;
  vartable[varID].ltype2         = -1;
  vartable[varID].lbounds        = 0;
  vartable[varID].level_sf       = 0;
  vartable[varID].level_unit     = 0;
  vartable[varID].recordTable    = NULL;
  vartable[varID].nsubtypes_alloc= 0;
  vartable[varID].nsubtypes      = 0;
  vartable[varID].instID         = CDI_UNDEFID;
  vartable[varID].modelID        = CDI_UNDEFID;
  vartable[varID].tableID        = CDI_UNDEFID;
  vartable[varID].typeOfGeneratingProcess   = CDI_UNDEFID;
  vartable[varID].productDefinitionTemplate = CDI_UNDEFID;
  vartable[varID].comptype       = CDI_COMPRESS_NONE;
  vartable[varID].complevel      = 1;
  vartable[varID].lmissval       = false;
  vartable[varID].missval        = 0;
  vartable[varID].name           = NULL;
  vartable[varID].stdname        = NULL;
  vartable[varID].longname       = NULL;
  vartable[varID].units          = NULL;
  vartable[varID].ensdata        = NULL;
  vartable[varID].tiles          = NULL;
}

/* Test if a variable specified by the given meta-data has already
 * been registered in "vartable". */
static unsigned
varGetEntry(int param, int zaxistype, int ltype1, int tsteptype, const char *name, const var_tile_t *tiles)
{
  for ( unsigned varID = 0; varID < varTablesize; varID++ )
    {
      /* testing for "param" implicitly checks if we are beyond the
       * current vartable size: */
      if ( vartable[varID].param == param )
        {
          int no_of_tiles = -1;
          if ( tiles ) no_of_tiles = tiles->numberOfTiles;
          int vt_no_of_tiles = -1;
          if ( vartable[varID].tiles )
            vt_no_of_tiles = subtypeGetGlobalDataP(vartable[varID].tiles,
                                                   SUBTYPE_ATT_NUMBER_OF_TILES);
          if ( (vartable[varID].zaxistype  == zaxistype)               &&
               (vartable[varID].ltype1     == ltype1   )               &&
               (vartable[varID].tsteptype  == tsteptype)               &&
               (vt_no_of_tiles == no_of_tiles) )
            {
              if ( name && name[0] && vartable[varID].name && vartable[varID].name[0] )
                {
                  if ( strcmp(name, vartable[varID].name) == 0 ) return varID;
                }
              else
                {
                  return varID;
                }
            }
        }
    }

  return (unsigned)-1;
}

static
void varFree(void)
{
  if ( CDI_Debug ) Message("call to varFree");

  for ( size_t varID = 0; varID < varTableUsed; varID++ )
    {
      if ( vartable[varID].recordTable )
        {
          for (int isub=0; isub<vartable[varID].nsubtypes_alloc; isub++)
            Free(vartable[varID].recordTable[isub].levelTable);
          Free(vartable[varID].recordTable);
        }

      if ( vartable[varID].name )     Free(vartable[varID].name);
      if ( vartable[varID].stdname )  Free(vartable[varID].stdname);
      if ( vartable[varID].longname ) Free(vartable[varID].longname);
      if ( vartable[varID].units )    Free(vartable[varID].units);
      if ( vartable[varID].ensdata )  Free(vartable[varID].ensdata);
      if ( vartable[varID].tiles )    subtypeDestroyPtr(vartable[varID].tiles);

      if ( vartable[varID].opt_grib_kvpair )
        {
          for (int i=0; i<vartable[varID].opt_grib_nentries; i++) {
            if ( vartable[varID].opt_grib_kvpair[i].keyword )
              Free(vartable[varID].opt_grib_kvpair[i].keyword);
          }
          Free(vartable[varID].opt_grib_kvpair);
        }
      vartable[varID].opt_grib_nentries    = 0;
      vartable[varID].opt_grib_kvpair_size = 0;
      vartable[varID].opt_grib_kvpair      = NULL;
    }

  if ( vartable )
    Free(vartable);

  vartable = NULL;
  varTablesize = 0;
  varTableUsed = 0;

  if ( Vct )
    Free(Vct);

  Vct = NULL;
  Vctsize = 0;
}

/* Search for a tile subtype with subtypeIndex == tile_index. */
static int tileGetEntry(unsigned varID, int tile_index)
{
  for (int isub=0; isub<vartable[varID].nsubtypes; isub++)
    if (vartable[varID].recordTable[isub].subtypeIndex == tile_index)
      return isub;
  return CDI_UNDEFID;
}


/* Resizes vartable:recordTable data structure, if necessary. */
static int tileNewEntry(int varID)
{
  int tileID = 0;
  if (vartable[varID].nsubtypes_alloc == 0)
    {
      /* create table for the first time. */
      vartable[varID].nsubtypes_alloc = 2;
      vartable[varID].nsubtypes       = 0;
      vartable[varID].recordTable     =
        (subtypetable_t *) Malloc((size_t)vartable[varID].nsubtypes_alloc * sizeof (subtypetable_t));
      if( vartable[varID].recordTable == NULL )
        SysError("Allocation of leveltable failed!");

      for (int isub = 0; isub<vartable[varID].nsubtypes_alloc; isub++) {
	vartable[varID].recordTable[isub].levelTable     = NULL;
        vartable[varID].recordTable[isub].levelTableSize = 0;
        vartable[varID].recordTable[isub].nlevels        = 0;
        vartable[varID].recordTable[isub].subtypeIndex   = CDI_UNDEFID;
      }
    }
  else
    {
      /* data structure large enough; find a free entry. */
      while(tileID <  vartable[varID].nsubtypes_alloc)
	{
	  if (vartable[varID].recordTable[tileID].levelTable == NULL) break;
	  tileID++;
	}
    }

  /* If the table overflows, double its size. */
  if (tileID == vartable[varID].nsubtypes_alloc)
    {
      tileID = vartable[varID].nsubtypes_alloc;
      vartable[varID].nsubtypes_alloc *= 2;
      vartable[varID].recordTable   =
        (subtypetable_t *) Realloc(vartable[varID].recordTable,
                                   (size_t)vartable[varID].nsubtypes_alloc * sizeof (subtypetable_t));
      if (vartable[varID].recordTable == NULL)
        SysError("Reallocation of leveltable failed");
      for(int isub=tileID; isub<vartable[varID].nsubtypes_alloc; isub++) {
	vartable[varID].recordTable[isub].levelTable     = NULL;
        vartable[varID].recordTable[isub].levelTableSize = 0;
        vartable[varID].recordTable[isub].nlevels        = 0;
        vartable[varID].recordTable[isub].subtypeIndex   = CDI_UNDEFID;
      }
    }

  return tileID;
}


static int levelNewEntry(unsigned varID, int level1, int level2, int tileID)
{
  int levelID = 0;
  int levelTableSize = vartable[varID].recordTable[tileID].levelTableSize;
  leveltable_t *levelTable = vartable[varID].recordTable[tileID].levelTable;

  /*
    Look for a free slot in levelTable.
    (Create the table the first time through).
  */
  if ( ! levelTableSize )
    {
      levelTableSize = 2;
      levelTable = (leveltable_t *) Malloc((size_t)levelTableSize
                                           * sizeof (leveltable_t));
      for ( int i = 0; i < levelTableSize; i++ )
        levelTable[i].recID = CDI_UNDEFID;
    }
  else
    {
      while( levelID < levelTableSize
             && levelTable[levelID].recID != CDI_UNDEFID )
        ++levelID;
    }
  /*
    If the table overflows, double its size.
  */
  if ( levelID == levelTableSize )
    {
      levelTable = (leveltable_t *) Realloc(levelTable,
                                            (size_t)(levelTableSize *= 2)
                                            * sizeof (leveltable_t));
      for( int i = levelID; i < levelTableSize; i++ )
        levelTable[i].recID = CDI_UNDEFID;
    }

  levelTable[levelID].level1   = level1;
  levelTable[levelID].level2   = level2;
  levelTable[levelID].lindex   = levelID;

  vartable[varID].recordTable[tileID].nlevels        = (unsigned)levelID+1;
  vartable[varID].recordTable[tileID].levelTableSize = levelTableSize;
  vartable[varID].recordTable[tileID].levelTable     = levelTable;

  return levelID;
}

#define  UNDEF_PARAM  -4711

static unsigned
paramNewEntry(int param)
{
  unsigned varID = 0;

  /*
    Look for a free slot in vartable.
    (Create the table the first time through).
  */
  if ( ! varTablesize )
    {
      varTablesize = 2;
      vartable = (vartable_t *) Malloc((size_t)varTablesize
                                       * sizeof (vartable_t));
      if( vartable == NULL )
	{
          Message("varTablesize = %d", varTablesize);
	  SysError("Allocation of vartable failed");
	}

      for( unsigned i = 0; i < varTablesize; i++ )
        {
          vartable[i].param = UNDEF_PARAM;
          vartable[i].opt_grib_kvpair      = NULL;
          vartable[i].opt_grib_kvpair_size = 0;
          vartable[i].opt_grib_nentries    = 0;
        }
    }
  else
    {
      while( varID < varTablesize )
	{
	  if ( vartable[varID].param == UNDEF_PARAM ) break;
	  varID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == varTablesize )
    {
      vartable = (vartable_t *) Realloc(vartable, (size_t)(varTablesize *= 2)
                                        * sizeof (vartable_t));
      for ( size_t i = varID; i < varTablesize; i++ )
        {
          vartable[i].param = UNDEF_PARAM;
          vartable[i].opt_grib_kvpair      = NULL;
          vartable[i].opt_grib_kvpair_size = 0;
          vartable[i].opt_grib_nentries    = 0;
        }
    }

  paramInitEntry(varID, param);

  return varID;
}


/* Append tile set to a subtype. Return index of the new tile (i.e.
 * the "entry->self" value). */
static
int varInsertTileSubtype(vartable_t *vptr, const var_tile_t *tiles)
{
  if ( tiles == NULL ) return -1;

  /* first, generate a subtype based on the info in "tiles". */

  subtype_t *subtype_ptr;
  subtypeAllocate(&subtype_ptr, SUBTYPE_TILES);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS, tiles->totalno_of_tileattr_pairs);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TILE_CLASSIFICATION      , tiles->tileClassification);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_NUMBER_OF_TILES          , tiles->numberOfTiles);

  /*
   * Here, we create a tile set for comparison that contains only one
   * tile/attribute pair (based on "tiles").
   */
  struct subtype_entry_t *entry = subtypeEntryInsert(subtype_ptr);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_NUMBER_OF_ATTR,            tiles->numberOfAttributes);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEINDEX,                 tiles->tileindex);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEATTRIBUTE,             tiles->attribute);

  if (vptr->tiles == NULL) {
    vptr->tiles = subtype_ptr;
    return 0;
  }
  else {
    tilesetInsertP(vptr->tiles, subtype_ptr);
    subtypeDestroyPtr(subtype_ptr);
    return vptr->tiles->nentries - 1;
  }

  return CDI_UNDEFID;
}


void varAddRecord(int recID, int param, int gridID, int zaxistype, int lbounds,
		  int level1, int level2, int level_sf, int level_unit, int prec,
		  int *pvarID, int *plevelID, int tsteptype, int numavg, int ltype1, int ltype2,
		  const char *name, const char *stdname, const char *longname, const char *units,
                  const var_tile_t *tiles, int *tile_index)
{
  unsigned varID = (cdiSplitLtype105 != 1 || zaxistype != ZAXIS_HEIGHT) ?
    varGetEntry(param, zaxistype, ltype1, tsteptype, name, tiles) : (unsigned) CDI_UNDEFID;

  if ( varID == (unsigned) CDI_UNDEFID )
    {
      varTableUsed++;
      varID = paramNewEntry(param);
      vartable[varID].gridID     = gridID;
      vartable[varID].zaxistype  = zaxistype;
      vartable[varID].ltype1     = ltype1;
      vartable[varID].ltype2     = ltype2;
      vartable[varID].lbounds    = lbounds;
      vartable[varID].level_sf   = level_sf;
      vartable[varID].level_unit = level_unit;
      vartable[varID].tsteptype  = tsteptype;

      if ( numavg ) vartable[varID].timave = 1;

      if ( name && name[0] )         vartable[varID].name     = strdup(name);
      if ( stdname && stdname[0] )   vartable[varID].stdname  = strdup(stdname);
      if ( longname && longname[0] ) vartable[varID].longname = strdup(longname);
      if ( units && units[0] )       vartable[varID].units    = strdup(units);
    }
  else
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( vartable[varID].gridID != gridID )
	{
	  Message("param = %s gridID = %d", paramstr, gridID);
	  Error("horizontal grid must not change for same parameter!");
	}
      if ( vartable[varID].zaxistype != zaxistype )
	{
	  Message("param = %s zaxistype = %d", paramstr, zaxistype);
	  Error("zaxistype must not change for same parameter!");
	}
    }

  if ( prec > vartable[varID].prec ) vartable[varID].prec = prec;

  /* append current tile to tile subtype info. */
  int this_tile = varInsertTileSubtype(&vartable[varID], tiles);
  int tileID = tileGetEntry(varID, this_tile);
  if ( tile_index ) (*tile_index) = this_tile;
  if ( tileID == CDI_UNDEFID )
    {
      tileID = tileNewEntry((int)varID);
      vartable[varID].recordTable[tileID].subtypeIndex = this_tile;
      vartable[varID].nsubtypes++;
    }

  /* append current level to level table info */
  int levelID = levelNewEntry(varID, level1, level2, tileID);
  if (CDI_Debug)
    Message("vartable[%d].recordTable[%d].levelTable[%d].recID = %d; level1,2=%d,%d",
            varID, tileID, levelID, recID, level1, level2);
  vartable[varID].recordTable[tileID].levelTable[levelID].recID = recID;

  *pvarID   = (int) varID;
  *plevelID = levelID;
}


/*
static
int dblcmp(const void *s1, const void *s2)
{
  int cmp = 0;

  if      ( *((double *) s1) < *((double *) s2) ) cmp = -1;
  else if ( *((double *) s1) > *((double *) s2) ) cmp =  1;

  return cmp;
}
*/
static
int cmpLevelTable(const void* s1, const void* s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t*) s1;
  const leveltable_t *y = (const leveltable_t*) s2;
  /*
  printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  */
  if      ( x->level1 < y->level1 ) cmp = -1;
  else if ( x->level1 > y->level1 ) cmp =  1;

  return cmp;
}

static
int cmpLevelTableInv(const void* s1, const void* s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t*) s1;
  const leveltable_t *y = (const leveltable_t*) s2;
  /*
  printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  */
  if      ( x->level1 < y->level1 ) cmp =  1;
  else if ( x->level1 > y->level1 ) cmp = -1;

  return cmp;
}


struct paraminfo
{
  int varid, param, ltype;
};

static
int cmp_param(const void* s1, const void* s2)
{
  const struct paraminfo *x = (const struct paraminfo*) s1;
  const struct paraminfo *y = (const struct paraminfo*) s2;

  int cmp = (( x->param > y->param ) - ( x->param < y->param )) * 2
           + ( x->ltype > y->ltype ) - ( x->ltype < y->ltype );

  return cmp;
}

struct varinfo
{
  int        varid;
  const char *name;
};
/*
static
int cmp_varname(const void *s1, const void *s2)
{
  const struct varinfo *x = (const struct varinfo *)s1,
                       *y = (const struct varinfo *)s2;
  return strcmp(x->name, y->name);
}
*/
static
int cmp_varname(const void *s1, const void *s2)
{
  const vartable_t *x = (const vartable_t *)s1,
                   *y = (const vartable_t *)s2;
  return strcmp(x->name, y->name);
}


void cdi_generate_vars(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;

  int *varids = (int *) Malloc(varTableUsed*sizeof(int));
  for ( size_t varID = 0; varID < varTableUsed; varID++ ) varids[varID] = (int)varID;
  /*
  if ( streamptr->sortparam )
    {
      struct paraminfo *varInfo = (struct paraminfo *) Malloc((size_t)varTableUsed * sizeof(struct paraminfo));

      for ( unsigned varID = 0; varID < varTableUsed; varID++ )
	{
	  varInfo[varID].varid = varids[varID];
	  varInfo[varID].param = vartable[varID].param;
	  varInfo[varID].ltype = vartable[varID].ltype1;
	}
      qsort(varInfo, (size_t)varTableUsed, sizeof(struct paraminfo), cmp_param);
      for ( unsigned varID = 0; varID < varTableUsed; varID++ )
	{
	  varids[varID] = varInfo[varID].varid;
	}
      Free(varInfo);
    }

  if ( streamptr->sortname )
    {
      qsort(vartable, (size_t)varTableUsed, sizeof(vartable_t), cmp_varname);
    }
  */
  for ( size_t index = 0; index < varTableUsed; index++ )
    {
      int varid      = varids[index];

      int gridID     = vartable[varid].gridID;
      int param      = vartable[varid].param;
      int ltype1     = vartable[varid].ltype1;
      int ltype2     = vartable[varid].ltype2;
      int zaxistype  = vartable[varid].zaxistype;
      if ( ltype1 == 0 && zaxistype == ZAXIS_GENERIC && cdiDefaultLeveltype != -1 )
	zaxistype = cdiDefaultLeveltype;
      int lbounds    = vartable[varid].lbounds;
      int prec       = vartable[varid].prec;
      int instID     = vartable[varid].instID;
      int modelID    = vartable[varid].modelID;
      int tableID    = vartable[varid].tableID;
      int tsteptype  = vartable[varid].tsteptype;
      int timave     = vartable[varid].timave;
      int timaccu    = vartable[varid].timaccu;
      int comptype   = vartable[varid].comptype;

      double level_sf = 1;
      if ( vartable[varid].level_sf != 0 ) level_sf = 1./vartable[varid].level_sf;

      /* consistency check: test if all subtypes have the same levels: */
      unsigned nlevels = vartable[varid].recordTable[0].nlevels;
      for ( int isub = 1; isub < vartable[varid].nsubtypes; isub++ ) {
        if ( vartable[varid].recordTable[isub].nlevels != nlevels )
          {
            fprintf(stderr, "var \"%s\": isub = %d / %d :: "
                    "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                    vartable[varid].name, isub, vartable[varid].nsubtypes,
                    nlevels, vartable[varid].recordTable[isub].nlevels);
            Error("zaxis size must not change for same parameter!");
          }

        leveltable_t *t1 = vartable[varid].recordTable[isub-1].levelTable;
        leveltable_t *t2 = vartable[varid].recordTable[isub  ].levelTable;
        for ( unsigned ilev = 0; ilev < nlevels; ilev++ )
          if ((t1[ilev].level1 != t2[ilev].level1)  ||
              (t1[ilev].level2 != t2[ilev].level2)  ||
              (t1[ilev].lindex != t2[ilev].lindex))
            {
              fprintf(stderr, "var \"%s\", varID=%d: isub = %d / %d :: "
                      "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                      vartable[varid].name, varid, isub, vartable[varid].nsubtypes,
                      nlevels, vartable[varid].recordTable[isub].nlevels);
              Message("t1[ilev].level1=%d / t2[ilev].level1=%d", t1[ilev].level1, t2[ilev].level1);
              Message("t1[ilev].level2=%d / t2[ilev].level2=%d", t1[ilev].level2, t2[ilev].level2);
              Message("t1[ilev].lindex=%d / t2[ilev].lindex=%d", t1[ilev].lindex, t2[ilev].lindex);
              Error("zaxis type must not change for same parameter!");
            }
      }
      leveltable_t *levelTable = vartable[varid].recordTable[0].levelTable;

      if ( ltype1 == 0 && zaxistype == ZAXIS_GENERIC && nlevels == 1 && levelTable[0].level1 == 0 )
	zaxistype = ZAXIS_SURFACE;

      double *dlevels = (double *) Malloc(nlevels*sizeof(double));

      if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
	for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = (level_sf*levelTable[levelID].level1 +
	                      level_sf*levelTable[levelID].level2)/2;
      else
	for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
	  dlevels[levelID] = level_sf*levelTable[levelID].level1;

      if ( nlevels > 1 )
	{
          bool linc = true, ldec = true, lsort = false;
          for ( unsigned levelID = 1; levelID < nlevels; levelID++ )
            {
              /* check increasing of levels */
              linc &= (dlevels[levelID] > dlevels[levelID-1]);
              /* check decreasing of levels */
              ldec &= (dlevels[levelID] < dlevels[levelID-1]);
            }
          /*
           * always sort pressure z-axis to ensure
           * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 > levelID2
           * unless already sorted in decreasing order
           */
          if ( (!linc && !ldec) && zaxistype == ZAXIS_PRESSURE )
            {
              qsort(levelTable, nlevels, sizeof(leveltable_t), cmpLevelTableInv);
              lsort = true;
            }
          /*
           * always sort hybrid and depth-below-land z-axis to ensure
           * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 < levelID2
           * unless already sorted in increasing order
           */
          else if ( (!linc && !ldec) ||
                    zaxistype == ZAXIS_HYBRID ||
                    zaxistype == ZAXIS_DEPTH_BELOW_LAND )
            {
              qsort(levelTable, nlevels, sizeof(leveltable_t), cmpLevelTable);
              lsort = true;
            }

          if ( lsort )
            {
              if ( lbounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
                for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
                  dlevels[levelID] = (level_sf*levelTable[levelID].level1 +
                                      level_sf*levelTable[levelID].level2)/2.;
              else
                for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
                  dlevels[levelID] = level_sf*levelTable[levelID].level1;
            }
	}

      double *dlevels1 = NULL;
      double *dlevels2 = NULL;
      if ( lbounds )
	{
	  dlevels1 = (double *) Malloc(nlevels*sizeof(double));
	  for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
	    dlevels1[levelID] = level_sf*levelTable[levelID].level1;
	  dlevels2 = (double *) Malloc(nlevels*sizeof(double));
	  for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
	    dlevels2[levelID] = level_sf*levelTable[levelID].level2;
        }

      const char **cvals = NULL;
      const char *unitptr = cdiUnitNamePtr(vartable[varid].level_unit);
      int zaxisID = varDefZaxis(vlistID, zaxistype, (int)nlevels, dlevels, cvals, 0, lbounds, dlevels1, dlevels2,
                                (int)Vctsize, Vct, NULL, NULL, unitptr, 0, 0, ltype1);

      if ( CDI_cmor_mode && nlevels == 1 && zaxistype != ZAXIS_HYBRID ) zaxisDefScalar(zaxisID);

      if ( ltype1 != ltype2 && ltype2 != -1 ) zaxisDefLtype2(zaxisID, ltype2);

      if ( zaxisInqType(zaxisID) == ZAXIS_REFERENCE )
        {
          if ( numberOfVerticalLevels > 0 ) zaxisDefNlevRef(zaxisID, numberOfVerticalLevels);
          if ( numberOfVerticalGrid > 0 ) zaxisDefNumber(zaxisID, numberOfVerticalGrid);
          if ( !cdiUUIDIsNull(uuidVGrid) ) zaxisDefUUID(zaxisID, uuidVGrid);
        }

      if ( lbounds ) Free(dlevels1);
      if ( lbounds ) Free(dlevels2);
      Free(dlevels);

      /* define new subtype for tile set */
      int tilesetID = CDI_UNDEFID;
      if ( vartable[varid].tiles ) tilesetID = vlistDefTileSubtype(vlistID, vartable[varid].tiles);

      /* generate new variable */
      int varID = stream_new_var(streamptr, gridID, zaxisID, tilesetID);
      varID = vlistDefVarTiles(vlistID, gridID, zaxisID, TIME_VARYING, tilesetID);

      vlistDefVarTsteptype(vlistID, varID, tsteptype);
      vlistDefVarParam(vlistID, varID, param);
      vlistDefVarDatatype(vlistID, varID, prec);
      vlistDefVarTimave(vlistID, varID, timave);
      vlistDefVarTimaccu(vlistID, varID, timaccu);
      vlistDefVarCompType(vlistID, varID, comptype);

      if ( vartable[varid].typeOfGeneratingProcess != CDI_UNDEFID )
        vlistDefVarTypeOfGeneratingProcess(vlistID, varID, vartable[varid].typeOfGeneratingProcess);

      if ( vartable[varid].productDefinitionTemplate != CDI_UNDEFID )
        vlistDefVarProductDefinitionTemplate(vlistID, varID, vartable[varid].productDefinitionTemplate);

      if ( vartable[varid].lmissval ) vlistDefVarMissval(vlistID, varID, vartable[varid].missval);
      if ( vartable[varid].name )     vlistDefVarName(vlistID, varID, vartable[varid].name);
      if ( vartable[varid].stdname )  vlistDefVarStdname(vlistID, varID, vartable[varid].stdname);
      if ( vartable[varid].longname ) vlistDefVarLongname(vlistID, varID, vartable[varid].longname);
      if ( vartable[varid].units )    vlistDefVarUnits(vlistID, varID, vartable[varid].units);

      if ( vartable[varid].ensdata )  vlistDefVarEnsemble(vlistID, varID, vartable[varid].ensdata->ens_index,
	                                                  vartable[varid].ensdata->ens_count,
							  vartable[varid].ensdata->forecast_init_type);

      vlist_t *vlistptr = vlist_to_pointer(vlistID);
      for ( int i = 0; i < vartable[varid].opt_grib_nentries; i++ )
        {
          resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries+1);
          vlistptr->vars[varID].opt_grib_nentries += 1;
          int idx = vlistptr->vars[varID].opt_grib_nentries-1;

          vlistptr->vars[varID].opt_grib_kvpair[idx] = vartable[varid].opt_grib_kvpair[i];
          vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = NULL;
	  if ( vartable[varid].opt_grib_kvpair[i].keyword )
	    vlistptr->vars[varID].opt_grib_kvpair[idx].keyword =
	      strdupx(vartable[varid].opt_grib_kvpair[i].keyword);
          vlistptr->vars[varID].opt_grib_kvpair[i].update = true;
        }
      /* note: if the key is not defined, we do not throw an error! */

      if ( cdiDefaultTableID != CDI_UNDEFID )
	{
	  int pdis, pcat, pnum;
	  cdiDecodeParam(param, &pnum, &pcat, &pdis);
          char name[CDI_MAX_NAME]; name[0] = 0;
          char longname[CDI_MAX_NAME]; longname[0] = 0;
          char units[CDI_MAX_NAME]; units[0] = 0;
          tableInqEntry(cdiDefaultTableID, pnum, -1, name, longname, units);
	  if ( name[0] )
	    {
	      if ( tableID != CDI_UNDEFID )
		{
		  vlistDefVarName(vlistID, varID, name);
                  if ( longname[0] ) vlistDefVarLongname(vlistID, varID, longname);
                  if ( units[0] ) vlistDefVarUnits(vlistID, varID, units);
                }
	      else
		tableID = cdiDefaultTableID;
	    }
	  if ( cdiDefaultModelID != CDI_UNDEFID ) modelID = cdiDefaultModelID;
	  if ( cdiDefaultInstID  != CDI_UNDEFID )  instID = cdiDefaultInstID;
	}

      if ( instID  != CDI_UNDEFID ) vlistDefVarInstitut(vlistID, varID, instID);
      if ( modelID != CDI_UNDEFID ) vlistDefVarModel(vlistID, varID, modelID);
      if ( tableID != CDI_UNDEFID ) vlistDefVarTable(vlistID, varID, tableID);
    }

  for ( size_t index = 0; index < varTableUsed; index++ )
    {
      int varid = varids[index];
      unsigned nlevels = vartable[varid].recordTable[0].nlevels;

      unsigned nsub = vartable[varid].nsubtypes >= 0 ? (unsigned)vartable[varid].nsubtypes : 0U;
      for ( size_t isub = 0; isub < nsub; isub++ )
        {
          sleveltable_t *restrict streamRecordTable
            = streamptr->vars[index].recordTable + isub;
          leveltable_t *restrict vartableLevelTable
            = vartable[varid].recordTable[isub].levelTable;
          for ( unsigned levelID = 0; levelID < nlevels; levelID++ )
            {
              streamRecordTable->recordID[levelID] = vartableLevelTable[levelID].recID;
              unsigned lindex;
              for ( lindex = 0; lindex < nlevels; lindex++ )
                if ( levelID == (unsigned)vartableLevelTable[lindex].lindex )
                  break;
              if ( lindex == nlevels )
                Error("Internal problem! lindex not found.");
              streamRecordTable->lindex[levelID] = (int)lindex;
            }
        }
    }

  Free(varids);

  varFree();
}


void varDefVCT(size_t vctsize, double *vctptr)
{
  if ( Vct == NULL && vctptr != NULL && vctsize > 0 )
    {
      Vctsize = vctsize;
      Vct = (double *) Malloc(vctsize*sizeof(double));
      memcpy(Vct, vctptr, vctsize*sizeof(double));
    }
}


void varDefZAxisReference(int nhlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE])
{
  numberOfVerticalLevels = nhlev;
  numberOfVerticalGrid = nvgrid;
  memcpy(uuidVGrid, uuid, CDI_UUID_SIZE);
}


bool zaxisCompare(int zaxisID, int zaxistype, int nlevels, bool lbounds, const double *levels, const char *longname, const char *units, int ltype1)
{
  bool differ = true;

  bool ltype_is_equal = (ltype1 == zaxisInqLtype(zaxisID));

  if ( ltype_is_equal && (zaxistype == zaxisInqType(zaxisID) || zaxistype == ZAXIS_GENERIC) )
    {
      bool zlbounds = (zaxisInqLbounds(zaxisID, NULL) > 0);
      if ( nlevels == zaxisInqSize(zaxisID) && zlbounds == lbounds )
	{
	  const double *dlevels = zaxisInqLevelsPtr(zaxisID);
          if ( dlevels && levels )
            {
              int levelID;
              for ( levelID = 0; levelID < nlevels; levelID++ )
                {
                  if ( fabs(dlevels[levelID] - levels[levelID]) > 1.e-9 )
                    break;
                }
              if ( levelID == nlevels ) differ = false;
            }

	  if ( ! differ )
	    {
              if ( longname && longname[0] )
                {
                  char zlongname[CDI_MAX_NAME]; zlongname[0] = 0;
                  zaxisInqLongname(zaxisID, zlongname);
                  if ( zlongname[0] && (strcmp(longname, zlongname) != 0) ) differ = true;
                }
              if ( units && units[0] )
                {
                  char zunits[CDI_MAX_NAME]; zunits[0] = 0;
                  zaxisInqUnits(zaxisID, zunits);
                  if ( zunits[0] && (strcmp(units, zunits) != 0) ) differ = true;
                }
            }
	}
    }

  return differ;
}

struct varDefZAxisSearchState
{
  int resIDValue;
  int zaxistype;
  int nlevels;
  bool lbounds;
  const double *levels;
  const char *longname;
  const char *units;
  int ltype;
};

static enum cdiApplyRet
varDefZAxisSearch(int id, void *res, void *data)
{
  struct varDefZAxisSearchState *state = (struct varDefZAxisSearchState *)data;
  (void)res;
  if ( zaxisCompare(id, state->zaxistype, state->nlevels, state->lbounds,
                    state->levels, state->longname, state->units, state->ltype)
      == false)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}


int varDefZaxis(int vlistID, int zaxistype, int nlevels, const double *levels, const char **cvals, size_t clength, bool lbounds,
                const double *levels1, const double *levels2, int vctsize, const double *vct, char *name,
		const char *longname, const char *units, int prec, int mode, int ltype1)
{
  /*
    mode: 0 search in vlist and zaxis table
          1 search in zaxis table
   */
  int zaxisID = CDI_UNDEFID;
  bool zaxisdefined = false;
  bool zaxisglobdefined = false;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int nzaxis = vlistptr->nzaxis;

  if ( mode == 0 )
    for ( int index = 0; index < nzaxis; index++ )
      {
	zaxisID = vlistptr->zaxisIDs[index];

	if ( !zaxisCompare(zaxisID, zaxistype, nlevels, lbounds, levels, longname, units, ltype1) )
	  {
	    zaxisdefined = true;
	    break;
	  }
      }

  if ( ! zaxisdefined )
    {
      struct varDefZAxisSearchState query;
      query.zaxistype = zaxistype;
      query.nlevels = nlevels;
      query.levels = levels;
      query.lbounds = lbounds;
      query.longname = longname;
      query.units = units;
      query.ltype = ltype1;

      if ((zaxisglobdefined
           = (cdiResHFilterApply(getZaxisOps(), varDefZAxisSearch, &query)
              == CDI_APPLY_STOP)))
        zaxisID = query.resIDValue;

      if ( mode == 1 && zaxisglobdefined)
	for (int index = 0; index < nzaxis; index++ )
	  if ( vlistptr->zaxisIDs[index] == zaxisID )
	    {
	      zaxisglobdefined = false;
	      break;
	    }
    }

  if ( ! zaxisdefined )
    {
      if ( ! zaxisglobdefined )
	{
	  zaxisID = zaxisCreate(zaxistype, nlevels);
	  if ( levels ) zaxisDefLevels(zaxisID, levels);
	  if ( lbounds )
	    {
	      zaxisDefLbounds(zaxisID, levels1);
	      zaxisDefUbounds(zaxisID, levels2);
	    }

          if ( cvals != NULL && nlevels != 0 && clength != 0 ) zaxisDefCvals(zaxisID, cvals, clength);

	  if ( (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && vctsize > 0 )
            zaxisDefVct(zaxisID, vctsize, vct);

	  if ( name && name[0] ) zaxisDefName(zaxisID, name);
	  if ( longname && longname[0] ) zaxisDefLongname(zaxisID, longname);
          else                           zaxisDefLongname(zaxisID, "");
	  if ( units && units[0] ) zaxisDefUnits(zaxisID, units);
          else                     zaxisDefUnits(zaxisID, "");
	  zaxisDefDatatype(zaxisID, prec);
	  zaxisDefLtype(zaxisID, ltype1);
	}

      vlistptr->zaxisIDs[nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  return zaxisID;
}


void varDefMissval(int varID, double missval)
{
  vartable[varID].lmissval = true;
  vartable[varID].missval = missval;
}


void varDefCompType(int varID, int comptype)
{
  if ( vartable[varID].comptype == CDI_COMPRESS_NONE )
    vartable[varID].comptype = comptype;
}


void varDefCompLevel(int varID, int complevel)
{
  vartable[varID].complevel = complevel;
}


int varInqInst(int varID)
{
  return vartable[varID].instID;
}


void varDefInst(int varID, int instID)
{
  vartable[varID].instID = instID;
}


int varInqModel(int varID)
{
  return vartable[varID].modelID;
}


void varDefModel(int varID, int modelID)
{
  vartable[varID].modelID = modelID;
}


int varInqTable(int varID)
{
  return vartable[varID].tableID;
}


void varDefTable(int varID, int tableID)
{
  vartable[varID].tableID = tableID;
}


void varDefEnsembleInfo(int varID, int ens_idx, int ens_count, int forecast_type)
{
  if ( vartable[varID].ensdata == NULL )
      vartable[varID].ensdata = (ensinfo_t *) Malloc( sizeof( ensinfo_t ) );

  vartable[varID].ensdata->ens_index = ens_idx;
  vartable[varID].ensdata->ens_count = ens_count;
  vartable[varID].ensdata->forecast_init_type = forecast_type;
}


void varDefTypeOfGeneratingProcess(int varID, int typeOfGeneratingProcess)
{
  vartable[varID].typeOfGeneratingProcess = typeOfGeneratingProcess;
}


void varDefProductDefinitionTemplate(int varID, int productDefinitionTemplate)
{
  vartable[varID].productDefinitionTemplate = productDefinitionTemplate;
}

#if  defined  (HAVE_LIBGRIB_API)
/* Resizes and initializes opt_grib_kvpair data structure. */
static
void resize_vartable_opt_grib_entries(vartable_t *var, int nentries)
{
  if (var->opt_grib_kvpair_size >= nentries)
    {
      return;   /* nothing to do; array is still large enough */
    }
  else
    {
      if ( CDI_Debug )
        Message("resize data structure, %d -> %d", var->opt_grib_kvpair_size, nentries);

      int i, new_size;
      new_size = (2*var->opt_grib_kvpair_size) > nentries ? (2*var->opt_grib_kvpair_size) : nentries;
      if (CDI_Debug)
        Message("resize vartable opt_grib_entries array to size %d", new_size);
      opt_key_val_pair_t *tmp = (opt_key_val_pair_t *) Malloc((size_t)new_size * sizeof (opt_key_val_pair_t));
      for (i=0; i<var->opt_grib_kvpair_size; i++) {
        tmp[i] = var->opt_grib_kvpair[i];
      }
      for (i=var->opt_grib_kvpair_size; i<new_size; i++) {
        tmp[i].int_val =     0;
        tmp[i].dbl_val =     0;
        tmp[i].update  = false;
        tmp[i].keyword =  NULL;
      } // for
      var->opt_grib_kvpair_size = new_size;
      Free(var->opt_grib_kvpair);
      var->opt_grib_kvpair = tmp;
    }
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
void varDefOptGribInt(int varID, int tile_index, long lval, const char *keyword)
{
  int idx = -1;
  for (int i=0; i<vartable[varID].opt_grib_nentries; i++)
    {
      if ( (strcmp(keyword, vartable[varID].opt_grib_kvpair[i].keyword) == 0 ) &&
           (vartable[varID].opt_grib_kvpair[i].data_type == t_int)             &&
           (vartable[varID].opt_grib_kvpair[i].subtype_index == tile_index) )
        idx = i;
    }

  if (idx == -1)
    {
      resize_vartable_opt_grib_entries(&vartable[varID], vartable[varID].opt_grib_nentries+1);
      vartable[varID].opt_grib_nentries += 1;
      idx = vartable[varID].opt_grib_nentries -1;
    }
  else
    {
      if (vartable[varID].opt_grib_kvpair[idx].keyword)
        Free(vartable[varID].opt_grib_kvpair[idx].keyword);
    }
  vartable[varID].opt_grib_kvpair[idx].data_type     = t_int;
  vartable[varID].opt_grib_kvpair[idx].int_val       = (int) lval;
  vartable[varID].opt_grib_kvpair[idx].keyword       = strdupx(keyword);
  vartable[varID].opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
void varDefOptGribDbl(int varID, int tile_index, double dval, const char *keyword)
{
  int idx = -1;
  for (int i=0; i<vartable[varID].opt_grib_nentries; i++)
    {
      if ( (strcmp(keyword, vartable[varID].opt_grib_kvpair[i].keyword) == 0 ) &&
           (vartable[varID].opt_grib_kvpair[i].data_type == t_double)          &&
           (vartable[varID].opt_grib_kvpair[i].subtype_index == tile_index) )
        idx = i;
    }

  if (idx == -1)
    {
      resize_vartable_opt_grib_entries(&vartable[varID], vartable[varID].opt_grib_nentries+1);
      vartable[varID].opt_grib_nentries += 1;
      idx = vartable[varID].opt_grib_nentries -1;
    }
  else
    {
      if (vartable[varID].opt_grib_kvpair[idx].keyword)
        Free(vartable[varID].opt_grib_kvpair[idx].keyword);
    }
  vartable[varID].opt_grib_kvpair[idx].data_type     = t_double;
  vartable[varID].opt_grib_kvpair[idx].dbl_val       = dval;
  vartable[varID].opt_grib_kvpair[idx].keyword       = strdupx(keyword);
  vartable[varID].opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
int varOptGribNentries(int varID)
{
  int nentries = vartable[varID].opt_grib_nentries;
  return nentries;
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
