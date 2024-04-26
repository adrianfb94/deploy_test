#include <cdi.h>
#include "cdi_uuid.h"
#include "cdo_int.h"

void cdo_print_attributes(FILE *fp, int cdiID, int varID, int nblanks);

static
void printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], int nbyte0, size_t n, const double vals[], size_t extbreak)
{
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80  || (i && i == extbreak) )
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static
void zaxis_print_kernel(int zaxisID, FILE *fp)
{
  char attstr[CDI_MAX_NAME];
  int type    = zaxisInqType(zaxisID);
  int nlevels = zaxisInqSize(zaxisID);
  int prec    = zaxisInqDatatype(zaxisID);
  size_t nvals = (size_t) zaxisInqLevels(zaxisID, NULL);

  int dig = (prec == CDI_DATATYPE_FLT64) ? CDO_dbl_digits : CDO_flt_digits;

  fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  fprintf(fp, "size      = %d\n", nlevels);

  if ( nlevels == 1 && zaxisInqScalar(zaxisID) ) fprintf(fp, "scalar    = true\n");

  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_NAME, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "name      = %s\n", attstr);
  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_LONGNAME, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "longname  = \"%s\"\n", attstr);
  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_UNITS, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "units     = \"%s\"\n", attstr);

  char **cvals = NULL;
  double *vals = nvals ? (double*) Malloc(nvals*sizeof(double)) : NULL;

  if ( nvals )
    {                
      zaxisInqLevels(zaxisID, vals);
      static const char prefix[] = "levels    = ";
      printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
    }
  else if ( type == ZAXIS_CHAR )
    {
      int clen = zaxisInqCLen(zaxisID);
      zaxisInqCVals(zaxisID, &cvals);
      fprintf(fp, "levels    = \n");
      for ( int i = 0; i < nlevels; i++ )
        {
          fprintf(fp, "     [%2d] = %.*s\n", i, clen, cvals[i]);
          Free(cvals[i]);
        }
    }
  


  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      {
        zaxisInqLbounds(zaxisID, vals);
        static const char prefix[] = "lbounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
      }

      {
        zaxisInqUbounds(zaxisID, vals);
        static const char prefix[] = "ubounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
      }
    }

  if ( vals ) Free(vals);
  if ( cvals ) Free(cvals);

  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int vctsize = zaxisInqVctSize(zaxisID);
      if ( vctsize )
        {
          fprintf(fp, "vctsize   = %d\n", vctsize);
          double *vct = (double*) Malloc(vctsize*sizeof(double));
          zaxisInqVct(zaxisID, vct);
          static const char prefix[] = "vct       = ";
          printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, vctsize, vct, vctsize/2);
          Free(vct);
        }
    }

  if ( type == ZAXIS_REFERENCE )
    {
      unsigned char uuid[CDI_UUID_SIZE];
      zaxisInqUUID(zaxisID, uuid);
      if ( !cdiUUIDIsNull(uuid) )
        {
          char uuidStr[37];
          cdiUUID2Str(uuid, uuidStr);
          if ( uuidStr[0] != 0 && strlen(uuidStr) == 36 )
            fprintf(fp, "uuid      = %s\n", uuidStr);
        }
    }

  cdo_print_attributes(fp, zaxisID, CDI_GLOBAL, 0);
}


void cdo_print_zaxis(int zaxisID)
{
  zaxis_print_kernel(zaxisID, stdout);
}
