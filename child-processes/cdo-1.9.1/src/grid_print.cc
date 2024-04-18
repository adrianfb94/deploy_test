#include <cdi.h>
#include "cdi_uuid.h"
#include "cdo_int.h"

static
void printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], int nbyte0, size_t n, const double vals[])
{
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static
void printIntsPrefixAutoBrk(FILE *fp, const char prefix[], int nbyte0, size_t n, const int vals[])
{
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%d ", vals[i]);
    }
  fputs("\n", fp);
}

static void
printBounds(FILE *fp, int dig, const char prefix[], size_t nbyte0,
            size_t n, size_t nvertex, const double bounds[])
{
  fputs(prefix, fp);
  for ( size_t i = 0; i < n; i++ )
    {
      if ( i > 0 ) fprintf(fp, "\n%*s", (int)nbyte0, "");
      for ( size_t iv = 0; iv < nvertex; iv++ )
        fprintf(fp, "%.*g ", dig, bounds[i*nvertex+iv]);
    }
  fputs("\n", fp);
}

static
void printMask(FILE *fp, const char prefix[], size_t nbyte0, size_t n, const int mask[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%d ", mask[i]);
    }
  fputs("\n", fp);
}

static inline
void *resizeBuffer(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static
void grid_print_attributes(FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  int atttype, attlen;
  char attname[CDI_MAX_NAME+1];
  void *attBuf = NULL;
  size_t attBufSize = 0;
  char fltstr[128];

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for ( int iatt = 0; iatt < natts; ++iatt )
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      if ( attlen == 0 ) continue;

      if ( strcmp(attname, "grid_mapping_name") == 0 ) continue;

      if ( atttype == CDI_DATATYPE_TXT )
        {
          size_t attSize = (size_t)(attlen+1)*sizeof(char);
          char *atttxt = (char *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          atttxt[attlen] = 0;
          fprintf(fp, "%s = \"%s\"\n", attname, atttxt);
        }
      else if ( atttype == CDI_DATATYPE_INT8  || atttype == CDI_DATATYPE_UINT8  ||
                atttype == CDI_DATATYPE_INT16 || atttype == CDI_DATATYPE_UINT16 ||
                atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32 )
        {
          size_t attSize = (size_t)attlen*sizeof(int);
          int *attint = (int *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          fprintf(fp, "%s =", attname);
          for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %d", attint[i]);
          fprintf(fp, "\n");
        }
      else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
        {
          size_t attSize = (size_t)attlen * sizeof(double);
          double *attflt = (double *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          fprintf(fp, "%s =", attname);
          if ( atttype == CDI_DATATYPE_FLT32 )
            for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %sf", double_to_attstr(CDO_flt_digits, fltstr, sizeof(fltstr), attflt[i]));
          else
            for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %s", double_to_attstr(CDO_dbl_digits, fltstr, sizeof(fltstr), attflt[i]));
          fprintf(fp, "\n");
        }
    }

  Free(attBuf);
}

static
void grid_print_kernel(int gridID, int opt, FILE *fp)
{
  int xdim, ydim;
  char attstr[CDI_MAX_NAME];
  char attstr2[CDI_MAX_NAME];
  size_t nxvals = (size_t) gridInqXvals(gridID, NULL);
  size_t nyvals = (size_t) gridInqYvals(gridID, NULL);
  size_t nxbounds = (size_t) gridInqXbounds(gridID, NULL);
  size_t nybounds = (size_t) gridInqYbounds(gridID, NULL);

  int type     = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  int xsize    = gridInqXsize(gridID);
  int ysize    = gridInqYsize(gridID);
  int nvertex  = gridInqNvertex(gridID);
  int prec     = gridInqDatatype(gridID);
  int xstrlen  = gridInqXIsc(gridID);
  int ystrlen  = gridInqYIsc(gridID);

  char **xcvals = NULL;
  if ( xstrlen )
    {
      xcvals = (char **) Malloc(xsize * sizeof(char *));
      for ( int i = 0; i < xsize; i++ ) xcvals[i] = (char*) Malloc((xstrlen+1) * sizeof(char));
      gridInqXCvals(gridID, xcvals);
      for ( int i = 0; i < xsize; i++ ) xcvals[i][xstrlen] = 0;
      for ( int i = 0; i < xsize; i++ ) 
        for ( int k = xstrlen-1; k; k-- ) if ( xcvals[i][k] == ' ' ) xcvals[i][k] = 0; else break;
    }

  char **ycvals = NULL;
  if ( ystrlen )
    {
      ycvals = (char **) Malloc(ysize * sizeof(char *));
      for ( int i = 0; i < ysize; i++ ) ycvals[i] = (char*) Malloc((ystrlen+1) * sizeof(char));
      gridInqYCvals(gridID, ycvals);
      for ( int i = 0; i < ysize; i++ ) ycvals[i][ystrlen] = 0;
      for ( int i = 0; i < ysize; i++ ) 
        for ( int k = ystrlen-1; k; k-- ) if ( ycvals[i][k] == ' ' ) ycvals[i][k] = 0; else break;
    }

  int dig = (prec == CDI_DATATYPE_FLT64) ? CDO_dbl_digits : CDO_flt_digits;

  fprintf(fp, "gridtype  = %s\n" "gridsize  = %d\n", gridNamePtr(type), gridsize);

  if ( type != GRID_GME )
    {
      if ( type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER )
        {
          if ( xsize > 0 ) fprintf(fp, "xsize     = %d\n", xsize);
          if ( ysize > 0 ) fprintf(fp, "ysize     = %d\n", ysize);
        }

      if ( nxvals > 0 || xcvals )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "xdimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xlongname = \"%s\"\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xunits    = \"%s\"\n", attstr);
        }

      if ( nyvals > 0 || ycvals )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "ydimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "ylongname = \"%s\"\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yunits    = \"%s\"\n", attstr);
        }

      if ( type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_VDIMNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] ) fprintf(fp, "vdimname  = %s\n", attstr);
        }
      if ( type == GRID_UNSTRUCTURED && nvertex > 0 ) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_CHARXY:
      {
        if ( type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED ) fprintf(fp, "np        = %d\n", gridInqNP(gridID));

	if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED )
	  {
	    xdim = gridsize;
	    ydim = gridsize;
	  }
        else if ( type == GRID_GAUSSIAN_REDUCED )
          {
	    xdim = 2;
	    ydim = ysize;
          }
	else
	  {
	    xdim = xsize;
	    ydim = ysize;
	  }

	if ( type == GRID_UNSTRUCTURED )
          {
            int number = gridInqNumber(gridID);
            int position = gridInqPosition(gridID);
            // const unsigned char *d;
            if ( number > 0 )
              {
                fprintf(fp, "number    = %d\n", number);
                if ( position >= 0 ) fprintf(fp, "position  = %d\n", position);
              }

            if ( gridInqReference(gridID, NULL) )
              {
                char reference_link[8192];
                gridInqReference(gridID, reference_link);
                fprintf(fp, "uri       = %s\n", reference_link);
              }
          }

	if ( nxvals > 0 )
	  {
	    double xfirst = 0.0, xinc = 0.0;

	    if ( type == GRID_LONLAT     || type == GRID_GAUSSIAN ||
		 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		xfirst = gridInqXval(gridID, 0);
		xinc   = gridInqXinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(xinc, 0) && opt )
	      {
                fprintf(fp, "xfirst    = %.*g\n"
                        "xinc      = %.*g\n", dig, xfirst, dig, xinc);
	      }
	    else
	      {
                double *xvals = (double*) Malloc(nxvals*sizeof(double));
                gridInqXvals(gridID, xvals);
                static const char prefix[] = "xvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nxvals, xvals);
                Free(xvals);
	      }
	  }

	if ( nxbounds )
	  {
            double *xbounds = (double*) Malloc(nxbounds*sizeof(double));
            gridInqXbounds(gridID, xbounds);
            static const char prefix[] = "xbounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, xdim, nvertex, xbounds);
            Free(xbounds);
	  }

        if ( xcvals )
          {
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, attstr);
            if ( attstr[0] )
              fprintf(fp, "x%ss  = \"%.*s\"", attstr, xstrlen, xcvals[0]);
            else
              fprintf(fp, "xstrings  = \"%.*s\"", xstrlen, xcvals[0]);
            for ( int i = 1; i < xsize; i++ )
              fprintf(fp, ", \"%.*s\"", xstrlen, xcvals[i]);
            fprintf(fp, "\n");

            for ( int i = 0; i < xsize; i++ ) Free(xcvals[i]);
            Free(xcvals);
          }

	if ( nyvals > 0 )
	  {
	    double yfirst = 0.0, yinc = 0.0;

	    if ( type == GRID_LONLAT || type == GRID_GENERIC ||
                 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		yfirst = gridInqYval(gridID, 0);
		yinc   = gridInqYinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(yinc, 0) && opt )
	      {
	  	fprintf(fp, "yfirst    = %.*g\n"
                        "yinc      = %.*g\n", dig, yfirst, dig, yinc);
	      }
	    else
	      {
                double *yvals = (double*) Malloc(nyvals*sizeof(double));
                gridInqYvals(gridID, yvals);
                static const char prefix[] = "yvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nyvals, yvals);
                Free(yvals);
	      }
	  }

	if ( nybounds )
	  {
            double *ybounds = (double*) Malloc(nybounds*sizeof(double));
            gridInqYbounds(gridID, ybounds);
            static const char prefix[] = "ybounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, ydim, nvertex, ybounds);
            Free(ybounds);
	  }

        if ( ycvals )
          {
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, attstr);
            if ( attstr[0] )
              fprintf(fp, "x%ss  = \"%.*s\"", attstr, ystrlen, ycvals[0]);
            else
              fprintf(fp, "ystrings  = \"%.*s\"", ystrlen, ycvals[0]);
            for ( int i = 1; i < ysize; i++ )
              fprintf(fp, ", \"%.*s\"", ystrlen, ycvals[i]);
            fprintf(fp, "\n");

            for ( int i = 0; i < ysize; i++ ) Free(ycvals[i]);
            Free(ycvals);
          }

	if ( gridHasArea(gridID) )
	  {
            double *area = (double*) Malloc(gridsize*sizeof(double));
            gridInqArea(gridID, area);
            static const char prefix[] = "area      = ";
            printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, gridsize, area);
            Free(area);
	  }

        if ( type == GRID_GAUSSIAN_REDUCED )
          {
            static const char prefix[] = "rowlon    = ";
            int *rowlon = (int *)Malloc((size_t)ysize*sizeof(int));
            gridInqRowlon(gridID, rowlon);
            printIntsPrefixAutoBrk(fp, prefix, (int)sizeof(prefix)-1,
                                   (size_t)(ysize > 0 ? ysize : 0), rowlon);
            Free(rowlon);
          }


        int uvRelativeToGrid = gridInqUvRelativeToGrid(gridID);
        if ( uvRelativeToGrid > 0 )
          fprintf(fp, "uvRelativeToGrid = %d\n", uvRelativeToGrid);

#ifdef HIRLAM_EXTENSIONS
        {
          int scanningMode = gridInqScanningMode(gridID);
          fprintf(fp, "scanningMode = %d\n", scanningMode);
        }
#endif // HIRLAM_EXTENSIONS

        if ( type == GRID_PROJECTION )
          {
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_MAPPING, CDI_MAX_NAME, attstr);
            if ( attstr[0] ) fprintf(fp, "grid_mapping = %s\n", attstr);
            attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, attstr);
            if ( attstr[0] ) fprintf(fp, "grid_mapping_name = %s\n", attstr);
            grid_print_attributes(fp, gridID);
          }

	break;
      }
    case GRID_SPECTRAL:
      {
        fprintf(fp, "truncation = %d\n"
                "complexpacking = %d\n", gridInqTrunc(gridID), gridInqComplexPacking(gridID) );
        break;
      }
    case GRID_FOURIER:
      {
	fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
	break;
      }
    case GRID_GME:
      {
        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
        fprintf(fp, "ni        = %d\n", ni);
        break;
      }
   default:
      {
	fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
        break;
      }
    }

  unsigned char uuid[CDI_UUID_SIZE];
  gridInqUUID(gridID, uuid);
  if ( !cdiUUIDIsNull(uuid) )
    {
      char uuidStr[37];
      cdiUUID2Str(uuid, uuidStr);
      if ( uuidStr[0] != 0 && strlen(uuidStr) == 36 )
        fprintf(fp, "uuid      = %s\n", uuidStr);
    }

  if ( gridInqMask(gridID, NULL) )
    {
      int *mask = (gridsize>0) ? (int*) Malloc((size_t)gridsize*sizeof(int)) : NULL;
      gridInqMask(gridID, mask);
      static const char prefix[] = "mask      = ";
      printMask(fp, prefix, sizeof(prefix)-1, (size_t)(gridsize > 0 ? gridsize : 0), mask);
      if ( mask ) Free(mask);
    }    

  int projID = gridInqProj(gridID);
  if ( projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION )
    grid_print_kernel(projID, opt, fp);
}


void cdo_print_grid(int gridID, int opt)
{
  grid_print_kernel(gridID, opt, stdout);
}
