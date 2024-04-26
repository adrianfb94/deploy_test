// This file is used in CDI and CDO !!!

#if defined (HAVE_CONFIG_H)
#include "../src/config.h"
#endif

#include <stdio.h>

#ifdef CDO
#define  streamInqFiletype        pstreamInqFiletype
#define  streamInqByteorder       pstreamInqByteorder
#define  streamInqTimestep        pstreamInqTimestep
#endif

#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

void my_reset_text_color(FILE *fp)
{
  (void)fp;
#ifdef CDO
  reset_text_color(fp);
#endif
}

void datetime2str(int date, int time, char *datetimestr, int maxlen)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);
  int hour, minute, second;
  cdiDecodeTime(time, &hour, &minute, &second);

  int len = sprintf(datetimestr, DATE_FORMAT "T" TIME_FORMAT, year, month, day, hour, minute, second);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void date2str(int date, char *datestr, int maxlen)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);

  int len = sprintf(datestr, DATE_FORMAT, year, month, day);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void time2str(int time, char *timestr, int maxlen)
{
  int hour, minute, second;
  cdiDecodeTime(time, &hour, &minute, &second);

  int len = sprintf(timestr, TIME_FORMAT, hour, minute, second);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void printFiletype(int streamID, int vlistID)
{
  int filetype = streamInqFiletype(streamID);

  // clang-format off
  switch ( filetype )
    {
    case CDI_FILETYPE_GRB:  printf("GRIB");  break;
    case CDI_FILETYPE_GRB2: printf("GRIB2");  break;
    case CDI_FILETYPE_NC:   printf("NetCDF");  break;
    case CDI_FILETYPE_NC2:  printf("NetCDF2");  break;
    case CDI_FILETYPE_NC4:  printf("NetCDF4");  break;
    case CDI_FILETYPE_NC4C: printf("NetCDF4 classic");  break;
    case CDI_FILETYPE_NC5:  printf("NetCDF5");  break;
    case CDI_FILETYPE_SRV:  printf("SERVICE");  break;
    case CDI_FILETYPE_EXT:  printf("EXTRA");  break;
    case CDI_FILETYPE_IEG:  printf("IEG");  break;
    default: printf("  File format: unsupported filetype %d" , filetype);  break;
    }

  if ( filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG )
    {
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:    printf("  BIGENDIAN");  break;
	case CDI_LITTLEENDIAN: printf("  LITTLEENDIAN");  break;
	default: printf("  byteorder: %d undefined", streamInqByteorder(streamID));  break;
	}
    }
  // clang-format on

  if ( filetype == CDI_FILETYPE_GRB || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C )
    {
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if      ( comptype == CDI_COMPRESS_SZIP ) printf(" SZIP");
	      else if ( comptype == CDI_COMPRESS_ZIP  ) printf(" ZIP");

	      break;
	    }
	}
    }

  if ( filetype == CDI_FILETYPE_GRB2 )
    {
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == CDI_COMPRESS_JPEG ) printf(" JPEG");
	      break;
	    }
	}
    }

  printf("\n");
}

static
void print_xvals(int gridID, int dig)
{
  int xsize = gridInqXsize(gridID);
  if ( xsize > 0 && gridInqXvals(gridID, NULL) )
    {
      char xname[CDI_MAX_NAME], xunits[CDI_MAX_NAME];
      gridInqXname(gridID, xname);
      gridInqXunits(gridID, xunits);

      double xfirst = gridInqXval(gridID, 0);
      double xlast  = gridInqXval(gridID, xsize-1);
      double xinc   = gridInqXinc(gridID);
      fprintf(stdout, "%33s : %.*g", xname, dig, xfirst);
      if ( xsize > 1 )
        {
          fprintf(stdout, " to %.*g", dig, xlast);
          if ( IS_NOT_EQUAL(xinc, 0) )
            fprintf(stdout, " by %.*g", dig, xinc);
        }
      fprintf(stdout, " %s", xunits);
      if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
      fprintf(stdout, "\n");
    }
}

static
void print_yvals(int gridID, int dig)
{
  int ysize = gridInqYsize(gridID);
  if ( ysize > 0 && gridInqYvals(gridID, NULL) )
    {
      char yname[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      gridInqYname(gridID, yname);
      gridInqYunits(gridID, yunits);

      double yfirst = gridInqYval(gridID, 0);
      double ylast  = gridInqYval(gridID, ysize-1);
      double yinc   = gridInqYinc(gridID);
      fprintf(stdout, "%33s : %.*g", yname, dig, yfirst);
      if ( ysize > 1 )
        {
          int gridtype = gridInqType(gridID);
          fprintf(stdout, " to %.*g", dig, ylast);
          if ( IS_NOT_EQUAL(yinc, 0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED )
            fprintf(stdout, " by %.*g", dig, yinc);
        }
      fprintf(stdout, " %s", yunits);
      fprintf(stdout, "\n");
    }
}

static
void print_xyvals2D(int gridID, int dig)
{
  if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
    {
      char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME], xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      gridInqXname(gridID, xname);
      gridInqYname(gridID, yname);
      gridInqXunits(gridID, xunits);
      gridInqYunits(gridID, yunits);

      size_t gridsize = gridInqSize(gridID);
      double *xvals2D = (double*) malloc(gridsize*sizeof(double));
      double *yvals2D = (double*) malloc(gridsize*sizeof(double));

      gridInqXvals(gridID, xvals2D);
      gridInqYvals(gridID, yvals2D);

      double xfirst = xvals2D[0];
      double xlast  = xvals2D[0];
      double yfirst = yvals2D[0];
      double ylast  = yvals2D[0];
      for ( size_t i = 1; i < gridsize; i++ )
        {
          if ( xvals2D[i] < xfirst ) xfirst = xvals2D[i];
          if ( xvals2D[i] > xlast  ) xlast  = xvals2D[i];
          if ( yvals2D[i] < yfirst ) yfirst = yvals2D[i];
          if ( yvals2D[i] > ylast  ) ylast  = yvals2D[i];
        }

      double xinc = 0;
      double yinc = 0;
      int gridtype = gridInqType(gridID);
      if ( gridtype == GRID_CURVILINEAR )
        {
          int xsize = gridInqXsize(gridID);
          if ( xsize > 1 )
            {
              double *xvals = (double*) malloc((size_t)xsize*sizeof(double));
              for ( int i = 0; i < xsize; ++i ) xvals[i] = xvals2D[i];
              xinc = fabs(xvals[xsize-1] - xvals[0])/(xsize-1);
              for ( int i = 2; i < xsize; i++ )
                if ( fabs(fabs(xvals[i-1] - xvals[i]) - xinc) > 0.01*xinc ) { xinc = 0; break; }
              free(xvals);
            }
          int ysize = gridInqYsize(gridID);
          if ( ysize > 1 )
            {
              double *yvals = (double*) malloc((size_t)ysize*sizeof(double));
              for ( int i = 0; i < ysize; ++i ) yvals[i] = yvals2D[i*xsize];
              yinc = fabs(yvals[ysize-1] - yvals[0])/(ysize-1);
              for ( int i = 2; i < ysize; i++ )
                if ( fabs(fabs(yvals[i-1] - yvals[i]) - yinc) > 0.01*yinc ) { yinc = 0; break; }
              free(yvals);
            }
        }

      fprintf(stdout, "%33s : %.*g", xname, dig, xfirst);
      if ( gridsize > 1 ) fprintf(stdout, " to %.*g", dig, xlast);
      if ( IS_NOT_EQUAL(xinc, 0) ) fprintf(stdout, " by %.*g", dig, xinc);
      fprintf(stdout, " %s", xunits);
      if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
      fprintf(stdout, "\n");
      fprintf(stdout, "%33s : %.*g", yname, dig, yfirst);
      if ( gridsize > 1 ) fprintf(stdout, " to %.*g", dig, ylast);
      if ( IS_NOT_EQUAL(yinc, 0) ) fprintf(stdout, " by %.*g", dig, yinc);
      fprintf(stdout, " %s", yunits);
      fprintf(stdout, "\n");

      free(xvals2D);
      free(yvals2D);
    }
}

static
void printGridInfoKernel(int gridID, int index, bool lproj)
{
  int gridtype = gridInqType(gridID);

  if ( lproj && gridtype != GRID_PROJECTION )
    fprintf(stderr, "Internal problem (%s): sub grid not equal GRID_PROJECTION!\n", __func__);

  int trunc    = gridInqTrunc(gridID);
  size_t gridsize = gridInqSize(gridID);
  size_t xsize    = gridInqXsize(gridID);
  size_t ysize    = gridInqYsize(gridID);
  size_t xysize   = xsize*ysize;

  // int prec     = gridInqDatatype(gridID);
  // int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;
  int dig = 7;
#ifdef CDO
  extern int CDO_flt_digits;
  dig = CDO_flt_digits;
#endif

  if ( !lproj )
    {
      fprintf(stdout, "  %4d : ", index+1);
#ifdef CDO
      set_text_color(stdout, RESET, BLUE);
#endif
      fprintf(stdout, "%-24s", gridNamePtr(gridtype));
      my_reset_text_color(stdout);
      fprintf(stdout, " : ");
    }

  if ( gridtype == GRID_LONLAT     ||
       gridtype == GRID_PROJECTION ||
       gridtype == GRID_GENERIC    ||
       gridtype == GRID_CHARXY    ||
       gridtype == GRID_GAUSSIAN   ||
       gridtype == GRID_GAUSSIAN_REDUCED )
    {
      if ( !lproj )
        {
#ifdef CDO
          set_text_color(stdout, RESET, GREEN);
#endif
          fprintf(stdout, "points=%zu", gridsize);
          if ( gridtype == GRID_GAUSSIAN_REDUCED )
            fprintf(stdout, "  nlat=%zu", ysize);
          else if ( xysize )
            fprintf(stdout, " (%zux%zu)", xsize, ysize);

          if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
            fprintf(stdout, "  np=%d", gridInqNP(gridID));
          my_reset_text_color(stdout);

          fprintf(stdout, "\n");
        }

      char name[CDI_MAX_NAME]; name[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, name);
      if ( gridtype == GRID_PROJECTION || name[0] )
        {
          if ( name[0] == 0 ) strcpy(name, "undefined");
#ifdef CDO
          set_text_color(stdout, RESET, BLUE);
#endif
          fprintf(stdout, "         %24s", "mapping");
          my_reset_text_color(stdout);
          fprintf(stdout, " : ");
#ifdef CDO
          set_text_color(stdout, RESET, GREEN);
#endif
          fprintf(stdout, "%s\n", name);
          my_reset_text_color(stdout);
        }

      print_xvals(gridID, dig);
      print_yvals(gridID, dig);

      if ( gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
        {
          fprintf(stdout, "%33s :", "available");
          if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
          if ( gridHasArea(gridID) )          fprintf(stdout, " area");
          if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
          fprintf(stdout, "\n");
        }
    }
  else if ( gridtype == GRID_SPECTRAL )
    {
#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      fprintf(stdout, "points=%zu  nsp=%zu  truncation=%d", gridsize, gridsize/2, trunc);
      if ( gridInqComplexPacking(gridID) ) fprintf(stdout, "  complexPacking");
      my_reset_text_color(stdout);
      fprintf(stdout, "\n");
    }
  else if ( gridtype == GRID_FOURIER )
    {
#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      fprintf(stdout, "points=%zu  nfc=%zu  truncation=%d\n", gridsize, gridsize/2, trunc);
      my_reset_text_color(stdout);
    }
  else if ( gridtype == GRID_GME )
    {
      int nd, ni, ni2, ni3;
      gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      fprintf(stdout, "points=%zu  nd=%d  ni=%d\n", gridsize, nd, ni);
      my_reset_text_color(stdout);
    }
  else if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
    {
#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      if ( gridtype == GRID_CURVILINEAR )
        fprintf(stdout, "points=%zu (%zux%zu)", gridsize, xsize, ysize);
      else
        fprintf(stdout, "points=%zu", gridsize);

      if ( gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0 )
        fprintf(stdout, "  nvertex=%d", gridInqNvertex(gridID));
      my_reset_text_color(stdout);

      fprintf(stdout, "\n");

      if ( gridtype == GRID_UNSTRUCTURED )
        {
          int number   = gridInqNumber(gridID);
          int position = gridInqPosition(gridID);
          if ( number > 0 )
            fprintf(stdout, "%33s : number=%d  position=%d\n", "grid", number, position);

          if ( gridInqReference(gridID, NULL) )
            {
              char reference_link[8192];
              gridInqReference(gridID, reference_link);
              fprintf(stdout, "%33s : %s\n", "uri", reference_link);
            }
        }

      print_xyvals2D(gridID, dig);
    }
  else /* if ( gridtype == GRID_GENERIC ) */
    {
#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      if ( ysize == 0 )
        fprintf(stdout, "points=%zu\n", gridsize);
      else
        fprintf(stdout, "points=%zu (%zux%zu)\n", gridsize, xsize, ysize);
      my_reset_text_color(stdout);
    }

  if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
    {
      if ( gridHasArea(gridID) ||
           gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
        {
          fprintf(stdout, "%33s :", "available");
          if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
          if ( gridHasArea(gridID) )          fprintf(stdout, " area");
          if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
          fprintf(stdout, "\n");
        }
    }

  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  gridInqUUID(gridID, uuidOfHGrid);
  if ( !cdiUUIDIsNull(uuidOfHGrid) )
    {
      char uuidOfHGridStr[37];
      cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
      if ( uuidOfHGridStr[0] != 0  && strlen(uuidOfHGridStr) == 36 )
        {
          fprintf(stdout, "%33s : %s\n", "uuid", uuidOfHGridStr);
        }
    }
}

static
void printGridInfo(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for ( int index = 0; index < ngrids; index++ )
    {
      int gridID = vlistGrid(vlistID, index);
      printGridInfoKernel(gridID, index, false);
      int projID = gridInqProj(gridID);
      if ( projID != CDI_UNDEFID )
        printGridInfoKernel(projID, index, true);
    }
}

static
void printZaxisInfo(int vlistID)
{
  char zaxisname[CDI_MAX_NAME], zname[CDI_MAX_NAME], zunits[CDI_MAX_NAME];

  int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; index++ )
    {
      double zinc = 0;
      int zaxisID   = vlistZaxis(vlistID, index);
      int zaxistype = zaxisInqType(zaxisID);
      int ltype     = zaxisInqLtype(zaxisID);
      int levelsize = zaxisInqSize(zaxisID);
      // int prec      = zaxisInqDatatype(zaxisID);
      // int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;
      int dig = 7;
#ifdef CDO
      extern int CDO_flt_digits;
      dig = CDO_flt_digits;
#endif

      zaxisName(zaxistype, zaxisname);
      zaxisInqName(zaxisID, zname);
      zaxisInqUnits(zaxisID, zunits);
      zunits[12] = 0;

      fprintf(stdout, "  %4d : ", vlistZaxisIndex(vlistID, zaxisID)+1);
#ifdef CDO
      set_text_color(stdout, RESET, BLUE);
#endif
      if ( zaxistype == ZAXIS_GENERIC && ltype != 0 )
        fprintf(stdout, "%-12s (ltype=%3d)", zaxisname, ltype);
      else
        fprintf(stdout, "%-24s", zaxisname);
      my_reset_text_color(stdout);

      fprintf(stdout, " :");

#ifdef CDO
      set_text_color(stdout, RESET, GREEN);
#endif
      fprintf(stdout, " levels=%d", levelsize);
      bool zscalar = (levelsize == 1) ? zaxisInqScalar(zaxisID) : false;
      if ( zscalar ) fprintf(stdout, "  scalar");
      my_reset_text_color(stdout);
      fprintf(stdout, "\n");

      if ( zaxisInqLevels(zaxisID, NULL) )
        {
          double *levels = (double*) malloc((size_t)levelsize*sizeof(double));
          zaxisInqLevels(zaxisID, levels);

          if ( !(zaxistype == ZAXIS_SURFACE && levelsize == 1 && !(fabs(levels[0]) > 0)) )
            {
              double zfirst = levels[0];
              double zlast  = levels[levelsize-1];
              if ( levelsize > 2 )
                {
                  zinc = (levels[levelsize-1] - levels[0]) / (levelsize-1);
                  for ( int levelID = 2; levelID < levelsize; ++levelID )
                    if ( fabs(fabs(levels[levelID] - levels[levelID-1]) - zinc) > 0.001*zinc )
                      {
                        zinc = 0;
                        break;
                      }
                }

              fprintf(stdout, "%33s : %.*g", zname, dig, zfirst);
              if ( levelsize > 1 )
                {
                  fprintf(stdout, " to %.*g", dig, zlast);
                  if ( IS_NOT_EQUAL(zinc, 0) )
                    fprintf(stdout, " by %.*g", dig, zinc);
                }
              fprintf(stdout, " %s", zunits);
              fprintf(stdout, "\n");
            }

          free(levels);
        }

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
        {
          double level1, level2;
          fprintf(stdout, "%33s : ", "bounds");

          level1 = zaxisInqLbound(zaxisID, 0);
          level2 = zaxisInqUbound(zaxisID, 0);
          fprintf(stdout, "%.*g-%.*g", dig, level1, dig, level2);
          if ( levelsize > 1 )
            {
              level1 = zaxisInqLbound(zaxisID, levelsize-1);
              level2 = zaxisInqUbound(zaxisID, levelsize-1);
              fprintf(stdout, " to %.*g-%.*g", dig, level1, dig, level2);
              if ( IS_NOT_EQUAL(zinc, 0) )
                fprintf(stdout, " by %.*g", dig, zinc);
            }
          fprintf(stdout, " %s", zunits);
          fprintf(stdout, "\n");
        }

      if ( zaxistype == ZAXIS_HYBRID )
        {
          char psname[CDI_MAX_NAME]; psname[0] = 0;
          cdiZaxisInqKeyStr(zaxisID, CDI_KEY_PSNAME, CDI_MAX_NAME, psname);
          int vctsize = zaxisInqVctSize(zaxisID);
          if ( vctsize || psname[0] )
            {
	      fprintf(stdout, "%33s :", "available");
              if ( vctsize   ) fprintf(stdout, " vct");
              if ( psname[0] ) fprintf(stdout, "  ps: %s", psname);
              fprintf(stdout, "\n");
            }
        }

      if ( zaxistype == ZAXIS_REFERENCE )
        {
          int number   = zaxisInqNumber(zaxisID);

          if ( number > 0 )
            {
              fprintf(stdout, "%33s : ", "zaxis");
              fprintf(stdout, "number=%d\n", number);
            }

          unsigned char uuidOfVGrid[CDI_UUID_SIZE];
          zaxisInqUUID(zaxisID, uuidOfVGrid);
          if ( !cdiUUIDIsNull(uuidOfVGrid) )
            {
              char uuidOfVGridStr[37];
              cdiUUID2Str(uuidOfVGrid, uuidOfVGridStr);
              if ( uuidOfVGridStr[0] != 0  && strlen(uuidOfVGridStr) == 36 )
                {
                  fprintf(stdout, "%33s : ", "uuid");
                  fprintf(stdout, "%s\n", uuidOfVGridStr);
                }
            }
        }
    }
}

static
void printSubtypeInfo(int vlistID)
{
  int nsubtypes = vlistNsubtypes(vlistID);
  for ( int index = 0; index < nsubtypes; index++)
    {
      int subtypeID = vlistSubtype(vlistID, index);
      int subtypesize = subtypeInqSize(subtypeID);
      // subtypePrint(subtypeID);
      fprintf(stdout, "  %4d : %-24s :", vlistSubtypeIndex(vlistID, subtypeID)+1, "tiles");
      fprintf(stdout, " ntiles=%d", subtypesize);
      fprintf(stdout, "\n");
    }
}

static
int printDateTime(int ntimeout, int vdate, int vtime)
{
  char vdatestr[32], vtimestr[32];

  if ( ntimeout == 4 )
    {
      ntimeout = 0;
      fprintf(stdout, "\n");
    }

  date2str(vdate, vdatestr, sizeof(vdatestr));
  time2str(vtime, vtimestr, sizeof(vtimestr));

  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  return ++ntimeout;
}

#define NUM_TIMESTEP 60
#define MAX_DOTS     80

static
int printDot(int ndotout, int *nfact, int *ncout)
{
  //printf("ncout %d %d %d\n",*ncout, (*ncout)%(*nfact), *nfact);
  if ( (*ncout)%(*nfact) == 0 )
    {
      if ( ndotout == MAX_DOTS )
	{
	  *ncout = 0;
	  ndotout = 0;
	  fprintf(stdout, "\n   ");
	  (*nfact) *= 10;
	}

      fprintf(stdout, ".");
      fflush(stdout);
      ndotout++;
    }

  (*ncout)++;

  return ndotout;
}

static
void printTimesteps(int streamID, int taxisID, int verbose)
{
  int nrecs;
  struct datetime {
    int vdate;
    int vtime;
    struct datetime *next;
  };
  struct datetime vdatetime[NUM_TIMESTEP];
  struct datetime *next_vdatetime = vdatetime;

  for ( int i = 0; i < NUM_TIMESTEP-1; ++i ) vdatetime[i].next = &vdatetime[i+1];
  vdatetime[NUM_TIMESTEP-1].next = &vdatetime[0];

  int ntimeout = 0;
  int ndotout = 0;
  int nvdatetime = 0;
  int ncout = 0;
  int nfact = 1;
  int tsID = 0;

#ifdef CDO
  dtlist_type *dtlist = dtlist_new();
#endif
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
#ifdef CDO
      dtlist_taxisInqTimestep(dtlist, taxisID, 0);
      int vdate = dtlist_get_vdate(dtlist, 0);
      int vtime = dtlist_get_vtime(dtlist, 0);
#else
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);
#endif

      if ( verbose || tsID < NUM_TIMESTEP )
	{
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	}
      else
	{
	  if ( tsID == 2*NUM_TIMESTEP ) fprintf(stdout, "\n   ");
	  if ( tsID >= 2*NUM_TIMESTEP ) ndotout = printDot(ndotout, &nfact, &ncout);

	  if ( nvdatetime < NUM_TIMESTEP )
	    {
	      vdatetime[nvdatetime].vdate = vdate;
	      vdatetime[nvdatetime].vtime = vtime;
	      nvdatetime++;
	    }
	  else
	    {
	      next_vdatetime->vdate = vdate;
	      next_vdatetime->vtime = vtime;
	      next_vdatetime = next_vdatetime->next;
	    }
	}

      tsID++;
    }

#ifdef CDO
  dtlist_delete(dtlist);
#endif
  if ( nvdatetime )
    {
      fprintf(stdout, "\n");

      ntimeout = 0;
      int toff = 0;
      if ( tsID > 2*NUM_TIMESTEP )
        {
          toff = tsID%4;
          if ( toff > 0 ) toff = 4 - toff;
          for ( int i = 0; i < toff; ++i ) next_vdatetime = next_vdatetime->next;
        }
      for ( int i = toff; i < nvdatetime; ++i )
	{
	  int vdate = next_vdatetime->vdate;
	  int vtime = next_vdatetime->vtime;
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	  next_vdatetime = next_vdatetime->next;
	}
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
