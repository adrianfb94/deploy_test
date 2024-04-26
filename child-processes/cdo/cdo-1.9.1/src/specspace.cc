#include <cdi.h>
#include "cdo_int.h"
#include "specspace.h"
#include "grid.h"


void geninx(long ntr, double *f, double *g);
void scaluv(double *fu, double *rclat, int nlat, int lot);
void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt);
void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
	   int nt, int nsp, int nlev);

void after_legini_full(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict pdev,
		       double *restrict pol2, double *restrict pol3, double *restrict coslat);
void after_legini(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict coslat);


void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDout);
  int nlon  = gridInqXsize(gridIDin);
  int nlat  = gridInqYsize(gridIDin);
  int waves = ntr + 1;
  int nfc   = waves * 2;

  double *fpwork = (double*) Malloc(nlat*nfc*nlev*sizeof(double));

  gp2fc(sptrans->trig, sptrans->ifax, arrayIn, fpwork, nlat, nlon, nlev, nfc);
  fc2sp(fpwork, arrayOut, sptrans->pold, nlev, nlat, nfc, ntr);

  Free(fpwork);
}
	   
   
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDin);
  int nlon  = gridInqXsize(gridIDout);
  int nlat  = gridInqYsize(gridIDout);
  int waves = ntr + 1;
  int nfc   = waves * 2;

  double *fpwork = (double*) Malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(arrayIn, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, arrayOut, nlat, nlon, nlev, nfc);

  Free(fpwork);
}


void four2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  (void)gridIDin;
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDout);
  int nlat  = sptrans->nlat;
  int waves = ntr + 1;
  int nfc   = waves * 2;

  fc2sp(arrayIn, arrayOut, sptrans->pold, nlev, nlat, nfc, ntr);
}


void spec2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDin);
  size_t nfc = gridInqSize(gridIDout);
  int nlat  = nfc_to_nlat(nfc, ntr);
  int waves = ntr + 1;
  nfc   = waves * 2;

  sp2fc(arrayIn, arrayOut, sptrans->poli, nlev, nlat, nfc, ntr);
}


void four2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDin);
  int nlon  = gridInqXsize(gridIDout);
  int nlat  = gridInqYsize(gridIDout);
  int waves = ntr + 1;
  int nfc   = waves * 2;

  fc2gp(sptrans->trig, sptrans->ifax, arrayIn, arrayOut, nlat, nlon, nlev, nfc);
}


void grid2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int nlev  = 1;
  int ntr   = gridInqTrunc(gridIDout);
  int nlon  = gridInqXsize(gridIDin);
  int nlat  = gridInqYsize(gridIDin);
  int waves = ntr + 1;
  int nfc   = waves * 2;

  gp2fc(sptrans->trig, sptrans->ifax, arrayIn, arrayOut, nlat, nlon, nlev, nfc);
}


void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int ntrIn  = gridInqTrunc(gridIDin);
  int ntrOut = gridInqTrunc(gridIDout);

  sp2sp(arrayIn, ntrIn, arrayOut, ntrOut);
}


void speccut(int gridIDin, double *arrayIn, double *arrayOut, int *waves)
{
  int ntr = gridInqTrunc(gridIDin);

  spcut(arrayIn, arrayOut, ntr, waves);
}


SPTRANS *sptrans_new(int nlon, int nlat, int ntr, int flag)
{
  SPTRANS *sptrans = (SPTRANS*) Malloc(sizeof(SPTRANS));

  sptrans->nlon = nlon;
  sptrans->nlat = nlat;
  sptrans->ntr  = ntr;

  int nsp = (ntr + 1)*(ntr + 2);
  sptrans->poldim = nsp / 2 * nlat;

  sptrans->trig = (double*) Malloc(nlon * sizeof(double));
  fft_set(sptrans->trig, sptrans->ifax, nlon);

  sptrans->poli = (double*) Malloc(sptrans->poldim * sizeof(double));
  sptrans->pold = (double*) Malloc(sptrans->poldim * sizeof(double));
  
  if ( flag )
    {
      sptrans->pol2 = (double*) Malloc(sptrans->poldim * sizeof(double));
      sptrans->pol3 = (double*) Malloc(sptrans->poldim * sizeof(double));
    }
  else
    {
      sptrans->pol2 = NULL;
      sptrans->pol3 = NULL;
    }

  sptrans->coslat  = (double*) Malloc(nlat * sizeof(double));
  sptrans->rcoslat = (double*) Malloc(nlat * sizeof(double));

  if ( flag )
    after_legini_full(ntr, nlat, sptrans->poli, sptrans->pold, NULL,
		      sptrans->pol2, sptrans->pol3, sptrans->coslat);
  else
    after_legini(ntr, nlat, sptrans->poli, sptrans->pold, sptrans->coslat);

  for ( int jgl = 0; jgl < nlat; ++jgl )
    sptrans->rcoslat[jgl] = 1.0 / sptrans->coslat[jgl];

  return sptrans;
}


void sptrans_delete(SPTRANS *sptrans)
{
  if ( sptrans )
    {
      if ( sptrans->trig ) { Free(sptrans->trig);  sptrans->trig = NULL; }
      if ( sptrans->poli ) { Free(sptrans->poli);  sptrans->poli = NULL; }
      if ( sptrans->pold ) { Free(sptrans->pold);  sptrans->pold = NULL; }
      if ( sptrans->pol2 ) { Free(sptrans->pol2);  sptrans->pol2 = NULL; }
      if ( sptrans->pol3 ) { Free(sptrans->pol3);  sptrans->pol3 = NULL; }
      if ( sptrans->coslat  ) { Free(sptrans->coslat);   sptrans->coslat = NULL; }
      if ( sptrans->rcoslat ) { Free(sptrans->rcoslat);  sptrans->rcoslat = NULL; }

      Free(sptrans); sptrans = NULL;
    }
}


DVTRANS *dvtrans_new(int ntr)
{
  DVTRANS *dvtrans = (DVTRANS*) Malloc(sizeof(DVTRANS));

  dvtrans->ntr = ntr;

  int dimsp = (ntr + 1)*(ntr + 2);
  dvtrans->fdim = dimsp / 2;

  dvtrans->f1 = (double*) Malloc(dvtrans->fdim * sizeof(double));
  dvtrans->f2 = (double*) Malloc(dvtrans->fdim * sizeof(double));

  geninx(ntr, dvtrans->f1, dvtrans->f2);

  return dvtrans;
}


void dvtrans_delete(DVTRANS *dvtrans)
{
  if ( dvtrans )
    {
      if ( dvtrans->f1 ) { Free(dvtrans->f1);  dvtrans->f1 = NULL; }
      if ( dvtrans->f2 ) { Free(dvtrans->f2);  dvtrans->f2 = NULL; }

      Free(dvtrans); dvtrans = NULL;
    }
}


void trans_uv2dv(SPTRANS *sptrans, int nlev,
		 int gridID1, double *gu, double *gv,
		 int gridID2, double *sd, double *svo)
{
  if ( gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("unexpected grid1 type: %s instead of Gaussian", gridNamePtr(gridInqType(gridID1)));

  if ( gridInqType(gridID2) != GRID_SPECTRAL )
    cdoAbort("unexpected grid2 type: %s instead of spectral", gridNamePtr(gridInqType(gridID2)));
    
  int ntr   = gridInqTrunc(gridID2);
  int nlon  = gridInqXsize(gridID1);
  int nlat  = gridInqYsize(gridID1);
  int waves = ntr + 1;
  int nfc   = waves * 2;

  double *fpwork1 = (double*) Malloc(nlat*nfc*nlev*sizeof(double));
  double *fpwork2 = (double*) Malloc(nlat*nfc*nlev*sizeof(double));

  gp2fc(sptrans->trig, sptrans->ifax, gu, fpwork1, nlat, nlon, nlev, nfc);
  gp2fc(sptrans->trig, sptrans->ifax, gv, fpwork2, nlat, nlon, nlev, nfc);

  scaluv(fpwork1, sptrans->coslat, nlat, nfc*nlev);
  scaluv(fpwork2, sptrans->coslat, nlat, nfc*nlev);

  uv2dv(fpwork1, fpwork2, sd, svo, sptrans->pol2, sptrans->pol3, nlev, nlat, ntr);

  Free(fpwork1);
  Free(fpwork2);
}


void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev,
		 int gridID1, double *sd, double *svo,
		 int gridID2, double *gu, double *gv)
{
  if ( gridInqType(gridID1) != GRID_SPECTRAL )
    Warning("unexpected grid1 type: %s", gridNamePtr(gridInqType(gridID1)));

  if ( gridInqType(gridID2) != GRID_GAUSSIAN )
    Warning("unexpected grid2 type: %s", gridNamePtr(gridInqType(gridID2)));

  int ntr   = gridInqTrunc(gridID1);
  int nlon  = gridInqXsize(gridID2);
  int nlat  = gridInqYsize(gridID2);
  int waves = ntr + 1;
  int nfc   = waves * 2;
  int dimsp = (ntr + 1)*(ntr + 2);

  double *su = gu;
  double *sv = gv;

  dv2uv(sd, svo, su, sv, dvtrans->f1, dvtrans->f2, ntr, dimsp, nlev);

  double *fpwork = (double*) Malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(su, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gu, nlat, nlon, nlev, nfc);

  sp2fc(sv, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gv, nlat, nlon, nlev, nfc);

  Free(fpwork);
}
