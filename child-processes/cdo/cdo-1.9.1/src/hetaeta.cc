/*
#define  NGP    100000
#define  OUTPUT 1
#define  OPOINT 34179
*/

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include "cdo_int.h"
#include "hetaeta.h"

static const double apr  = 101325.0;   /* reference pressure */
static const double aipr = 1.0/101325.0;

static const double p_firef = 40000.0; /* pressure of reference geopotential */

static const double epsilon = 0.622;
static const double rair    = 287.04;
static const double cpair   = 1004.6;

static const double eta_pbl = 0.8;     /* upper limit of BPL eta coordiante */

#if defined(OUTPUT)
static FILE *old, *new;
#endif


static
long int_index(long n, double *x1, double x2)
{
  long k;
  long klo = 0;
  long khi = n-1;

  while ( khi-klo > 1 )
    {
      k = (khi+klo)/2;
      if (x1[k] > x2)
	khi = k;
      else
	klo = k;
    }

  /*
  for ( klo = 0; klo < n-1; klo++ )
    if ( x2 >= x1[klo] && x2 < x1[klo+1] ) break;
  */
  /* printf("%d %d %g %g %g\n", klo, khi, x1[klo], x1[klo+1], x2); */

  return klo;
}

static
double esat(double temperature)
{
  double zes;
  double tc = temperature-273.16;

  if ( tc < 0.0 )
    zes = 21.8745584*tc / (temperature-7.66);
  else
    zes = 17.2693882*tc / (temperature-35.86);

  double es = 610.78*exp(zes);

  return es;
}


#define MAX_VARS  512

/* Source from INTERA */

static
void hetaeta_sc(bool ltq, int lpsmod, long ij, long ngp, long nlev1, long nlev2, long nvars,
		const double *restrict af1, const double *restrict bf1, const double *restrict etah2,
		const double *restrict af2, const double *restrict bf2, const double *restrict w1, 
		const double *restrict w2, const long *restrict jl1, const long *restrict jl2,
		const double *restrict ah1, const double *restrict bh1, const double *restrict ps1, 
		double epsm1i, const double *restrict q1, const double *restrict t1, 
		const double *restrict fis1, const double *restrict fis2, double *restrict ps2,
		const double *restrict ah2, const double *restrict bh2, double *restrict *restrict vars1,
		double *restrict *restrict vars2, double *restrict t2, double *restrict q2,
		double *restrict tscor, double *restrict pscor, double *restrict secor, long jblt,
		double *restrict ph1, double *restrict lnph1, double *restrict fi1, double *restrict pf1,
		double *restrict lnpf1, double *restrict tv1, double *restrict theta1, double *restrict rh1,
		double *restrict zvar, double *restrict ph2, double *restrict lnph2, double *restrict fi2,
		double *restrict pf2, double *restrict rh2, double *restrict wgt, long *restrict idx,
		double *restrict rh_pbl, double *restrict theta_pbl, double *restrict *restrict vars_pbl, 
		double *restrict zt2, double *restrict zq2)
{
  long ijk, ijk1, ijk2;
  long jlev = 0, jlevr, jnop;
  long klo;
  long jjblt;
  double zq1, zt1;
  double dfi, fiadj = 0, dteta = 0;

  double rair_d_cpair = rair/cpair;

  long nlev1p1 = nlev1+1;
  long nlev2p1 = nlev2+1;

  /******* initialise atmospheric fields in old system */
  
  /* pressure */
  ph1[0]   =  0.0;
  lnph1[0] = -1.0; 
  for ( int k = 1; k < nlev1p1; ++k )
    {
      ph1[k]   = ah1[k]+bh1[k]*ps1[ij];
      lnph1[k] = log(ph1[k]);
    } 

  for ( int k = 0; k < nlev1; ++k )
    {
      pf1[k]   = af1[k]+bf1[k]*ps1[ij];
      lnpf1[k] = log(pf1[k]);
    }

  /* virtual temperature, relative humidity, potential temperature */
  if ( ltq )
    for ( int k = 0; k < nlev1; ++k )
      {
	ijk = k*ngp+ij;
	zq1 = q1[ijk];
	zt1 = t1[ijk];
	tv1[k]    = (1.0+epsm1i*zq1)*zt1;
	rh1[k]    = zq1*pf1[k]/(epsilon*esat(zt1));
	theta1[k] = zt1*pow(apr/pf1[k], rair_d_cpair);
      }

  /* ****** integrate hydrostatic equation, using interpolated orography */
  if ( ltq )
    {
      fi1[0] = 0.0;
      fi1[nlev1] = fis1[ij];
      for ( int k = nlev1-1; k > 0; --k )
	{
	  fi1[k] = fi1[k+1]+rair*tv1[k]*(lnph1[k+1]-lnph1[k]);
	}
    }
#if defined(OUTPUT)
  if ( ij == OPOINT )
    for ( int k = nlev1-1; k >= 0; --k )
      { 
	ijk = k*ngp+ij;
	if ( ltq ) { t = t1[ijk]; q = q1[ijk]; fi = fi1[k]; }
	else       { t = 0; q = 0; fi = 0; }
	fprintf(old, "%3d %18.10f %18.10f %18.10f %18.10f", k, fi/g, pf1[k], t, q);
	for ( int iv = 0; iv < nvars; ++iv ) fprintf(old, " %18.10f", vars1[iv][ijk]);
	fprintf(old, "\n");
      }
#endif

  /******* find new surface pressure
	   extra-/interpolate to new orography
	   linear regression works not well for extrapolation
	   separation necessary
  */
  if ( ltq )
    {
      for ( int k = nlev1-2; k > 0; --k )
	{
	  /* find index for regression, 1 <= jlev <= nlevec-2 */
	  jlev = k;
	  if (fis2[ij] < fi1[k]) break;
	}

      double zdff = fi1[jlev+1]-fis2[ij];
      
      /* get the number of points used for estimation of regression coefficients */
      jlevr = 0;
      for ( int k = jlev-1; k > 0; --k )
	{
	  jlevr = k;
	  double zdffl = fi1[k]-fi1[jlev+1];
	  if ( zdffl >= zdff ) break;
	}

      jnop = jlev+1-jlevr+1;

      /*
	get coefficients of regression between Tv and lnP ::Tv = B*lnP + C
	using three levels surounding new orography geopotential
      */
      double zsumt  = 0.0;
      double zsump  = 0.0;
      double zsumpp = 0.0;
      double zsumtp = 0.0;

      for ( int k = jlevr; k <= jlev+1; ++k )
	{
	  zsumt  = zsumt  + tv1[k];
	  zsump  = zsump  + lnpf1[k];
	  zsumpp = zsumpp + lnpf1[k]*lnpf1[k];
	  zsumtp = zsumtp + tv1[k]*lnpf1[k];
	}

      /* final regression coefficients */
      double zb = jnop*zsumpp - zsump*zsump;
      double zc = (zsumt*zsumpp-zsump*zsumtp)/zb;
      zb = (jnop*zsumtp-zsump*zsumt)/zb;

      /* calculate preliminary surface pressure, adjust to middle level */
      double zps = lnph1[jlev];

      /* calculate preliminary pressure */
      if ( fabs(zb) < 1.0e-20 )
	{
	  /* constant virtual temperature near new surface */
	  ps2[ij] = exp(zps+(fi1[jlev]-fis2[ij])/(zc*rair));
	}
      else
	{
	  /* virtual temperatur not constant near new surface */
	  double zbb = zc*zc + zb*(zps*(zb*zps+2.0*zc)+2.0*(fi1[jlev]-fis2[ij])/rair);
	  ps2[ij] = exp((sqrt(zbb)-zc)/zb);
	}
    }
  else
    {
      ps2[ij] = ps1[ij];
    }


  ph2[0]   =  0.0;
  lnph2[0] = -1.0; 
  for ( int k = 1; k < nlev2p1; ++k )
    {
      ph2[k]   = ah2[k]+bh2[k]* ps2[ij];
      lnph2[k] = log(ph2[k]);
    } 

  for ( int k = 0; k < nlev2; ++k )
    {
      pf2[k]   = af2[k]+bf2[k]* ps2[ij];
      /* lnpf2[k] = log(pf2[k]); */
    }

  /******* find reference geopotential,  */

  if ( lpsmod && ltq )
    {
      /* using old pressure at half levels
	 find first level below reference pressure */
      for ( int k = 1; k < nlev1p1; ++k )
	{
	  jlev = k;
	  if ( ph1[k] > p_firef ) break;
	}
	  
      fiadj = fi1[jlev]+(fi1[jlev-1]-fi1[jlev])*
	log(ph1[jlev]/p_firef)/log(ph1[jlev]/ph1[jlev-1]);
    }

  /******* find the new boundary layer top */

  /* using the pressure from the old system */
  double pbl_lim = ps1[ij]*eta_pbl;
  jjblt = nlev2-1;
  for ( int k = nlev2-1; k > 0; --k )
    {
      /* find the next upper level in new system */
      double pbl_lim_need = ps2[ij] *etah2[k];
      if (pbl_lim > pbl_lim_need) break;
      jjblt = jjblt-1;
    }

  /* correct the merging level */
  if ( jblt < jjblt ) jjblt = jblt;

  /******* PBL profile interpolation */
  /* tension spline interpolation with full eta levels */

  if ( ltq )
    for ( int k = jjblt; k < nlev2; ++k )
      {
	theta_pbl[k] = w1[k]*theta1[jl1[k]] + w2[k]*theta1[jl2[k]];
	rh_pbl[k]    = w1[k]*rh1[jl1[k]]    + w2[k]*rh1[jl2[k]];
      }

  for ( int iv = 0; iv < nvars; ++iv )
    for ( int k = jjblt; k < nlev2; ++k )
      {
	ijk1 = jl1[k]*ngp+ij;
	ijk2 = jl2[k]*ngp+ij;
	vars_pbl[iv][k] = w1[k]*vars1[iv][ijk1] + w2[k]*vars1[iv][ijk2];
      }

  /******* linear interpolation using pressure in free atmosphere
	   pressure in new system using preliminary pressure */

  for ( int k = 0; k <= jjblt; ++k )
    {
      idx[k] = int_index(nlev1, pf1, pf2[k]);
    }

  for ( int k = 0; k <= jjblt; ++k )
    {
      wgt[k] = (pf1[idx[k]+1]-pf2[k])/(pf1[idx[k]+1]-pf1[idx[k]]);
    }

  if ( ltq )
    {
      for ( int k = 0; k < nlev1; ++k )
	{
	  ijk = k*ngp+ij;
	  zvar[k] = t1[ijk];
	}

      for ( int k = 0; k <= jjblt; ++k )
	{
	  klo = idx[k];
	  zt2[k] = wgt[k]*zvar[klo] + (1-wgt[k])*zvar[klo+1];
	  rh2[k] = wgt[k]*rh1[klo]  + (1-wgt[k])*rh1[klo+1];
	}
    }

  for ( int iv = 0; iv < nvars; ++iv )
    {
      for ( int k = 0; k < nlev1; ++k )
	{
	  ijk = k*ngp+ij;
	  zvar[k] = vars1[iv][ijk];
	}
      for ( int k = 0; k <= jjblt; ++k )
	{
	  ijk = k*ngp+ij;
	  klo = idx[k];
	  vars2[iv][ijk] = wgt[k]*zvar[klo] + (1-wgt[k])*zvar[klo+1];
	}
    }

  /******* merge boundary layer and free atmosphere */

  if ( ltq )
    {
      /* correction of potential temperature at top of PBL */
      dteta = zt2[jjblt]*pow(apr/pf2[jjblt], rair_d_cpair)-theta_pbl[jjblt];
	  
      /* merge top layer values */
      rh2[jjblt] = 0.5*(rh2[jjblt]+rh_pbl[jjblt]);
    }

  ijk = jjblt*ngp+ij;
  for ( int iv = 0; iv < nvars; ++iv )
    {
      vars2[iv][ijk] = 0.5*(vars2[iv][ijk]+vars_pbl[iv][jjblt]);
    }

  /* correct boundary profile values */
  if ( ltq )
    for ( int k = jjblt+1; k < nlev2; ++k ) 
      {
	zt2[k] = (theta_pbl[k]+dteta)*pow(pf2[k]/apr, rair_d_cpair);
	rh2[k] = rh_pbl[k];
      }

  for ( int iv = 0; iv < nvars; ++iv )
    for ( int k = jjblt+1; k < nlev2; ++k ) 
      {
	ijk = k*ngp+ij;
	vars2[iv][ijk] = vars_pbl[iv][k];
      }

  if ( ltq )
    for ( int k = 0; k < nlev2; ++k )
      {
	zq2[k] = rh2[k]*epsilon*esat(zt2[k])/pf2[k];
      }

  /******* reference level correction */
  if ( lpsmod && ltq )
    {
      /* integrate hydrostatic equation with preliminary temperature and pressure */
      fi2[nlev2] = fis2[ij];
      fi2[0]     = -1.0; /* top not defined, infinity */
    
      /* problem at top level, top pressure is zero per definition */
      for ( int k = nlev2-1; k > 0; --k )
	{
	  fi2[k] = fi2[k+1]+rair*zt2[k]*(lnph2[k+1]-lnph2[k])*(1.0+epsm1i*zq2[k]);
	}

      /* search next level above reference level in new system */
      for ( int k = nlev2-1; k > 0; --k )
	{ 
	  jlev = k; 
	  if (ph2[k] < p_firef) break;
	}
	  
      /* correct surface pressure */
      dfi = fiadj-(fi2[jlev+1]+(fi2[jlev]-fi2[jlev+1])* 
		   log(ph2[jlev+1]/p_firef)/log(ph2[jlev+1]/ph2[jlev]));
      double ztv = (1.0+epsm1i*zq2[nlev2-1])*zt2[nlev2-1];
      ps2[ij] = ps2[ij] *exp(dfi/(rair*ztv));
    }

  /******* final calculation of specific humidity profiles */
  if ( ltq )
    {
      for ( int k = 0; k < nlev2; ++k )
	pf2[k] = af2[k] + bf2[k]* ps2[ij];

      for ( int k = 0; k < nlev2; ++k )
	zq2[k] = rh2[k]*epsilon*esat(zt2[k])/pf2[k];
    }

#if defined(OUTPUT)
  if ( ij == OPOINT )
    for ( int k = nlev2-1; k >= 0; --k )
      { 
	ijk = k*ngp+ij;
	if ( ltq ) { t = t2[k]; q = zq2[k]; fi = fi2[k]; }
	else       { t = 0; q = 0; fi = 0; }
	fprintf(new, "%3d %18.10f %18.10f %18.10f %18.10f", k, fi/g, pf2[k], t, q);
	for ( int iv = 0; iv < nvars; ++iv ) fprintf(new, " %18.10f", vars2[iv][ijk]);
	fprintf(new, "\n");
      }
#endif

  if ( ltq )
    {
      /* calculate surface temperature correction (old version) */
      tscor[ij] = dteta*pow(ps2[ij]/apr, rair_d_cpair);
      pscor[ij] = pow(ps2[ij]/ps1[ij], rair_d_cpair);

      /* correction term of static energy of lowest layer */
      secor[ij] = tv1[nlev1-1]*(cpair+rair*(1.0-ph1[nlev1-1]
					    /(ps1[ij]-ph1[nlev1-1])*log(ps1[ij]/ph1[nlev1-1])));
    }

  if ( ltq )
    for ( int k = 0; k < nlev2; ++k ) 
      {
	ijk = k*ngp+ij;
	t2[ijk] = zt2[k];
	q2[ijk] = zq2[k];
      }
}


void hetaeta(bool ltq, int ngp, const int *imiss,
	     int nlev1, const double *ah1, const double *bh1,
             const double *fis1, const double *ps1, 
             const double *t1, const double *q1,
             int nlev2, const double *ah2, const double *bh2, 
             const double *fis2, double *ps2, 
             double *t2, double *q2,
	     int nvars, double **vars1, double **vars2,
	     double *tscor, double *pscor, double *secor)
{
  double epsm1i;
  long jblt;
  long jlev = 0;
  int lpsmod = 1;
  double *vars_pbl[MAX_VARS];
#if defined(OUTPUT)
  double t, q, fi;

  old = fopen("old.dat","w");
  new = fopen("new.dat","w");
#endif

  long nlev1p1 = nlev1+1;
  long nlev2p1 = nlev2+1;

  NEW_2D(double, ph1, ompNumThreads, nlev1p1);
  NEW_2D(double, lnph1, ompNumThreads, nlev1p1);
  NEW_2D(double, fi1, ompNumThreads, nlev1p1);

  NEW_2D(double, pf1, ompNumThreads, nlev1);
  NEW_2D(double, lnpf1, ompNumThreads, nlev1);
  NEW_2D(double, tv1, ompNumThreads, nlev1);
  NEW_2D(double, theta1, ompNumThreads, nlev1);
  NEW_2D(double, rh1, ompNumThreads, nlev1);
  NEW_2D(double, zvar, ompNumThreads, nlev1);

  NEW_2D(double, ph2, ompNumThreads, nlev2p1);
  NEW_2D(double, lnph2, ompNumThreads, nlev2p1);
  NEW_2D(double, fi2, ompNumThreads, nlev2p1);

  NEW_2D(double, pf2, ompNumThreads, nlev2);
  NEW_2D(double, rh2, ompNumThreads, nlev2);
  NEW_2D(double, wgt, ompNumThreads, nlev2);
  NEW_2D(long, idx, ompNumThreads, nlev2);

  NEW_2D(double, zt2, ompNumThreads, ltq?nlev2:0);
  NEW_2D(double, zq2, ompNumThreads, ltq?nlev2:0);
  NEW_2D(double, rh_pbl, ompNumThreads, ltq?nlev2:0);
  NEW_2D(double, theta_pbl, ompNumThreads, ltq?nlev2:0);

  if ( nvars > MAX_VARS )
    {
      fprintf(stderr, "Too many vars (max = %d)!\n", MAX_VARS);
      exit(-1);
    }
  // if ( nvars > 0 )
  double ***vars_pbl_2 = (double***) Malloc(ompNumThreads*sizeof(double**));

  if ( nvars > 0 )
    {
      for ( int i = 0; i < ompNumThreads; i++ )
        {
	  vars_pbl_2[i]  = (double **) Malloc(nvars*sizeof(double *));
	  for ( int iv = 0; iv < nvars; ++iv )
	    vars_pbl_2[i][iv] = (double*) Malloc(nlev2*sizeof(double));
	}
    }
  
  double *af1 = new double[nlev1];
  double *bf1 = new double[nlev1];
  double *etaf1 = new double[nlev1];

  double *etah2 = new double[nlev2p1];

  double *af2 = new double[nlev2];
  double *bf2 = new double[nlev2];
  double *etaf2 = new double[nlev2];

  double *w1 = new double[nlev2];
  double *w2 = new double[nlev2];
  long *jl1 = new long[nlev2];
  long *jl2 = new long[nlev2];

  /******* set coordinate system ETA's, A's, B's
	   calculate half and full level ETA
	   set the boundary layer index */

  /* input system */

  for ( int k = 0; k < nlev1; ++k )
    {
      af1[k] = 0.5*(ah1[k]+ah1[k+1]);
      bf1[k] = 0.5*(bh1[k]+bh1[k+1]);
    }

  // etah1[nlev1] = ah1[nlev1]*aipr+bh1[nlev1];
  for ( int k = 0; k < nlev1; ++k )
    { 
      // etah1[k] = ah1[k]*aipr+bh1[k];
      etaf1[k] = af1[k]*aipr+bf1[k];
    }

  /* output system */

  /* calculates full level VCT */
  for ( int k = 0; k < nlev2; ++k )
    {
      af2[k] = 0.5*(ah2[k]+ah2[k+1]);
      bf2[k] = 0.5*(bh2[k]+bh2[k+1]);
    }

  etah2[nlev2] = ah2[nlev2]*aipr+bh2[nlev2];
  jblt = nlev2;
  for ( int k = nlev2-1; k >= 0; --k )
    { 
      etah2[k] = ah2[k]*aipr+bh2[k];
      etaf2[k] = af2[k]*aipr+bf2[k];
      if (etah2[k] > eta_pbl) jblt = k;
    }

  /* calculate weights for PBL interpolation */
  for ( int k = 0; k < nlev2; ++k )
    {
      /* scan through new vertical levels
	 set changes outside the full level eta's of old system to constant */
      if ( etaf2[k] <= etaf1[0] )
	{
	  /* at top of atmosphere */
	  jl1[k] = 0;
	  jl2[k] = 1;
	  w2[k]  = 0.0;
	}
      else if ( etaf2[k] >= etaf1[nlev1-1])
	{
	  /* at surface of atmosphere */
	  jl1[k] = nlev1-2;
	  jl2[k] = nlev1-1;
	  w2[k]  = 1.0;
	}
      else
	{
	  for ( jlev = nlev1-2; jlev >= 1; jlev-- )
	    {
	      jl1[k] = jlev; /* find nearest eta level below */
	      if (etaf2[k] > etaf1[jlev]) break;
	    }
	  jl2[k] = jl1[k]+1;
	  w2[k]  = log(etaf2[k]/etaf1[jl1[k]]) / log(etaf1[jl2[k]]/etaf1[jl1[k]]);
	}
      w1[k] = 1.0-w2[k];
    }

  epsm1i = 1.0/epsilon-1.0;


#if defined(_OPENMP)
#pragma omp parallel for default(none) firstprivate(lpsmod)  private(vars_pbl)  schedule(dynamic,1) \
  shared(ngp, ph1, lnph1, fi1, pf1, lnpf1, tv1, theta1, rh1, zvar, ph2, lnph2, \
	 fi2, pf2, rh_pbl, zt2, zq2, theta_pbl, rh2, wgt, idx, vars_pbl_2, \
	 af1, bf1, etah2, af2, bf2, w1, w2, jl1, jl2,	\
	 ltq, nvars, imiss, ah1, bh1, ps1, nlev1, epsm1i, q1, t1, fis1, fis2, ps2, \
	 ah2, bh2, nlev2, vars1, vars2, t2, q2, tscor, pscor, secor, jblt)
#endif
  for ( int ij = 0; ij < ngp; ++ij )
    {
      int ompthID = cdo_omp_get_thread_num();

      for ( int iv = 0; iv < nvars; ++iv )
        vars_pbl[iv] = vars_pbl_2[ompthID][iv];

      if ( imiss && imiss[ij] ) continue;

      hetaeta_sc(ltq, lpsmod, ij, ngp, nlev1, nlev2, nvars,
		 af1, bf1, etah2, af2, bf2, w1, w2, jl1, jl2,
		 ah1, bh1, ps1, epsm1i, q1, t1, fis1, fis2, ps2,
		 ah2, bh2, vars1, vars2, t2, q2, tscor, pscor, secor, jblt,
		 ph1[ompthID], lnph1[ompthID], fi1[ompthID], pf1[ompthID], lnpf1[ompthID], tv1[ompthID], theta1[ompthID], rh1[ompthID], zvar[ompthID], 
                 ph2[ompthID], lnph2[ompthID], fi2[ompthID], pf2[ompthID], rh2[ompthID], wgt[ompthID], idx[ompthID],
		 rh_pbl[ompthID], theta_pbl[ompthID], vars_pbl, zt2[ompthID], zq2[ompthID]);

    } /* end for ij */

#if defined(OUTPUT)
  fclose(old);
  fclose(new);
#endif

  DELETE_2D(ph1);
  DELETE_2D(lnph1);
  DELETE_2D(fi1);

  DELETE_2D(pf1);
  DELETE_2D(lnpf1);
  DELETE_2D(tv1);
  DELETE_2D(theta1);
  DELETE_2D(rh1);
  DELETE_2D(zvar);

  DELETE_2D(ph2);
  DELETE_2D(lnph2);
  DELETE_2D(fi2);

  DELETE_2D(pf2);
  // DELETE_2D(lnpf2);
  DELETE_2D(rh2);
  DELETE_2D(wgt);
  DELETE_2D(idx);

  DELETE_2D(zt2); 
  DELETE_2D(zq2); 
  DELETE_2D(rh_pbl); 
  DELETE_2D(theta_pbl); 

  if ( nvars > 0 )
    {
      for ( int i = 0; i < ompNumThreads; i++ )
        {
	  for ( int iv = 0; iv < nvars; ++iv )
	    Free(vars_pbl_2[i][iv]);
	  
	  Free(vars_pbl_2[i]);
	}
    }
  Free(vars_pbl_2);

  delete[] af1;
  delete[] bf1;
  delete[] etaf1;
  delete[] etah2;
  delete[] af2;
  delete[] bf2;
  delete[] etaf2;
  delete[] w1;
  delete[] w2;
  delete[] jl1;
  delete[] jl2;

  return;
}

#if defined(TEST_HETAETA)
#define NGP  (512*1024)
int ompNumThreads = 2;

int main (int argc, char *argv[])
{
  printf("NGP = %d\n", NGP);

  double a2[41] = {
       0.00000000000000000,       2000.00000000000000000,       4000.00000000000000000,
    6000.00000000000000000,       8000.00000000000000000,       9976.13671875000000000,
   11902.14453125000000000,      13722.03125000000000000,      15379.80468750000000000,
   16819.47265625000000000,      18045.18359375000000000,      19027.69531250000000000,
   19755.10937500000000000,      20222.20312500000000000,      20429.86328125000000000,
   20384.48046875000000000,      20097.40234375000000000,      19584.32812500000000000,
   18864.75000000000000000,      17961.35937500000000000,      16899.46875000000000000,
   15706.44921875000000000,      14411.12500000000000000,      13043.21875000000000000,
   11632.75781250000000000,      10209.50000000000000000,       8802.35546875000000000,
    7438.80468750000000000,       6144.31640625000000000,       4941.77734375000000000,
    3850.91333007812500000,       2887.69653320312500000,       2063.77978515625000000,
    1385.91259765625000000,        855.36181640625000000,        467.33349609375000000,
     210.39390563964843750,         65.88919067382812500,          7.36769962310791016,
       0.00000000000000000,          0.00000000000000000 };

  double b2[41] = {
    0.00000000000000000,    0.00000000000000000,    0.00000000000000000,
    0.00000000000000000,    0.00000000000000000,    0.00039085815660655,
    0.00182679994031787,    0.00513499975204468,    0.01114289835095406,
    0.02067789807915688,    0.03412120044231415,    0.05169039964675903,
    0.07353377342224121,    0.09967470169067383,    0.13002246618270874,
    0.16438430547714233,    0.20247590541839600,    0.24393308162689209,
    0.28832298517227173,    0.33515489101409912,    0.38389205932617188,
    0.43396288156509399,    0.48477149009704590,    0.53570991754531860,
    0.58616840839385986,    0.63554751873016357,    0.68326860666275024,
    0.72878581285476685,    0.77159661054611206,    0.81125342845916748,
    0.84737491607666016,    0.87965691089630127,    0.90788388252258301,
    0.93194031715393066,    0.95182150602340698,    0.96764522790908813,
    0.97966271638870239,    0.98827010393142700,    0.99401938915252686,
    0.99763011932373047,    1.00000000000000000 };

  double a1[20] = {
       0.00000000000000000,       2000.00000000000000000,       4000.00000000000000000,   
    6046.10937500000000000,       8267.92968750000000000,      10609.51171875000000000,   
   12851.10156250000000000,      14698.50000000000000000,      15861.12890625000000000,   
   16116.23828125000000000,      15356.92187500000000000,      13621.46093750000000000,  
   11101.55859375000000000,       8127.14453125000000000,       5125.14062500000000000,   
    2549.96899414062500000,        783.19506835937500000,          0.00000000000000000,   
       0.00000000000000000,          0.00000000000000000 };  

  double b1[20] = {
    0.00000000000000000,    0.00000000000000000,    0.00000000000000000,
    0.00033899326808751,    0.00335718691349030,    0.01307003945112228,
    0.03407714888453484,    0.07064980268478394,    0.12591671943664551,
    0.20119541883468628,    0.29551959037780762,    0.40540921688079834,
    0.52493220567703247,    0.64610791206359863,    0.75969839096069336,
    0.85643762350082397,    0.92874687910079956,    0.97298520803451538,
    0.99228149652481079,    1.00000000000000000 };

  double iu1[19] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  double iv1[19] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  double it1[19] = {
    224.2257,  212.7505,  209.553,  208.7785,  211.7619,  220.2336,  221.2698,
    220.3876,  227.1461,  237.6735, 248.1776,  258.1013,  264.4792,  269.1322,
    271.9017,  275.6761,  279.819,  282.2512,  284.141 };

  double iq1[19] = {
    2.512447e-06,  2.176736e-06,  2.170464e-06,  2.01653e-06,  1.805185e-06,
    1.726813e-06,  3.75322e-06,   8.901303e-06,  3.285719e-05, 0.0001270178,
    0.0003347051,  0.0007223329,  0.001228461,   0.001733165,  0.002967748,
    0.004558741,   0.004706143,   0.004668835,   0.004677606 };

  double icl1[19] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  -4.987459e-40,  -4.791847e-39,  -3.970467e-23,
    -1.902515e-23,  -1.694066e-21,  -3.705769e-22,  -1.799945e-21,  -4.632211e-22,
    2.072752e-05,  0.000149563,  -1.482308e-20,  -2.541099e-21,  5.033612e-05 };

  double ici1[19] = {
    -4.408104e-37, 0.0,  0.0,  -2.003328e-25,  -9.305782e-24,  -2.15067e-23,  -9.926167e-23,
    -1.958764e-21,  -8.735027e-22,  -2.779327e-22,  -2.117582e-21,  -1.323489e-21,  -8.470329e-22,
    -4.102816e-22,  -1.429368e-21,  -2.646978e-21,  -5.029258e-22,  -8.205632e-22,  -1.588187e-21 };

  double icc1[19] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.4445496,  0.0,  0.0,  0.001098633} ;

  double ifis1 = 9.121094, ips1 = 102511.8;
  
  double *vars1[5];
  double *vars2[5];

  int ij, k;
  bool ltq = true;

  double *fis1 = (double*) Malloc(NGP*sizeof(double));
  double *ps1  = (double*) Malloc(NGP*sizeof(double));
  double *fis2 = (double*) Malloc(NGP*sizeof(double));
  double *ps2  = (double*) Malloc(NGP*sizeof(double));

  double *tscor  = (double*) Malloc(NGP*sizeof(double));
  double *pscor  = (double*) Malloc(NGP*sizeof(double));
  double *secor  = (double*) Malloc(NGP*sizeof(double));

  double *t1  = (double*) Malloc(NGP*19*sizeof(double));
  double *q1  = (double*) Malloc(NGP*19*sizeof(double));
  double *u1  = (double*) Malloc(NGP*19*sizeof(double));
  double *v1  = (double*) Malloc(NGP*19*sizeof(double));
  double *cl1 = (double*) Malloc(NGP*19*sizeof(double));
  double *ci1 = (double*) Malloc(NGP*19*sizeof(double));
  double *cc1 = (double*) Malloc(NGP*19*sizeof(double));

  double *t2  = (double*) Malloc(NGP*40*sizeof(double));
  double *q2  = (double*) Malloc(NGP*40*sizeof(double));
  double *u2  = (double*) Malloc(NGP*40*sizeof(double));
  double *v2  = (double*) Malloc(NGP*40*sizeof(double));
  double *cl2 = (double*) Malloc(NGP*40*sizeof(double));
  double *ci2 = (double*) Malloc(NGP*40*sizeof(double));
  double *cc2 = (double*) Malloc(NGP*40*sizeof(double));

  for ( ij = 0; ij < NGP; ++ij )
    {
      ps1[ij]  = ips1;
      fis1[ij] = ifis1;
    }

  for ( k = 0; k < 19; ++k )
    for ( ij = 0; ij < NGP; ++ij )
      {
	t1[k*NGP+ij] = it1[k];
	q1[k*NGP+ij] = iq1[k];
	u1[k*NGP+ij] = iu1[k];
	v1[k*NGP+ij] = iv1[k];
	cl1[k*NGP+ij] = icl1[k];
	ci1[k*NGP+ij] = ici1[k];
	cc1[k*NGP+ij] = icc1[k];
      }

  vars1[0] = u1;
  vars1[1] = v1;
  vars1[2] = cl1;
  vars1[3] = ci1;
  vars1[4] = cc1;

  vars2[0] = u2;
  vars2[1] = v2;
  vars2[2] = cl2;
  vars2[3] = ci2;
  vars2[4] = cc2;

  for ( ij = 0; ij < NGP; ++ij ) fis2[ij] = fis1[ij];

  if ( ltq )
    hetaeta(ltq, NGP, NULL,
	    19, a1, b1,
	    fis1, ps1,
	    t1, q1,
	    40, a2, b2,
	    fis2, ps2,
	    t2, q2,
	    5, vars1, vars2,
	    tscor, pscor, secor);
  else
    hetaeta(ltq, NGP, NULL,
	    19, a1, b1,
	    fis1, ps1,
	    NULL, NULL,
	    40, a2, b2,
	    fis2, ps2,
	    NULL, NULL,
	    5, vars1, vars2,
	    NULL, NULL, NULL);

  double ot2[40] = {
    224.226, 212.75, 209.616, 208.895, 210.468, 213.454, 218.091, 220.49, 220.977, 221.114,
    220.731, 220.803, 223.595, 226.547, 230.566, 234.99, 239.411, 243.722, 248.125, 252.164,
    256.221, 259.499, 262.082, 264.584, 266.487, 268.319, 269.761, 270.879, 271.881, 273.501,
    274.932, 276.493, 278.199, 279.638, 280.811, 281.746, 282.613, 283.479, 284.087, 284.343};

  int lerror = 0;
  for ( k = 0; k < 40; ++k )
    if ( fabs(t2[k*NGP]-ot2[k]) > 0.001 ) lerror=1;

  if ( lerror )
    for ( k = 0; k < 40; ++k )
      printf("%d %g %g %g\n", k, t2[k*NGP], ot2[k], fabs(t2[k*NGP]-ot2[k]));

  return 0;
}
#endif
