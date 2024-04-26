#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "cdo_int.h"
#include "dmemory.h"
#include "util.h"
#include "grid.h"
#include "grid_search.h"


#define  PI       M_PI
#define  PI2      (2.0*PI)


static int gridsearch_method_nn = GS_KDTREE;

static
double cdo_default_search_radius(void)
{
  extern double gridsearch_radius;

  double search_radius = gridsearch_radius;

  if ( search_radius <    0. ) search_radius = 0.;
  if ( search_radius >  180. ) search_radius = 180.;

  search_radius = search_radius*DEG2RAD;

  return search_radius;
}

static inline
void LLtoXYZ_f(double lon, double lat, float *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = cos_lat * cos(lon);
   xyz[1] = cos_lat * sin(lon);
   xyz[2] = sin(lat);
}

static inline
void LLtoXYZ_kd(double lon, double lat, kdata_t *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = KDATA_SCALE(cos_lat * cos(lon));
   xyz[1] = KDATA_SCALE(cos_lat * sin(lon));
   xyz[2] = KDATA_SCALE(sin(lat));
}

static constexpr
float square(const float x)
{
  return x*x;
}

static constexpr
float distance(const float *restrict a, const float *restrict b)
{
  return (square((a[0]-b[0]))+square((a[1]-b[1]))+square((a[2]-b[2])));
}


void gridsearch_set_method(const char *methodstr)
{
  if      ( strcmp(methodstr, "kdtree")  == 0 ) gridsearch_method_nn = GS_KDTREE;
  else if ( strcmp(methodstr, "nearpt3") == 0 ) gridsearch_method_nn = GS_NEARPT3;
  else if ( strcmp(methodstr, "full")    == 0 ) gridsearch_method_nn = GS_FULL;
  else
    cdoAbort("gridsearch method %s not available!", methodstr);
}


struct gridsearch *gridsearch_create_reg2d(bool lcyclic, size_t nx, size_t ny, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->nx = nx;
  gs->ny = ny;

  unsigned nxm = nx;
  if ( lcyclic ) nxm++;

  double *reg2d_center_lon = (double *) Malloc(nxm*sizeof(double));
  double *reg2d_center_lat = (double *) Malloc(ny*sizeof(double));

  memcpy(reg2d_center_lon, lons, nxm*sizeof(double));
  memcpy(reg2d_center_lat, lats, ny*sizeof(double));

  double *coslon = (double *) Malloc(nx*sizeof(double));
  double *sinlon = (double *) Malloc(nx*sizeof(double));
  double *coslat = (double *) Malloc(ny*sizeof(double));
  double *sinlat = (double *) Malloc(ny*sizeof(double));

  for ( size_t n = 0; n < nx; ++n )
    {
      double rlon = lons[n];
      if ( rlon > PI2 ) rlon -= PI2;
      if ( rlon < 0   ) rlon += PI2;
      coslon[n] = cos(rlon);
      sinlon[n] = sin(rlon);
    }
  for ( size_t n = 0; n < ny; ++n )
    {
      coslat[n] = cos(lats[n]);
      sinlat[n] = sin(lats[n]);
    }

  gs->reg2d_center_lon = reg2d_center_lon;
  gs->reg2d_center_lat = reg2d_center_lat;

  gs->coslon = coslon;
  gs->sinlon = sinlon;
  gs->coslat = coslat;
  gs->sinlat = sinlat;

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


struct kdNode *gs_create_kdtree(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n * sizeof(struct kd_point));  
  // see  example_cartesian.c
  if ( cdoVerbose ) printf("kdtree lib init 3D: n=%zu  nthreads=%d\n", n, ompNumThreads);
  kdata_t min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
  kdata_t *restrict point;
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ ) 
    {
      point = pointlist[i].point;
      LLtoXYZ_kd(lons[i], lats[i], point);
      for ( size_t j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  struct kdNode *kdt = kd_buildTree(pointlist, n, min, max, 3, ompNumThreads);
  if ( pointlist ) Free(pointlist);

  return kdt;
}


void gs_destroy_nearpt3(struct gsNear *near)
{
  if ( near )
    {
#if defined(ENABLE_NEARPT3)
      if ( near->nearpt3 ) nearpt3_destroy(near->nearpt3);
#endif
      if ( near->pts )
        {
          Free(near->pts[0]);
          Free(near->pts);
        }

      Free(near);
    }
}


struct gsNear *gs_create_nearpt3(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gsNear *near = (struct gsNear *) Calloc(1, sizeof(struct gsNear));

  Coord_T **p = (Coord_T **) Malloc(n*sizeof(Coord_T *));
  p[0] = (Coord_T *) Malloc(3*n*sizeof(Coord_T));
  for ( size_t i = 1; i < n; i++ ) p[i] = p[0] + i*3;

  float point[3];

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ )
    {
      LLtoXYZ_f(lons[i], lats[i], point);

      p[i][0] = NPT3SCALE(point[0]);
      p[i][1] = NPT3SCALE(point[1]);
      p[i][2] = NPT3SCALE(point[2]);
    }

  near->n = n;
  near->plons = lons;
  near->plats = lats;
  near->pts = p;
#if defined(ENABLE_NEARPT3)
  near->nearpt3 = nearpt3_preprocess(n, p);
#else
  cdoAbort("nearpt3 support not compiled in!");
#endif
  
  return near;
}


void gs_destroy_full(struct gsFull *full)
{
  if ( full )
    {
      if ( full->pts )
        {
          Free(full->pts[0]);
          Free(full->pts);
        }

      Free(full);
    }
}


struct gsFull *gs_create_full(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gsFull *full = (struct gsFull *) Calloc(1, sizeof(struct gsFull));

  float **p = (float **) Malloc(n*sizeof(float *));
  p[0] = (float *) Malloc(3*n*sizeof(float));
  for ( size_t i = 1; i < n; i++ ) p[i] = p[0] + i*3;

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ )
    {
      LLtoXYZ_f(lons[i], lats[i], p[i]);
    }
  
  full->n = n;
  full->plons = lons;
  full->plats = lats;
  full->pts = p;

  return full;
}


struct gridsearch *gridsearch_create(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->n = n;
  if ( n == 0 ) return gs;

  gs->kdt  = gs_create_kdtree(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


struct gridsearch *gridsearch_create_nn(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if ( n == 0 ) return gs;

  if      ( gs->method_nn == GS_KDTREE  ) gs->kdt  = gs_create_kdtree(n, lons, lats);
  else if ( gs->method_nn == GS_NEARPT3 ) gs->near = gs_create_nearpt3(n, lons, lats);
  else if ( gs->method_nn == GS_FULL    ) gs->full = gs_create_full(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


void gridsearch_delete(struct gridsearch *gs)
{
  if ( gs )
    {
      if ( gs->kdt ) kd_destroyTree(gs->kdt);
      
      if ( gs->reg2d_center_lon ) Free(gs->reg2d_center_lon);
      if ( gs->reg2d_center_lat ) Free(gs->reg2d_center_lat);

      if ( gs->coslat ) Free(gs->coslat);
      if ( gs->coslon ) Free(gs->coslon);
      if ( gs->sinlat ) Free(gs->sinlat);
      if ( gs->sinlon ) Free(gs->sinlon);

      if ( gs->near ) gs_destroy_nearpt3(gs->near);
      if ( gs->full ) gs_destroy_full(gs->full);

      Free(gs);
    }
}

static
double gs_set_range(double *prange)
{
  double range;

  if ( prange )
    range = *prange;
  else
    range = SQR(2 * M_PI);     /* This has to be bigger than the presumed
                                * maximum distance to the NN but smaller
                                * than once around the sphere. The content
                                * of this variable is replaced with the
                                * distance to the NN squared. */
  return range;
}


kdNode *gs_nearest_kdtree(kdNode *kdt, double lon, double lat, double *prange)
{
  if ( kdt == NULL ) return NULL;
  
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);

  kdata_t point[3];
  LLtoXYZ_kd(lon, lat, point);

  kdNode *node = kd_nearest(kdt, point, &range, 3);

  float frange = KDATA_INVSCALE(range);
  if ( !(frange < range0) ) node = NULL;
  if ( prange ) *prange = frange;

  return node;
}


unsigned gs_nearest_nearpt3(struct gsNear *near, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  if ( near == NULL ) return index;
  
#if defined(ENABLE_NEARPT3)
  float range0 = gs_set_range(prange);

  float point[3];
  LLtoXYZ_f(lon, lat, point);

  Coord_T q[3];
  q[0] = NPT3SCALE(point[0]);
  q[1] = NPT3SCALE(point[1]);
  q[2] = NPT3SCALE(point[2]);

  int closestpt = nearpt3_query(near->nearpt3, q);

  if ( closestpt >= 0 )
    {
      float point0[3];
      LLtoXYZ_f(near->plons[closestpt], near->plats[closestpt], point0);
      
      float range = distance(point, point0);
      if ( range < range0 )
        {
           index = (unsigned) closestpt;
           *prange = range;
        }
    }
#else
  UNUSED(lon);
  UNUSED(lat);
  UNUSED(prange);
#endif

  return index;
}


size_t gs_nearest_full(struct  gsFull *full, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  if ( full == NULL ) return index;
  
  float range0 = gs_set_range(prange);

  float point[3];
  LLtoXYZ_f(lon, lat, point);

  size_t n = full->n;
  float **pts = full->pts;
  size_t closestpt = n;
  float dist = FLT_MAX;
  for ( size_t i = 0; i < n; i++ )
    {
      float d = distance(point, pts[i]);
      if ( closestpt >=n || d < dist || (d<=dist && i < closestpt) )
        {
          dist = d;
          closestpt = i;
        }
    }

  if ( closestpt < n )
    {
      if ( dist < range0 )
        {
          *prange = dist;
          index = closestpt;
        }
    }
  
  return index;
}


size_t gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;

  if ( gs )
    {
      if ( gs->method_nn == GS_KDTREE )
        {
          kdNode *node = gs_nearest_kdtree(gs->kdt, lon, lat, prange);
          if ( node ) index = (int) node->index;
        }
      else if ( gs->method_nn == GS_NEARPT3 )
        {
          index = gs_nearest_nearpt3(gs->near, lon, lat, prange);
        }
      else if ( gs->method_nn == GS_FULL )
        {
          index = gs_nearest_full(gs->full, lon, lat, prange);
        }
      else
        {
          cdoAbort("gridsearch_nearest::method_nn undefined!");
        }
    }

  return index;
}


struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn)
{
  if ( gs->kdt == NULL ) return NULL;
  
  kdata_t point[3];
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);
  struct pqueue *result = NULL;

  LLtoXYZ_kd(lon, lat, point);

  if ( gs )
    {
      result = kd_qnearest(gs->kdt, point, &range, nnn, 3);
      // printf("range %g %g %g %p\n", lon, lat, range, node);

      float frange = KDATA_INVSCALE(range);
      /*
      if ( !(frange < range0) )
        {
          if ( result )
            {
              struct resItem *p;
              while ( pqremove_min(result, &p) ) Free(p); // Free the result node taken from the heap
              Free(result->d); // free the heap
              Free(result);    // and free the heap information structure
            }
          result = NULL;
        }
      */
      if ( prange ) *prange = frange;
    }
  
  return result;
}

#define  BIGNUM   1.e+20
#define  TINY     1.e-14

static
void knn_store_distance(size_t nadd, double distance, size_t num_neighbors, size_t *restrict nbr_add, double *restrict nbr_dist)
{
  if ( num_neighbors == 1 )
    {
      if ( distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]) )
	{
	  nbr_add[0]  = nadd;
	  nbr_dist[0] = distance;
	}
    }
  else
    {
      for ( size_t nchk = 0; nchk < num_neighbors; ++nchk )
	{
	  if ( distance < nbr_dist[nchk] || (distance <= nbr_dist[nchk] && nadd < nbr_add[nchk]) )
	    {
	      for ( size_t n = num_neighbors-1; n > nchk; --n )
		{
		  nbr_add[n]  = nbr_add[n-1];
		  nbr_dist[n] = nbr_dist[n-1];
		}
	      nbr_add[nchk]  = nadd;
	      nbr_dist[nchk] = distance;
	      break;
	    }
	}
    }
}

static
void knn_check_distance(size_t num_neighbors, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // If distance is zero, set to small number
  for ( size_t nchk = 0; nchk < num_neighbors; ++nchk )
    if ( nbr_add[nchk] != GS_NOT_FOUND && nbr_dist[nchk] <= 0. ) nbr_dist[nchk] = TINY;
}


void gridsearch_knn_init(struct gsknn *knn)
{
  unsigned ndist = knn->ndist;
  size_t *restrict add = knn->add;
  double *restrict dist = knn->dist;

  for ( size_t i = 0; i < ndist; ++i )
    {
      add[i]  = GS_NOT_FOUND;
      dist[i] = BIGNUM;
    }
}


struct gsknn *gridsearch_knn_new(size_t size)
{
  struct gsknn *knn = (struct gsknn *) Malloc(sizeof(struct gsknn));
  
  knn->ndist   = size;
  knn->size    = size;
  knn->mask    = (bool*) Malloc(size*sizeof(bool));     // mask at nearest neighbors
  knn->add     = (size_t*) Malloc(size*sizeof(size_t)); // source address at nearest neighbors
  knn->dist    = (double*) Malloc(size*sizeof(double)); // angular distance of the nearest neighbors
  knn->tmpadd  = NULL;
  knn->tmpdist = NULL;

  gridsearch_knn_init(knn);

  return knn;
}


void gridsearch_knn_delete(struct gsknn *knn)
{
  if ( knn )
    {
      knn->size = 0;
      if ( knn->dist    ) Free(knn->dist);
      if ( knn->add     ) Free(knn->add);
      if ( knn->tmpdist ) Free(knn->tmpdist);
      if ( knn->tmpadd  ) Free(knn->tmpadd);
      if ( knn->mask    ) Free(knn->mask);
      Free(knn);
    }
}


size_t gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  double search_radius = gs->search_radius;

  // Initialize distance and address arrays
  gridsearch_knn_init(knn);

  size_t num_neighbors = knn->size;
  size_t *restrict nbr_add = knn->add;
  double *restrict nbr_dist = knn->dist;

  size_t ndist = num_neighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if ( ndist > 8 ) ndist += 8;
  else             ndist *= 2; 
  if ( ndist > gs->n ) ndist = gs->n;

  if ( knn->tmpadd  == NULL ) knn->tmpadd  = (size_t*) Malloc(ndist*sizeof(size_t));
  if ( knn->tmpdist == NULL ) knn->tmpdist = (double*) Malloc(ndist*sizeof(double));

  size_t *adds = knn->tmpadd;
  double *dist = knn->tmpdist;
  
  const double range0 = SQR(search_radius);
  double range = range0;

  size_t j = 0;

  if ( num_neighbors == 1 )
    {
      size_t nadd = gridsearch_nearest(gs, plon, plat, &range);
      if ( nadd != GS_NOT_FOUND )
        {
          //if ( range < range0 )
            {
              dist[j] = sqrt(range);
              adds[j] = nadd;
              j++;
            }
        }
    }
  else
    {
      struct pqueue *gs_result = gridsearch_qnearest(gs, plon, plat, &range, ndist);
      if ( gs_result )
        {
          size_t nadd;
          struct resItem *p;
          while ( pqremove_min(gs_result, &p) )
            {
              nadd  = p->node->index;
              range = p->dist_sq;
              Free(p); // Free the result node taken from the heap

              if ( range < range0 )
                {
                  dist[j] = sqrt(range);
                  adds[j] = nadd;
                  j++;
                }
            }
          Free(gs_result->d); // free the heap
          Free(gs_result);    // and free the heap information structure
        }
    }

  ndist = j;
  size_t max_neighbors = (ndist < num_neighbors) ? ndist : num_neighbors;

  for ( j = 0; j < ndist; ++j )
    knn_store_distance(adds[j], dist[j], max_neighbors, nbr_add, nbr_dist);

  knn_check_distance(max_neighbors, nbr_add, nbr_dist);

  if ( ndist > num_neighbors ) ndist = num_neighbors;

  knn->ndist = ndist;

  return ndist;
} // gridsearch_knn
