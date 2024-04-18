#if defined(_OPENMP)
#include <omp.h> // omp_get_wtime
#endif

#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"
#include "grid_search.h"

//  Interpolation using a distance-weighted average

static
void nbr_store_distance(size_t nadd, double distance, size_t num_neighbors, size_t *restrict nbr_add, double *restrict nbr_dist)
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
void nbr_check_distance(size_t num_neighbors, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // If distance is zero, set to small number
  for ( size_t nchk = 0; nchk < num_neighbors; ++nchk )
    if ( nbr_add[nchk] < ULONG_MAX && nbr_dist[nchk] <= 0. ) nbr_dist[nchk] = TINY;
}


double nbr_compute_weights(size_t num_neighbors, const int *restrict src_grid_mask, bool *restrict nbr_mask, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // Compute weights based on inverse distance if mask is false, eliminate those points

  double dist_tot = 0.; // sum of neighbor distances (for normalizing)

  if ( src_grid_mask )
    {
      for ( size_t n = 0; n < num_neighbors; ++n )
        {
          nbr_mask[n] = false;
          if ( nbr_add[n] < ULONG_MAX )
            if ( src_grid_mask[nbr_add[n]] )
              {
                nbr_dist[n] = 1./nbr_dist[n];
                dist_tot += nbr_dist[n];
                nbr_mask[n] = true;
              }
        }
    }
  else
    {
      for ( size_t n = 0; n < num_neighbors; ++n )
        {
          nbr_mask[n] = false;
          if ( nbr_add[n] < ULONG_MAX )
            {
              nbr_dist[n] = 1./nbr_dist[n];
              dist_tot += nbr_dist[n];
              nbr_mask[n] = true;
            }
        }
    }
  

  return dist_tot;
}


size_t nbr_normalize_weights(size_t num_neighbors, double dist_tot, const bool *restrict nbr_mask, size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // Normalize weights and store the link

  size_t nadds = 0;

  for ( size_t n = 0; n < num_neighbors; ++n )
    {
      if ( nbr_mask[n] )
        {
          nbr_dist[nadds] = nbr_dist[n]/dist_tot;
          nbr_add[nadds]  = nbr_add[n];
          nadds++;
        }
    }

  return nadds;
}


//  This routine finds the closest num_neighbor points to a search point and computes a distance to each of the neighbors.

#define MAX_SEARCH_CELLS 25
static
void grid_search_nbr_reg2d(struct gridsearch *gs, size_t num_neighbors, remapgrid_t *src_grid, size_t *restrict nbr_add, double *restrict nbr_dist, 
			   double plon, double plat, const int *restrict src_grid_dims)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  size_t n, nadd;
  size_t ii, jj;
  long i, j, ix;
  size_t src_add[MAX_SEARCH_CELLS];
  size_t *src_add_tmp = NULL;
  size_t *psrc_add = src_add;
  size_t num_add = 0;
  double distance;   //  Angular distance
  double cos_search_radius = cos(gs->search_radius);
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
  double *restrict coslon = gs->coslon;
  double *restrict sinlon = gs->sinlon;
  double *restrict coslat = gs->coslat;
  double *restrict sinlat = gs->sinlat;
  double *restrict src_center_lon = gs->reg2d_center_lon;
  double *restrict src_center_lat = gs->reg2d_center_lat;

  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];

  size_t nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( plon < src_center_lon[0]     ) plon += PI2;
  if ( plon > src_center_lon[nxm-1] ) plon -= PI2;

  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      if ( src_grid->is_cyclic && ii == (nxm-1) ) ii = 0;

      long k;
      for ( k = 3; k < 10000; k+=2 )
        if ( num_neighbors <= ((k-2)*(k-2)) ) break;

      if ( (k*k) > MAX_SEARCH_CELLS ) psrc_add = src_add_tmp = (size_t *) Malloc(k*k*sizeof(size_t));

      k /= 2;

      long j0 = jj-k;
      long jn = jj+k;
      long i0 = ii-k;
      long in = ii+k;
      if ( j0 < 0 ) j0 = 0;
      if ( jn >= ny ) jn = ny-1;
      if ( (in-i0) > nx ) { i0 = 0; in = nx-1; }
      
      for ( j = j0; j <= jn; ++j )
	for ( i = i0; i <= in; ++i )
	  {
	    ix = i;
	    
	    if ( src_grid->is_cyclic )
	      {
		if ( ix <   0 ) ix += nx;
		if ( ix >= nx ) ix -= nx;
	      }

	    if ( ix >= 0 && ix < nx && j >= 0 && j < ny )
	      psrc_add[num_add++] = j*nx+ix;
	  }
    }

  // Initialize distance and address arrays
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = ULONG_MAX;
      nbr_dist[n] = BIGNUM;
    }

  if ( lfound )
    {
      size_t ix, iy;

      for ( size_t na = 0; na < num_add; ++na )
	{
	  nadd = psrc_add[na];

	  iy = nadd/nx;
	  ix = nadd - iy*nx;

	  // Find distance to this point
	  distance =  sinlat_dst*sinlat[iy] + coslat_dst*coslat[iy]*
	             (coslon_dst*coslon[ix] + sinlon_dst*sinlon[ix]);

	  // Check that distance is inside the range of -1 to 1, otherwise the result of acos(distance) is NaN
	  if ( distance > 1. ) distance =  1.;

	  if ( distance >= cos_search_radius )
	    {
	      distance = acos(distance);

	      // Store the address and distance if this is one of the smallest four so far
	      nbr_store_distance(nadd, distance, num_neighbors, nbr_add, nbr_dist);
	    }
	}

      nbr_check_distance(num_neighbors, nbr_add, nbr_dist);

      if ( src_add_tmp ) Free(src_add_tmp);
    }
  else if ( src_grid->lextrapolate )
    {
      int search_result;

      if ( num_neighbors < 4 )
	{
	  size_t nbr_add4[4];
	  double nbr_dist4[4];
	  for ( n = 0; n < num_neighbors; ++n ) nbr_add4[n] = ULONG_MAX;
	  search_result = grid_search_reg2d_nn(nx, ny, nbr_add4, nbr_dist4, plat, plon, src_center_lat, src_center_lon);
	  if ( search_result < 0 )
	    {
	      for ( n = 0; n < num_neighbors; ++n ) nbr_add[n]  = nbr_add4[n];
	      for ( n = 0; n < num_neighbors; ++n ) nbr_dist[n] = nbr_dist4[n];
	    }
	}
      else
	{
	  search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
	}

      if ( search_result >= 0 )
	for ( n = 0; n < num_neighbors; ++n ) nbr_add[n] = ULONG_MAX;
    }
} // grid_search_nbr_reg2d


int grid_search_nbr(struct gridsearch *gs, size_t num_neighbors, size_t *restrict nbr_add, double *restrict nbr_dist, double plon, double plat)
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
  for ( size_t n = 0; n < num_neighbors; ++n ) nbr_add[n]  = ULONG_MAX;
  for ( size_t n = 0; n < num_neighbors; ++n ) nbr_dist[n] = BIGNUM;

  size_t ndist = num_neighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if ( ndist > 8 ) ndist += 8;
  else             ndist *= 2; 
  if ( ndist > gs->n ) ndist = gs->n;

  double zdist[32];
  size_t zadds[32];
  double *dist = zdist;
  size_t *adds = zadds;
  if ( num_neighbors > 16 )
    {
      dist = (double*) Malloc(ndist*sizeof(double));
      adds = (size_t*) Malloc(ndist*sizeof(size_t));
    }
  
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
    nbr_store_distance(adds[j], dist[j], max_neighbors, nbr_add, nbr_dist);

  nbr_check_distance(max_neighbors, nbr_add, nbr_dist);

  if ( num_neighbors > 16 )
    {
      Free(dist);
      Free(adds);
    }

  if ( ndist > num_neighbors ) ndist = num_neighbors;

  return ndist;
} // grid_search_nbr

//  This routine computes the inverse-distance weights for a nearest-neighbor interpolation.

void remap_distwgt_weights(size_t num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;
  size_t nx = src_grid->dims[0];
  size_t ny = src_grid->dims[1];
  bool lcyclic = src_grid->is_cyclic;

  weightlinks_t *weightlinks = (weightlinks_t *) Malloc(tgt_grid_size*sizeof(weightlinks_t));
  weightlinks[0].addweights = (addweight_t *) Malloc(num_neighbors*tgt_grid_size*sizeof(addweight_t));
  for ( size_t tgt_cell_add = 1; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    weightlinks[tgt_cell_add].addweights = weightlinks[0].addweights + num_neighbors*tgt_cell_add;

  NEW_2D(bool, nbr_mask, ompNumThreads, num_neighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, num_neighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, num_neighbors); // angular distance four nearest neighbors

#if defined(_OPENMP)
  double start = 0;
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  struct gridsearch *gs = NULL;
  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(lcyclic, nx, ny, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);
  else
    gs = gridsearch_create(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

#if defined(_OPENMP)
  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", omp_get_wtime()-start);
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(gs, weightlinks, num_neighbors, remap_grid_type, src_grid, tgt_grid, tgt_grid_size, findex) \
  shared(nbr_mask, nbr_add, nbr_dist)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      int ompthID = cdo_omp_get_thread_num();

      weightlinks[tgt_cell_add].nlinks = 0;	

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;
	
      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, num_neighbors, src_grid, nbr_add[ompthID], nbr_dist[ompthID], 
			      plon, plat, src_grid->dims);
      else
        grid_search_nbr(gs, num_neighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      double dist_tot = nbr_compute_weights(num_neighbors, src_grid->mask, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      // Normalize weights and store the link
      size_t nadds = nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      for ( size_t n = 0; n < nadds; ++n )
        if ( nbr_mask[ompthID][n] ) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      store_weightlinks(0, nadds, nbr_add[ompthID], nbr_dist[ompthID], tgt_cell_add, weightlinks);
    }

  progressStatus(0, 1, 1);

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

  weightlinks2remaplinks(0, tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) Free(weightlinks);

#if defined(_OPENMP)
  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime()-start);
#endif

} // remap_distwgt_weights

static
void distwgt_remap(double* restrict tgt_point, const double* restrict src_array, size_t nadds, const double wgts[4], const size_t src_add[4])
{
  *tgt_point = 0.;
  for ( size_t n = 0; n < nadds; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
}


void remap_distwgt(size_t num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  int src_remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid 

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;
  size_t nx = src_grid->dims[0];
  size_t ny = src_grid->dims[1];
  bool lcyclic = src_grid->is_cyclic;

  NEW_2D(bool, nbr_mask, ompNumThreads, num_neighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, num_neighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, num_neighbors); // angular distance four nearest neighbors

#if defined(_OPENMP)
  double start = 0;
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  struct gridsearch *gs = NULL;
  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(lcyclic, nx, ny, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);
  else
    gs = gridsearch_create(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

#if defined(_OPENMP)
  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", omp_get_wtime()-start);
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(gs, num_neighbors, src_remap_grid_type, src_grid, tgt_grid, tgt_grid_size, findex) \
  shared(src_array, tgt_array, missval, nbr_mask, nbr_add, nbr_dist)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      int ompthID = cdo_omp_get_thread_num();

      tgt_array[tgt_cell_add] = missval;

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, num_neighbors, src_grid, nbr_add[ompthID], nbr_dist[ompthID], 
			      plon, plat, src_grid->dims);
      else
        grid_search_nbr(gs, num_neighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      
      // Compute weights based on inverse distance if mask is false, eliminate those points
      double dist_tot = nbr_compute_weights(num_neighbors, src_grid->mask, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      // Normalize weights and store the link
      size_t nadds = nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      for ( size_t n = 0; n < nadds; ++n )
        if ( nbr_mask[ompthID][n] ) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      if ( nadds > 1 ) sort_add_and_wgts(nadds, nbr_add[ompthID], nbr_dist[ompthID]);

      if ( nadds ) distwgt_remap(&tgt_array[tgt_cell_add], src_array, nadds, nbr_dist[ompthID], nbr_add[ompthID]);
    }

  progressStatus(0, 1, 1);

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

#if defined(_OPENMP)
  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime()-start);
#endif

} // remap_distwgt
