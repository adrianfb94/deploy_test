#include "cdo_int.h"
#include "grid.h"
#include "remap.h"


static
void calc_bin_addr(size_t gridsize, size_t nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, size_t* restrict bin_addr)
{
  size_t n2;

  for ( size_t n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      bin_addr[n2  ] = gridsize;
      bin_addr[n2+1] = 0;
    }

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  private(n2)  shared(gridsize, nbins, bin_lats, bin_addr, cell_bound_box)
#endif
  for ( size_t nele = 0; nele < gridsize; ++nele )
    {
      size_t nele4 = nele<<2;
      restr_t cell_bound_box_lat1 = cell_bound_box[nele4  ];
      restr_t cell_bound_box_lat2 = cell_bound_box[nele4+1];
      for ( size_t n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  if ( cell_bound_box_lat1 <= bin_lats[n2+1] &&
	       cell_bound_box_lat2 >= bin_lats[n2  ] )
	    {
#if defined(_OPENMP)
#pragma omp critical
#endif
              {
                bin_addr[n2  ] = MIN(nele, bin_addr[n2  ]);
                bin_addr[n2+1] = MAX(nele, bin_addr[n2+1]);
              }
            }
	}
    }
}


void calc_lat_bins(remapgrid_t* src_grid, remapgrid_t* tgt_grid, int map_type)
{
  size_t nbins;
  size_t n;      /* Loop counter                  */
  size_t n2;
  double dlat;                /* lat/lon intervals for search bins  */
  restr_t *bin_lats = NULL;

  nbins = src_grid->num_srch_bins;
  dlat = PI/nbins;

  if ( cdoVerbose ) cdoPrint("Using %d latitude bins to restrict search.", nbins);

  if ( nbins > 0 )
    {
      bin_lats = src_grid->bin_lats = (restr_t*) Realloc(src_grid->bin_lats, 2*nbins*sizeof(restr_t));

      for ( n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  bin_lats[n2  ] = RESTR_SCALE((n  )*dlat - PIH);
	  bin_lats[n2+1] = RESTR_SCALE((n+1)*dlat - PIH);
	}

      src_grid->bin_addr = (size_t*) Realloc(src_grid->bin_addr, 2*nbins*sizeof(size_t));

      calc_bin_addr(src_grid->size, nbins, bin_lats, src_grid->cell_bound_box, src_grid->bin_addr);

      if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSERV_YAC )
	{
	  tgt_grid->bin_addr = (size_t*) Realloc(tgt_grid->bin_addr, 2*nbins*sizeof(size_t));

	  calc_bin_addr(tgt_grid->size, nbins, bin_lats, tgt_grid->cell_bound_box, tgt_grid->bin_addr);

	  Free(src_grid->bin_lats); src_grid->bin_lats = NULL;
	}
   }

  if ( map_type == MAP_TYPE_CONSERV_YAC )
    {
      Free(tgt_grid->cell_bound_box); tgt_grid->cell_bound_box = NULL;
    }
 
  if ( map_type == MAP_TYPE_DISTWGT )
    {
      Free(src_grid->cell_bound_box); src_grid->cell_bound_box = NULL;
    }
}


size_t get_srch_cells(size_t tgt_cell_add, size_t nbins, size_t *bin_addr1, size_t *bin_addr2,
                      restr_t *tgt_cell_bound_box, restr_t *src_cell_bound_box, size_t src_grid_size, size_t *srch_add)
{
  size_t num_srch_cells;  /* num cells in restricted search arrays   */
  size_t min_add;         /* addresses for restricting search of     */
  size_t max_add;         /* destination grid                        */
  size_t n, n2;           /* generic counters                        */
  size_t src_cell_add;    /* current linear address for src cell     */
  size_t src_cell_addm4;
  restr_t bound_box_lat1, bound_box_lat2, bound_box_lon1, bound_box_lon2;

  /* Restrict searches first using search bins */

  min_add = src_grid_size - 1;
  max_add = 0;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( tgt_cell_add >= bin_addr1[n2] && tgt_cell_add <= bin_addr1[n2+1] )
	{
	  if ( bin_addr2[n2  ] < min_add ) min_add = bin_addr2[n2  ];
	  if ( bin_addr2[n2+1] > max_add ) max_add = bin_addr2[n2+1];
	}
    }

  /* Further restrict searches using bounding boxes */

  bound_box_lat1 = tgt_cell_bound_box[0];
  bound_box_lat2 = tgt_cell_bound_box[1];
  bound_box_lon1 = tgt_cell_bound_box[2];
  bound_box_lon2 = tgt_cell_bound_box[3];

  num_srch_cells = 0;
  for ( src_cell_add = min_add; src_cell_add <= max_add; ++src_cell_add )
    {
      src_cell_addm4 = src_cell_add<<2;
      if ( (src_cell_bound_box[src_cell_addm4+2] <= bound_box_lon2)  &&
	   (src_cell_bound_box[src_cell_addm4+3] >= bound_box_lon1) )
	{
	  if ( (src_cell_bound_box[src_cell_addm4  ] <= bound_box_lat2)  &&
	       (src_cell_bound_box[src_cell_addm4+1] >= bound_box_lat1) )
	    {
	      srch_add[num_srch_cells] = src_cell_add;
	      num_srch_cells++;
	    }
	}
    }

  if ( bound_box_lon1 < RESTR_SCALE(0.) || bound_box_lon2 > RESTR_SCALE(PI2) )
    {
      if ( bound_box_lon1 < RESTR_SCALE(0.) )
	{
	  bound_box_lon1 += RESTR_SCALE(PI2);
	  bound_box_lon2 += RESTR_SCALE(PI2);
	}
      else
	{
	  bound_box_lon1 -= RESTR_SCALE(PI2);
	  bound_box_lon2 -= RESTR_SCALE(PI2);
	}

      for ( src_cell_add = min_add; src_cell_add <= max_add; ++src_cell_add )
	{
	  src_cell_addm4 = src_cell_add<<2;
	  if ( (src_cell_bound_box[src_cell_addm4+2] <= bound_box_lon2)  &&
	       (src_cell_bound_box[src_cell_addm4+3] >= bound_box_lon1) )
	    {
	      if ( (src_cell_bound_box[src_cell_addm4  ] <= bound_box_lat2)  &&
		   (src_cell_bound_box[src_cell_addm4+1] >= bound_box_lat1) )
		{
		  size_t ii;
		  for ( ii = 0; ii < num_srch_cells; ++ii )
		    if ( srch_add[ii] == src_cell_add ) break;
		  
		  if ( ii == num_srch_cells )
		    {
		      srch_add[num_srch_cells] = src_cell_add;
		      num_srch_cells++;
		    }
		}
	    }
	}
    }

  return num_srch_cells;
}

static
int grid_search_nn(size_t min_add, size_t max_add, size_t *restrict nbr_add, double *restrict nbr_dist, 
		   double plat, double plon,
		   const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  size_t n, srch_add;
  size_t i;
  double dist_min, distance; /* For computing dist-weighted avg */
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      distance = acos(coslat_dst*cos(src_center_lat[srch_add])*
		     (coslon_dst*cos(src_center_lon[srch_add]) +
                      sinlon_dst*sin(src_center_lon[srch_add]))+
		      sinlat_dst*sin(src_center_lat[srch_add]));

      if ( distance < dist_min )
	{
          for ( n = 0; n < 4; ++n )
	    {
	      if ( distance < nbr_dist[n] )
		{
		  for ( i = 3; i > n; --i )
		    {
		      nbr_add [i] = nbr_add [i-1];
		      nbr_dist[i] = nbr_dist[i-1];
		    }
		  search_result = -1;
		  nbr_add [n] = srch_add;
		  nbr_dist[n] = distance;
		  dist_min = nbr_dist[3];
		  break;
		}
	    }
        }
    }

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return search_result;
}


int grid_search(remapgrid_t *src_grid, size_t *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		const double *restrict src_center_lat, const double *restrict src_center_lon,
		const restr_t *restrict src_grid_bound_box, const size_t *restrict src_bin_add)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double src_lats[4]             ! latitudes  of the four corner points
    double src_lons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center 
    double src_center_lon[]        ! longitude of each src grid center

    restr_t src_grid_bound_box[][4] ! bound box for source grid

    int src_bin_add[][2]           ! latitude bins for restricting
  */
  /*  Local variables */
  size_t n, n2, next_n, srch_add, srch_add4;    /* dummy indices                    */
  /* Vectors for cross-product check */
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon;
  int scross[4], scross_last = 0;
  int search_result = 0;
  restr_t *bin_lats = src_grid->bin_lats;

  size_t nbins = src_grid->num_srch_bins;

  restr_t rlat = RESTR_SCALE(plat);
  restr_t rlon = RESTR_SCALE(plon);

  // restrict search first using bins

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  // addresses for restricting search
  size_t min_add = src_grid->size-1;
  size_t max_add = 0;

  for ( size_t  n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( rlat >= bin_lats[n2] && rlat <= bin_lats[n2+1] )
	{
	  if ( src_bin_add[n2  ] < min_add ) min_add = src_bin_add[n2  ];
	  if ( src_bin_add[n2+1] > max_add ) max_add = src_bin_add[n2+1];
	}
    }
 
  /* Now perform a more detailed search */

  size_t nx = src_grid_dims[0];
  size_t ny = src_grid_dims[1];

  /* srch_loop */
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      srch_add4 = srch_add<<2;
      /* First check bounding box */
      if ( rlon >= src_grid_bound_box[srch_add4+2] &&
	   rlon <= src_grid_bound_box[srch_add4+3] &&
	   rlat >= src_grid_bound_box[srch_add4  ] &&
	   rlat <= src_grid_bound_box[srch_add4+1])
	{
	  /* We are within bounding box so get really serious */

          /* Determine neighbor addresses */
          size_t j = srch_add/nx;
          size_t i = srch_add - j*nx;

          size_t ip1 = (i < (nx-1)) ? i + 1 : (src_grid->is_cyclic) ? 0 : i;
          size_t jp1 = (j < (ny-1)) ? j + 1 : j;

          size_t n_add  = jp1*nx + i;
          size_t e_add  = j  *nx + ip1;
	  size_t ne_add = jp1*nx + ip1;

          src_lons[0] = src_center_lon[srch_add];
          src_lons[1] = src_center_lon[e_add];
          src_lons[2] = src_center_lon[ne_add];
          src_lons[3] = src_center_lon[n_add];

          src_lats[0] = src_center_lat[srch_add];
          src_lats[1] = src_center_lat[e_add];
          src_lats[2] = src_center_lat[ne_add];
          src_lats[3] = src_center_lat[n_add];

	  /* For consistency, we must make sure all lons are in same 2pi interval */

          vec1_lon = src_lons[0] - plon;
          if      ( vec1_lon >  PI ) src_lons[0] -= PI2;
          else if ( vec1_lon < -PI ) src_lons[0] += PI2;

          for ( n = 1; n < 4; ++n )
	    {
	      vec1_lon = src_lons[n] - src_lons[0];
	      if      ( vec1_lon >  PI ) src_lons[n] -= PI2;
	      else if ( vec1_lon < -PI ) src_lons[n] += PI2;
	    }

          /* corner_loop */
          for ( n = 0; n < 4; ++n )
	    {
	      next_n = (n+1)%4;

	      /*
		Here we take the cross product of the vector making up each box side 
                with the vector formed by the vertex and search point.
                If all the cross products are positive, the point is contained in the box.
	      */
	      vec1_lat = src_lats[next_n] - src_lats[n];
	      vec1_lon = src_lons[next_n] - src_lons[n];
	      vec2_lat = plat - src_lats[n];
	      vec2_lon = plon - src_lons[n];

	      /* Check for 0,2pi crossings */

	      if      ( vec1_lon >  THREE*PIH ) vec1_lon -= PI2;
	      else if ( vec1_lon < -THREE*PIH ) vec1_lon += PI2;

	      if      ( vec2_lon >  THREE*PIH ) vec2_lon -= PI2;
	      else if ( vec2_lon < -THREE*PIH ) vec2_lon += PI2;

	      double cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;

	      /* If cross product is less than ZERO, this cell doesn't work    */
	      /* 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero */
	      scross[n] = cross_product < 0 ? -1 : cross_product > 0 ? 1 : 0;

	      if ( n == 0 ) scross_last = scross[n];

	      if ( (scross[n] < 0 && scross_last > 0) || (scross[n] > 0 && scross_last < 0) ) break;

	      scross_last = scross[n];
	    } /* corner_loop */

	  if ( n >= 4 )
	    {
	      n = 0;
	      if      ( scross[0]>=0 && scross[1]>=0 && scross[2]>=0 && scross[3]>=0 ) n = 4;
	      else if ( scross[0]<=0 && scross[1]<=0 && scross[2]<=0 && scross[3]<=0 ) n = 4;
	    }

	  /* If cross products all same sign, we found the location */
          if ( n >= 4 )
	    {
	      src_add[0] = srch_add;
	      src_add[1] = e_add;
	      src_add[2] = ne_add;
	      src_add[3] = n_add;

	      search_result = 1;

	      return search_result;
	    }

	  /* Otherwise move on to next cell */

        } /* Bounding box check */
    } /* srch_loop */

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return search_result;

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_nn(min_add, max_add, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return search_result;
}  /* grid_search */
