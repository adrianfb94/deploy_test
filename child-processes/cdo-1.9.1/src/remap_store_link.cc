#include "cdo.h"
#include "cdo_int.h"
#include "remap.h"
#include "remap_store_link.h"

static
int cmp_adds(const void *s1, const void *s2)
{
  int cmp = 0;
  const addweight_t* c1 = (const addweight_t*) s1;
  const addweight_t* c2 = (const addweight_t*) s2;

  if      ( c1->add < c2->add ) cmp = -1;
  else if ( c1->add > c2->add ) cmp =  1;

  return cmp;
}

static
int cmp_adds4(const void *s1, const void *s2)
{
  int cmp = 0;
  const addweight4_t* c1 = (const addweight4_t*) s1;
  const addweight4_t* c2 = (const addweight4_t*) s2;

  if      ( c1->add < c2->add ) cmp = -1;
  else if ( c1->add > c2->add ) cmp =  1;

  return cmp;
}

static
void sort_addweights(size_t num_weights, addweight_t *addweights)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( addweights[n].add < addweights[n-1].add ) break;
  if ( n == num_weights ) return;

  qsort(addweights, num_weights, sizeof(addweight_t), cmp_adds);
}

static
void sort_addweights4(size_t num_weights, addweight4_t *addweights)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( addweights[n].add < addweights[n-1].add ) break;
  if ( n == num_weights ) return;

  qsort(addweights, num_weights, sizeof(addweight4_t), cmp_adds);
}


void sort_add_and_wgts(size_t num_weights, size_t *src_add, double *wgts)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == num_weights ) return;

  if ( num_weights )
    {
      std::vector<addweight_t> addweights(num_weights);

      for ( n = 0; n < num_weights; ++n )
        {
          addweights[n].add    = src_add[n];
          addweights[n].weight = wgts[n];
        }

      qsort(addweights.data(), num_weights, sizeof(addweight_t), cmp_adds);

      for ( n = 0; n < num_weights; ++n )
        {
          src_add[n] = addweights[n].add;
          wgts[n]    = addweights[n].weight;
        }  
    }
}


void sort_add_and_wgts4(size_t num_weights, size_t *src_add, double wgts[4][4])
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == num_weights ) return;

  if ( num_weights )
    {
      std::vector<addweight4_t> addweights(num_weights);

      for ( n = 0; n < num_weights; ++n )
        {
          addweights[n].add = src_add[n];
          for ( unsigned k = 0; k < 4; ++k )
            addweights[n].weight[k] = wgts[n][k];
        }

      qsort(addweights.data(), num_weights, sizeof(addweight4_t), cmp_adds4);

      for ( n = 0; n < num_weights; ++n )
        {
          src_add[n] = addweights[n].add;
          for ( unsigned k = 0; k < 4; ++k )
            wgts[n][k] = addweights[n].weight[k];
        }
    }
}


void store_weightlinks(int lalloc, size_t num_weights, size_t *srch_add, double *weights, size_t cell_add, weightlinks_t *weightlinks)
{
  weightlinks[cell_add].nlinks = 0;
  weightlinks[cell_add].offset = 0;

  if ( num_weights )
    {
      addweight_t *addweights = NULL;
      if ( lalloc )
        addweights = (addweight_t *) Malloc(num_weights*sizeof(addweight_t));
      else
        addweights = weightlinks[cell_add].addweights;

      for ( size_t n = 0; n < num_weights; ++n )
	{
	  addweights[n].add    = srch_add[n];
	  addweights[n].weight = weights[n];
	}

      if ( num_weights > 1 ) sort_addweights(num_weights, addweights);

      weightlinks[cell_add].nlinks = num_weights;

      if ( lalloc ) weightlinks[cell_add].addweights = addweights;
    }
}


void store_weightlinks4(size_t num_weights, size_t *srch_add, double weights[4][4], size_t cell_add, weightlinks4_t *weightlinks)
{
  weightlinks[cell_add].nlinks = 0;
  weightlinks[cell_add].offset = 0;

  if ( num_weights )
    {
      addweight4_t *addweights = weightlinks[cell_add].addweights;

      for ( size_t n = 0; n < num_weights; ++n )
	{
	  addweights[n].add = srch_add[n];
	  for ( unsigned k = 0; k < 4; ++k )
	    addweights[n].weight[k] = weights[n][k];
	}

      sort_addweights4(num_weights, addweights);

      weightlinks[cell_add].nlinks = num_weights;
    }
}


void weightlinks2remaplinks(int lalloc, size_t tgt_grid_size, weightlinks_t *weightlinks, remapvars_t *rv)
{
  size_t nlinks = 0;

  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      if ( weightlinks[tgt_cell_add].nlinks )
	{
	  weightlinks[tgt_cell_add].offset = nlinks;
	  nlinks += weightlinks[tgt_cell_add].nlinks;
	}
    }

  rv->max_links = nlinks;
  rv->num_links = nlinks;
  if ( nlinks )
    {
      rv->src_cell_add = (size_t*) Malloc(nlinks*sizeof(size_t));
      rv->tgt_cell_add = (size_t*) Malloc(nlinks*sizeof(size_t));
      rv->wts          = (double*) Malloc(nlinks*sizeof(double));
      size_t *restrict src_cell_adds = rv->src_cell_add;
      size_t *restrict tgt_cell_adds = rv->tgt_cell_add;
      double *restrict wts = rv->wts;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(none) shared(src_cell_adds,tgt_cell_adds,wts,weightlinks,tgt_grid_size) 
#endif
      for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
	{
	  size_t num_links = weightlinks[tgt_cell_add].nlinks;
	  if ( num_links )
	    {
	      size_t offset = weightlinks[tgt_cell_add].offset;
              addweight_t *addweights = weightlinks[tgt_cell_add].addweights;
 	      for ( size_t ilink = 0; ilink < num_links; ++ilink )
		{
		  src_cell_adds[offset+ilink] = addweights[ilink].add;
		  tgt_cell_adds[offset+ilink] = tgt_cell_add;
		  wts[offset+ilink] = addweights[ilink].weight;
		}
            }
	}

      if ( lalloc )
        {
          for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
            {
              size_t num_links = weightlinks[tgt_cell_add].nlinks;
              if ( num_links ) Free(weightlinks[tgt_cell_add].addweights);
            }
        }
      else
        {
          Free(weightlinks[0].addweights);
        }
    }
}


void weightlinks2remaplinks4(size_t tgt_grid_size, weightlinks4_t *weightlinks, remapvars_t *rv)
{
  size_t nlinks = 0;

  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      if ( weightlinks[tgt_cell_add].nlinks )
	{
	  weightlinks[tgt_cell_add].offset = nlinks;
	  nlinks += weightlinks[tgt_cell_add].nlinks;
	}
    }

  rv->max_links = nlinks;
  rv->num_links = nlinks;
  if ( nlinks )
    {
      rv->src_cell_add = (size_t*) Malloc(nlinks*sizeof(size_t));
      rv->tgt_cell_add = (size_t*) Malloc(nlinks*sizeof(size_t));
      rv->wts          = (double*) Malloc(4*nlinks*sizeof(double));
      size_t *restrict src_cell_adds = rv->src_cell_add;
      size_t *restrict tgt_cell_adds = rv->tgt_cell_add;
      double *restrict wts = rv->wts;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(src_cell_adds,tgt_cell_adds,wts,weightlinks,tgt_grid_size)
#endif
      for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
	{
	  size_t num_links = weightlinks[tgt_cell_add].nlinks;
	  if ( num_links )
	    {
	      size_t offset = weightlinks[tgt_cell_add].offset;
	      addweight4_t *addweights = weightlinks[tgt_cell_add].addweights;
	      for ( size_t ilink = 0; ilink < num_links; ++ilink )
		{
		  src_cell_adds[offset+ilink] = addweights[ilink].add;
		  tgt_cell_adds[offset+ilink] = tgt_cell_add;
		  for ( size_t k = 0; k < 4; ++k )
		    wts[(offset+ilink)*4+k] = addweights[ilink].weight[k];
		}
            }
	}

      Free(weightlinks[0].addweights);
    }
}
