/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"


void farfun(field_type *field1, field_type field2, int function)
{
  // clang-format off
  switch (function)
    {
    case func_add:     faradd(field1, field2);   break;
    case func_min:     farmin(field1, field2);   break;
    case func_max:     farmax(field1, field2);   break;
    case func_sum:     farsum(field1, field2);   break;
    case func_mean:    farsum(field1, field2);   break;
    case func_avg:     faradd(field1, field2);   break;
    case func_sub:     farsub(field1, field2);   break;
    case func_mul:     farmul(field1, field2);   break;
    case func_div:     fardiv(field1, field2);   break;
    case func_atan2:   faratan2(field1, field2); break;
    case func_setmiss: farsetmiss(field1, field2); break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}

static
int farsetnmiss(int len, double *restrict array, double missval)
{
  int nmiss = 0;

  if ( DBL_IS_NAN(missval) )
    {
      for ( int i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array[i], missval) || array[i] < 0 )
          {
            array[i] = missval;
            nmiss++;
          }
    }
  else
    {
      for ( int i = 0; i < len; i++ )
        if ( IS_EQUAL(array[i], missval) || array[i] < 0 )
          {
            array[i] = missval;
            nmiss++;
          }
    }

  return nmiss;
}


void farcpy(field_type *field1, field_type field2)
{
  int nwpv = field1->nwpv;
  int gridsize1 = field1->size;
  int gridsize2 = field2.size;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;
  const float *restrict array2f = field2.ptrf;

  if ( gridsize1 == 0 ) gridsize1 = gridInqSize(field1->grid);
  if ( gridsize2 == 0 ) gridsize2 = gridInqSize(field2.grid);

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridsize1;

  if ( len != nwpv*gridsize2 )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( field2.memtype == MEMTYPE_FLOAT )
    for ( int i = 0; i < len; i++ ) array1[i] = array2f[i];
  else
    for ( int i = 0; i < len; i++ ) array1[i] = array2[i];
}


void faradd(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;
  const float *restrict array2f = field2.ptrf;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] = ADDMN(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      if ( field2.memtype == MEMTYPE_FLOAT )
        {
          for ( int i = 0; i < len; i++ ) array1[i] += array2f[i];
        }
      else
        {
          int fpeRaised = array_add_array(len, array1, array2);
        }
    }
}


void farsum(field_type *field1, field_type field2)
{
  int nwpv = field1->nwpv;
  int gridsize1 = field1->size;
  int gridsize2 = field2.size;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;
  const float *restrict array2f = field2.ptrf;

  if ( gridsize1 == 0 ) gridsize1 = gridInqSize(field1->grid);
  if ( gridsize2 == 0 ) gridsize2 = gridInqSize(field2.grid);

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridsize1;

  if ( len != nwpv*gridsize2 )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i];
	    else
	      array1[i] = array2[i];
	  }

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      if ( field2.memtype == MEMTYPE_FLOAT )
        {
          for ( int i = 0; i < len; i++ ) array1[i] += array2f[i];
        }
      else
        {
          array_add_array(len, array1, array2);
        }
    }
}


void farsumw(field_type *field1, field_type field2, double w)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += w*array2[i];
	    else
	      array1[i] = w*array2[i];
	  }

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(array1,array2,w,len)
#endif
      for ( int i = 0; i < len; i++ ) array1[i] += w*array2[i];
    }
}

/* 
 * Compute the occurrence of values in field, if they do not equal refval.
 * This can be used to compute the lengths of multiple periods in a timeseries.
 * Missing field values are handled like refval, i.e. they stop a running
 * period. If there is missing data in the occurence field, missing fields
 * values do not change anything (they do not start a non-period by setting
 * occurrence to zero).
 */
void farsumtr(field_type *occur, field_type field, const double refval)
{
  double omissval = occur->missval;
  double fmissval = field.missval;
  double *restrict oarray = occur->ptr;
  double *restrict farray = field.ptr;

  int len = gridInqSize(occur->grid);

  if ( len != gridInqSize(field.grid) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( occur->nmiss > 0 || field.nmiss > 0 )
    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(farray[i], fmissval) )
	  {
	    if ( !DBL_IS_EQUAL(oarray[i], omissval) )
	      oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
	    else
	      oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : 1.0;
	  }
	else
	{
	  if ( !DBL_IS_EQUAL(oarray[i], omissval) )
	    oarray[i] = 0.0;
	}

      occur->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(oarray[i], omissval) ) occur->nmiss++;
    }
  else
    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
      for ( int i = 0; i < len; i++ ) 
	oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
    }
}


void farsumq(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;
  const float *restrict array2f = field2.ptrf;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i]*array2[i];
	    else
	      array1[i] = array2[i]*array2[i];
	  }

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      if ( field2.memtype == MEMTYPE_FLOAT )
        {
          for ( int i = 0; i < len; i++ ) array1[i] += ((double)array2f[i])*array2f[i];
        }
      else
        {
          for ( int i = 0; i < len; i++ ) array1[i] += array2[i]*array2[i];
        }
    }
}


void farsumqw(field_type *field1, field_type field2, double w)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += w*array2[i]*array2[i];
	    else
	      array1[i] = w*array2[i]*array2[i];
	  }

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] += w*array2[i]*array2[i];
    }
}


void farsub(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] = SUBMN(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] -= array2[i];
    }
}


void farmul(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);
  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] = MULMN(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] *= array2[i];
    }
}


void fardiv(field_type *field1, field_type field2)
{
  int nwpv  = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *array1  = field1->ptr;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);
  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( int i = 0; i < len; i++ ) 
    array1[i] = DIVMN(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void faratan2(field_type *field1, field_type field2)
{
  int nwpv  = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( int i = 0; i < len; i++ ) 
    array1[i] = DBL_IS_EQUAL(array1[i],missval1) || DBL_IS_EQUAL(array2[i],missval2) ? missval1 : atan2(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void farsetmiss(field_type *field1, field_type field2)
{
  int nwpv  = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  int nmiss1 = field1->nmiss;
  double missval1 = field1->missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);
  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 )
    {
      for ( int i = 0; i < len; i++ ) 
        array1[i] = DBL_IS_EQUAL(array1[i],missval1) ? array2[i] : array1[i];

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
}


void farmin(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MIN(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ )
	array1[i] = MIN(array1[i], array2[i]);
    }
}


void farmax(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MAX(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ )
	array1[i] = MAX(array1[i], array2[i]);
    }
}


void farvar(field_type *field1, field_type field2, field_type field3, int divisor)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;
  const double *restrict array3  = field3.ptr;
  double temp;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( (nmiss1 || nmiss2) /*&& (DBL_IS_NAN(missval1) || DBL_IS_NAN(missval2))*/ )
    {
      for ( int i = 0; i < len; i++ )
        {
          temp      = DIV( MULMN(array1[i], array1[i]), array3[i]);
          array1[i] = DIV( SUBMN(array2[i], temp), array3[i]-divisor);
          if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
        }
    }
  else
    {
      for ( int i = 0; i < len; i++ )
        {
          temp      = DIV( MUL(array1[i], array1[i]), array3[i]);
          array1[i] = DIV( SUB(array2[i], temp), array3[i]-divisor);
          if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
        }
    }

  field1->nmiss = farsetnmiss(len, array1, missval1);
}


void farstd(field_type *field1, field_type field2, field_type field3, int divisor)
{
  int nwpv  = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double *restrict array1  = field1->ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farvar(field1, field2, field3, divisor);

  int nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
  field1->nmiss = nmiss;
}


void farcvar(field_type *field1, field_type field2, int nsets, int divisor)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  const int nsetx = nsets - divisor;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1  = field1->ptr;
  const double *restrict array2  = field2.ptr;
  double temp;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nsetx == 0 )
    {
      for ( int i = 0; i < len; i++ ) array1[i] = missval1;
    }
  else if ( (nmiss1 || nmiss2) /*&& (DBL_IS_NAN(missval1) || DBL_IS_NAN(missval2))*/ )
    {
      for ( int i = 0; i < len; i++ )
        {
          temp      = MULMN(array1[i], array1[i]) / nsets;
          array1[i] = SUBMN(array2[i], temp) / nsetx;
          if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
        }
    }
  else
    {
      for ( int i = 0; i < len; i++ )
        {
          temp      = MUL(array1[i], array1[i]) / nsets;
          array1[i] = SUB(array2[i], temp) / nsetx;
          if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
        }
    }

  field1->nmiss = farsetnmiss(len, array1, missval1);
}


void farcstd(field_type *field1, field_type field2, int nsets, int divisor)
{
  int nwpv  = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double *restrict array1 = field1->ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farcvar(field1, field2, nsets, divisor);

  int nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
  field1->nmiss = nmiss;
}


void farmoq(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = array2[i]*array2[i];
	else
	  array1[i] = missval1;

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] = array2[i]*array2[i];
    }
}


void farmoqw(field_type *field1, field_type field2, double w)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = w*array2[i]*array2[i];
	else
	  array1[i] = missval1;

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] = w*array2[i]*array2[i];
    }
}

/**
 * Counts the number of nonmissing values. The result of the operation
 * is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    miss
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */  
void farcount(field_type *field1, field_type field2)
{
  int nwpv   = field1->nwpv;
  int grid1  = field1->grid;
  int grid2  = field2.grid;
  int nmiss1 = field1->nmiss;
  int nmiss2 = field2.nmiss;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  int len = nwpv*gridInqSize(grid1);

  if ( len != (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( int i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += 1.0;
	    else
	      array1[i] = 1.0;
	  }

      field1->nmiss = 0;
      for ( int i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( int i = 0; i < len; i++ ) 
	array1[i] += 1.0;
    }
}
