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

#include "cdo.h"
#include "cdo_int.h"
#include <cdi.h>


void farcfun(field_type *field, double rconst, int function)
{
  switch (function)
    {
    case func_add: farcadd(field, rconst);  break;
    case func_sub: farcsub(field, rconst);  break;
    case func_mul: farcmul(field, rconst);  break;
    case func_div: farcdiv(field, rconst);  break;
    case func_mod: farmod(field, rconst);   break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
}

void farcmul(field_type *field, double rconst)
{
  int i, len;
  int    nwpv     = field->nwpv;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len    = nwpv*gridInqSize(grid);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array[i] = MULMN(array[i], rconst);
    }
  else
    {
      /*
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i)
#endif
      */
      for ( i = 0; i < len; i++ ) 
	array[i] *= rconst;
    }
}


void farcdiv(field_type *field, double rconst)
{
  int i, len;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  if ( nmiss > 0 || IS_EQUAL(rconst, 0) )
    {
      for ( i = 0; i < len; i++ )
	array[i] = DIVMN(array[i], rconst);

      if ( IS_EQUAL(rconst, 0) ) field->nmiss = len;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array[i] /= rconst;
    }
}


void farcadd(field_type *field, double rconst)
{
  int i, len;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array[i] = ADDMN(array[i], rconst);
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array[i] += rconst;
    }
}


void farcsub(field_type *field, double rconst)
{
  farcadd(field, -rconst);
}


void farinv(field_type *field)
{
  int    grid     = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  int len = gridInqSize(grid);

  for ( int i = 0; i < len; i++ ) 
    array[i] = DIVMN(1.0, array[i]);

  field->nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array[i], missval1) ) field->nmiss++;
}


void farround(field_type *field)
{
  int    grid     = field->grid;
  double missval1 = field->missval;
  double *array   = field->ptr;

  int len = gridInqSize(grid);

  for ( int i = 0; i < len; i++ ) 
    array[i] = round(array[i]);

  field->nmiss = 0;
  for ( int i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array[i], missval1) ) field->nmiss++;
}


void farmod(field_type *field, double divisor)
{
  int    grid     = field->grid;
  double missval1 = field->missval;
  double *array   = field->ptr;

  int len = gridInqSize(grid);

  for ( int i = 0; i < len; i++ )
    {
      array[i] = DBL_IS_EQUAL(array[i], missval1) ? missval1 : fmod(array[i], divisor);
    }
}
