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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include "cdo_int.h"
#include "dmemory.h"
#include "listarray.h"


#define  DEFAULT_ALLINC  1024


static void listaInit(lista_t *lista, int type)
{
  lista->array  = NULL;
  lista->nalloc = 0;
  lista->allinc = DEFAULT_ALLINC;
  lista->type   = type;
}


lista_t *lista_new(int type)
{
  lista_t *lista = NULL;

  if ( type != INT_LISTA && type != FLT_LISTA )
    {
      fprintf(stderr, "%s: type %d unsupported!\n", __func__, type);
    }
  else
    {
      lista = (lista_t*) Malloc(sizeof(lista_t));
      listaInit(lista, type);
    }

  return lista;
}


void lista_destroy(lista_t *lista)
{
  if ( lista )
    {
      if ( lista->array ) Free(lista->array);
      Free(lista);
    }
}


void *lista_dataptr(lista_t *lista)
{
  return lista->array;
}


static void listaCheck(lista_t *lista, int num)
{
  while ( lista->nalloc < (num+1) )
    {
      lista->nalloc += lista->allinc;
      if ( lista->type == INT_LISTA )
	lista->array = (int*) Realloc(lista->array, lista->nalloc*sizeof(int));
      else
	lista->array = (double*) Realloc(lista->array, lista->nalloc*sizeof(double));
    }
}


void lista_set_int(lista_t *lista, int num, int ival)
{
  listaCheck(lista, num);

  ((int *) lista->array)[num] = ival;
}


void lista_set_flt(lista_t *lista, int num, double fval)
{
  listaCheck(lista, num);

  ((double *) lista->array)[num] = fval;
}


int lista_get_int(lista_t *lista, int num)
{
  int ival = ((int *) lista->array)[num];

  return ival;
}


double lista_get_flt(lista_t *lista, int num)
{
  double fval = ((double *) lista->array)[num];

  return fval;
}

static
int get_ival(const char *intstr, int idefault, int istart, int iend, int *ilast)
{
  int i;
  int ival = idefault;

  for ( i = istart; i < iend; i++ )
    {
      if ( ! (isdigit(intstr[i]) || intstr[i] == '-') )
	{
	  if ( intstr[i] == '/' )
	    ival = parameter2intlist(intstr+i+1);
	  else
	    fprintf(stderr, "Syntax error in >%.*s<! Character %c not allowed.\n", iend, intstr, intstr[i]);
	  break;
	}
    }

  *ilast = i;

  return ival;
}


void split_intstring(const char *intstr, int *first, int *last, int *inc)
{
  int istrlen = strlen(intstr);
  *first = parameter2intlist(intstr);
  *last  = *first;
  *inc   = 1;

  int i, start = 1;
  *last = get_ival(intstr, *first, start, istrlen, &i);

  if ( i < istrlen )
    {
      start = i+1;
      *inc = get_ival(intstr, 1, start, istrlen, &i);
    }
}


int args2int_lista(int argc, char **argv, lista_t *lista)
{
  int nint = 0;
  int ival;
  int first, last, inc;
  int iarg;

  for ( iarg = 0; iarg < argc; iarg++ )
    {
      split_intstring(argv[iarg], &first, &last, &inc);

      if ( inc >= 0 )
	{
	  for ( ival = first; ival <= last; ival += inc )
	    lista_set_int(lista, nint++, ival);
	}
      else
	{
	  for ( ival = first; ival >= last; ival += inc )
	    lista_set_int(lista, nint++, ival);
	}
    }

  return nint;
}


int args2flt_lista(int argc, char **argv, lista_t *lista)
{
  int i, nint = 0;
  int ival;
  int first, last, inc;
  int iarg;
  int len;
  double tmp_val;

  for ( iarg = 0; iarg < argc; iarg++ )
    {
      len = (int) strlen(argv[iarg]);
      for ( i = 0; i < len; i++ )
	if ( argv[iarg][i] != '/' && argv[iarg][i] != '-' && ! isdigit(argv[iarg][i]) ) break;
      
      if ( i != len )
	{
	  /*
	  if      ( strcmp(argv[iarg],  "inf") == 0 )
	    tmp_val =  DBL_MAX;
	  else if ( strcmp(argv[iarg], "-inf") == 0 )
	    tmp_val = -DBL_MAX;
	  else  
	  */                                    
	    tmp_val = parameter2double(argv[iarg]);

	  lista_set_flt(lista, nint++, tmp_val);
	}
      else
	{
	  split_intstring(argv[iarg], &first, &last, &inc);

	  if ( inc >= 0 )
	    {
	      for ( ival = first; ival <= last; ival += inc )
		lista_set_flt(lista, nint++, (double) ival);
	    }
	  else
	    {
	      for ( ival = first; ival >= last; ival += inc )
		lista_set_flt(lista, nint++, (double) ival);
	    }
	}
    }

  return nint;
}
