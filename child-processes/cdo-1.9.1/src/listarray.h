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

#ifndef _LISTA_H
#define _LISTA_H

#define  INT_LISTA  1
#define  FLT_LISTA  2


typedef struct {
  void *array;
  int nalloc;
  int allinc;
  int type;
}
lista_t;


lista_t *lista_new(int type);
void lista_destroy(lista_t *lista);
void *lista_dataptr(lista_t *lista);
void lista_set_int(lista_t *lista, int num, int ival);
void lista_set_flt(lista_t *lista, int num, double fval);
int lista_get_int(lista_t *lista, int num);
double lista_get_flt(lista_t *lista, int num);
int args2int_lista(int argc, char **argv, lista_t *lista);
int args2flt_lista(int argc, char **argv, lista_t *lista);

#endif  /* _LISTA_H */
