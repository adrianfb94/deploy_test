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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "field.h"
#include "expr.h"
#include "expr_fun.h"
#include "expr_yacc.h"


static const char *ExIn[] = {"expr", "init"};
static const char *tmpvnm = "_tmp_";
int pointID = -1;
int surfaceID = -1;

enum {FT_STD, FT_CONST, FT_FLD, FT_VERT, FT_COORD, FT_1C};

#define    COMPLT(x,y)  ((x) < (y))
#define    COMPGT(x,y)  ((x) > (y))
#define    COMPLE(x,y)  ((x) <= (y))
#define    COMPGE(x,y)  ((x) >= (y))
#define    COMPNE(x,y)  IS_NOT_EQUAL(x,y)
#define    COMPEQ(x,y)  IS_EQUAL(x,y)
#define   COMPLEG(x,y)  ((x) < (y) ? -1. : ((x) > (y)))
#define   COMPAND(x,y)  (IS_NOT_EQUAL(x,0) && IS_NOT_EQUAL(y,0))
#define    COMPOR(x,y)  (IS_NOT_EQUAL(x,0) || IS_NOT_EQUAL(y,0))
#define  MVCOMPLT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLT(x,y))
#define  MVCOMPGT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGT(x,y))
#define  MVCOMPLE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLE(x,y))
#define  MVCOMPGE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGE(x,y))
#define  MVCOMPNE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPNE(x,y))
#define  MVCOMPEQ(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPEQ(x,y))
#define MVCOMPLEG(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLEG(x,y))
#define MVCOMPAND(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPAND(x,y))
#define  MVCOMPOR(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPOR(x,y))

static double f_int(double x)          { return (int)(x); }
static double f_nint(double x)         { return round(x); }
static double f_sqr(double x)          { return x*x;      }
static double f_rad(double x)          { return x*M_PI/180.; }
static double f_deg(double x)          { return x*180./M_PI; }
static double pt_ngp(paramType *p)     { return p->ngp;   }
static double pt_nlev(paramType *p)    { return p->nlev;  }
static double pt_size(paramType *p)    { return p->ngp*p->nlev; }
static double pt_missval(paramType *p) { return p->missval; }

typedef struct {
  int type;
  int flag;
  const char *name;  // function name
  double (*func)();    // pointer to function
}
func_t;

static func_t fun_sym_tbl[] =
{
  // scalar functions
  {FT_STD, 0, "abs",   (double (*)()) (double (*)(double)) fabs},
  {FT_STD, 0, "floor", (double (*)()) (double (*)(double)) floor},
  {FT_STD, 0, "ceil",  (double (*)()) (double (*)(double)) ceil},
  {FT_STD, 0, "sqrt",  (double (*)()) (double (*)(double)) sqrt},
  {FT_STD, 0, "exp",   (double (*)()) (double (*)(double)) exp},
  {FT_STD, 0, "erf",   (double (*)()) (double (*)(double)) erf},
  {FT_STD, 0, "log",   (double (*)()) (double (*)(double)) log},
  {FT_STD, 0, "ln",    (double (*)()) (double (*)(double)) log},
  {FT_STD, 0, "log10", (double (*)()) (double (*)(double)) log10},
  {FT_STD, 0, "sin",   (double (*)()) (double (*)(double)) sin},
  {FT_STD, 0, "cos",   (double (*)()) (double (*)(double)) cos},
  {FT_STD, 0, "tan",   (double (*)()) (double (*)(double)) tan},
  {FT_STD, 0, "sinh",  (double (*)()) (double (*)(double)) sinh},
  {FT_STD, 0, "cosh",  (double (*)()) (double (*)(double)) cosh},
  {FT_STD, 0, "tanh",  (double (*)()) (double (*)(double)) tanh},
  {FT_STD, 0, "asin",  (double (*)()) (double (*)(double)) asin},
  {FT_STD, 0, "acos",  (double (*)()) (double (*)(double)) acos},
  {FT_STD, 0, "atan",  (double (*)()) (double (*)(double)) atan},
  {FT_STD, 0, "asinh", (double (*)()) (double (*)(double)) asinh},
  {FT_STD, 0, "acosh", (double (*)()) (double (*)(double)) acosh},
  {FT_STD, 0, "atanh", (double (*)()) (double (*)(double)) atanh},
  {FT_STD, 0, "gamma", (double (*)()) (double (*)(double)) tgamma},
  {FT_STD, 0, "int",   (double (*)()) f_int},
  {FT_STD, 0, "nint",  (double (*)()) f_nint},
  {FT_STD, 0, "sqr",   (double (*)()) f_sqr},
  {FT_STD, 0, "rad",   (double (*)()) f_rad},
  {FT_STD, 0, "deg",   (double (*)()) f_deg},

  // constant functions
  {FT_CONST, 0, "ngp",     (double (*)()) pt_ngp},      // number of horizontal grid points
  {FT_CONST, 0, "nlev",    (double (*)()) pt_nlev},     // number of vertical levels
  {FT_CONST, 0, "size",    (double (*)()) pt_size},     // ngp*nlev
  {FT_CONST, 0, "missval", (double (*)()) pt_missval},  // Returns the missing value of a variable

  // cdo field functions (Reduce grid to point)
  {FT_FLD, 0, "fldmin",  (double (*)()) fldmin},
  {FT_FLD, 0, "fldmax",  (double (*)()) fldmax},
  {FT_FLD, 0, "fldsum",  (double (*)()) fldsum},
  {FT_FLD, 1, "fldmean", (double (*)()) fldmeanw},
  {FT_FLD, 1, "fldavg",  (double (*)()) fldavgw},
  {FT_FLD, 1, "fldstd",  (double (*)()) fldstdw},
  {FT_FLD, 1, "fldstd1", (double (*)()) fldstd1w},
  {FT_FLD, 1, "fldvar",  (double (*)()) fldvarw},
  {FT_FLD, 1, "fldvar1", (double (*)()) fldvar1w},

  // cdo field functions (Reduce level to point)
  {FT_VERT, 0, "vertmin",  (double (*)()) fldmin},
  {FT_VERT, 0, "vertmax",  (double (*)()) fldmax},
  {FT_VERT, 0, "vertsum",  (double (*)()) fldsum},
  {FT_VERT, 1, "vertmean", (double (*)()) fldmeanw},
  {FT_VERT, 1, "vertavg",  (double (*)()) fldavgw},
  {FT_VERT, 1, "vertstd",  (double (*)()) fldstdw},
  {FT_VERT, 1, "vertstd1", (double (*)()) fldstd1w},
  {FT_VERT, 1, "vertvar",  (double (*)()) fldvarw},
  {FT_VERT, 1, "vertvar1", (double (*)()) fldvar1w},
  
  {FT_COORD, 0, "clon",       NULL},
  {FT_COORD, 0, "clat",       NULL},
  {FT_COORD, 0, "clev",       NULL},
  {FT_COORD, 0, "gridarea",   NULL},
  {FT_COORD, 0, "gridweight", NULL},

  {FT_1C, 0, "sellevel",  NULL},
  {FT_1C, 0, "sellevidx", NULL},
};

static int NumFunc = sizeof(fun_sym_tbl) / sizeof(fun_sym_tbl[0]);

static
void node_data_delete(nodeType *p)
{
  if ( p )
    {
      if ( p->param.data ) { Free(p->param.data); p->param.data = NULL;}
    }
}

static
void node_delete(nodeType *p)
{
  if ( p )
    {
      if ( p->type == typeVar ) node_data_delete(p);
      Free(p);
    }
}

static
int get_funcID(const char *fun)
{
  int funcID = -1;
  for ( int i = 0; i < NumFunc; i++ )
    if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
      { 
	funcID = i;
	break;
      }

  if ( funcID == -1 ) cdoAbort("Function >%s< not available!", fun);

  if ( strcmp(fun_sym_tbl[funcID].name, "nint") == 0 ) cdo_check_round();

  return funcID;
}

static
void param_meta_copy(paramType *out, paramType *in)
{
  out->gridID   = in->gridID;
  out->zaxisID  = in->zaxisID;
  out->steptype = in->steptype;
  out->ngp      = in->ngp;
  out->nlev     = in->nlev;
  out->missval  = in->missval;
  out->nmiss    = 0;
  out->coord    = 0;
  out->lmiss    = true;
  out->name     = NULL;
  out->longname = NULL;
  out->units    = NULL;
  out->data     = NULL;
}

static
nodeType *expr_con_con(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type    = typeCon;
  p->ltmpobj = true;

  double cval1 = p1->u.con.value;
  double cval2 = p2->u.con.value;

  switch ( oper )
    {
    case '+':  cval1 = cval1 + cval2; break;
    case '-':  cval1 = cval1 - cval2; break;
    case '*':  cval1 = cval1 * cval2; break;
    case '/':  cval1 = cval1 / cval2; break;
    case '^':  cval1 = pow(cval1, cval2); break;
    default:   cdoAbort("%s: operator %c unsupported!", __func__, oper); break;
    }

  p->u.con.value = cval1;

  return p;
}

static
void oper_expr_con_var(int oper, bool nmiss, size_t n, double missval1, double missval2,
                       double *restrict odat, double cval, const double *restrict idat)
{
  size_t i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADDMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval + idat[i];
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUBMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval - idat[i];
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MULMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval * idat[i];
      break;
    case '/':
      for ( i=0; i<n; ++i ) odat[i] = DIVMN(cval, idat[i]);
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POWMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(cval, idat[i]);
      break;
    case LT:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(cval, idat[i]);
      break;
    case GT:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(cval, idat[i]);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(cval, idat[i]);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(cval, idat[i]);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(cval, idat[i]);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(cval, idat[i]);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(cval, idat[i]);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(cval, idat[i]);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(cval, idat[i]);
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
      break;
    }
}

static
void oper_expr_var_con(int oper, bool nmiss, size_t n, double missval1, double missval2,
                       double *restrict odat, const double *restrict idat, double cval)
{
  size_t i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADDMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] + cval;
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUBMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] - cval;
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MULMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] * cval;
      break;
    case '/':
      if ( nmiss || IS_EQUAL(cval, 0) ) for ( i=0; i<n; ++i ) odat[i] = DIVMN(idat[i], cval);
      else                              for ( i=0; i<n; ++i ) odat[i] = idat[i] / cval;
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POWMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(idat[i], cval);
      break;
    case LT:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(idat[i], cval);
      break;
    case GT:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(idat[i], cval);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(idat[i], cval);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(idat[i], cval);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(idat[i], cval);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(idat[i], cval);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(idat[i], cval);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(idat[i], cval);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(idat[i], cval);
      break;
    default:
      cdoAbort("%s: operator '%c' unsupported!", __func__, oper);
      break;
    }
}

static
void oper_expr_var_var(int oper, bool nmiss, size_t ngp, double missval1, double missval2,
                       double *restrict odat, const double *restrict idat1, const double *restrict idat2)
{
  size_t i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = ADDMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] + idat2[i];
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = SUBMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] - idat2[i];
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MULMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] * idat2[i];
      break;
    case '/':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = DIVMN(idat1[i], idat2[i]);
      else
        {
          for ( i = 0; i < ngp; ++i )
            {
              if ( IS_EQUAL(idat2[i], 0.) ) odat[i] = missval1;
              else                          odat[i] = idat1[i] / idat2[i];
            }
        }
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = POWMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = pow(idat1[i], idat2[i]);
      break;
    case LT:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLT(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLT(idat1[i], idat2[i]);
      break;
    case GT:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGT(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGT(idat1[i], idat2[i]);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLE(idat1[i], idat2[i]);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGE(idat1[i], idat2[i]);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPNE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPNE(idat1[i], idat2[i]);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPEQ(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPEQ(idat1[i], idat2[i]);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLEG(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLEG(idat1[i], idat2[i]);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPAND(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPAND(idat1[i], idat2[i]);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPOR(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPOR(idat1[i], idat2[i]);
      break;
    default:
      cdoAbort("%s: operator %d (%c) unsupported!", __func__, oper, oper);
      break;
    }
}

static
nodeType *expr_con_var(int init, int oper, nodeType *p1, nodeType *p2)
{
  size_t ngp   = p2->param.ngp;
  size_t nlev  = p2->param.nlev;
  size_t nmiss = p2->param.nmiss;
  double missval1 = p2->param.missval;
  double missval2 = p2->param.missval;

  size_t n = ngp*nlev;

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);
  param_meta_copy(&p->param, &p2->param);
  p->param.name = p->u.var.nm;

  if ( ! init )
    {
      p->param.data = (double*) Malloc(n*sizeof(double));
      double *restrict odat = p->param.data;
      const double *restrict idat = p2->param.data;
      double cval = p1->u.con.value;

      oper_expr_con_var(oper, nmiss>0, n, missval1, missval2, odat, cval, idat);

      nmiss = 0;
      for ( size_t i = 0; i < n; i++ )
        if ( DBL_IS_EQUAL(odat[i], missval1) ) nmiss++;

      p->param.nmiss = nmiss;
    }
  
  return p;
}

static
nodeType *expr_var_con(int init, int oper, nodeType *p1, nodeType *p2)
{
  size_t ngp   = p1->param.ngp;
  size_t nlev  = p1->param.nlev;
  size_t nmiss = p1->param.nmiss;
  double missval1 = p1->param.missval;
  double missval2 = p1->param.missval;

  size_t n = ngp*nlev;

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);
  param_meta_copy(&p->param, &p1->param);
  p->param.name = p->u.var.nm;

  if ( ! init )
    {
      p->param.data = (double*) Malloc(n*sizeof(double));
      double *restrict odat = p->param.data;
      const double *restrict idat = p1->param.data;
      double cval = p2->u.con.value;

      oper_expr_var_con(oper, nmiss>0, n, missval1, missval2, odat, idat, cval);

      nmiss = 0;
      for ( size_t i = 0; i < n; i++ )
        if ( DBL_IS_EQUAL(odat[i], missval1) ) nmiss++;

      p->param.nmiss = nmiss;
    }
  
  return p;
}

static
nodeType *expr_var_var(int init, int oper, nodeType *p1, nodeType *p2)
{
  nodeType *px = p1;
  size_t nmiss1 = p1->param.nmiss;
  size_t nmiss2 = p2->param.nmiss;
  double missval1 = p1->param.missval;
  double missval2 = p2->param.missval;

  size_t ngp1 = p1->param.ngp;
  size_t ngp2 = p2->param.ngp;

  size_t ngp = ngp1;

  if ( ngp1 != ngp2 )
    {
      if ( ngp1 == 1 || ngp2 == 1 )
        {
          if ( ngp1 == 1 ) { ngp = ngp2; px = p2; }
        }
      else 
        {
          cdoAbort("%s: Number of grid points differ (%s[%ld] <-> %s[%ld])",
                   __func__, p1->param.name, ngp1, p2->param.name, ngp2);
        }
    }

  size_t nlev1 = p1->param.nlev;
  size_t nlev2 = p2->param.nlev;

  size_t nlev = nlev1;
  if ( nlev1 != nlev2 )
    {
      if ( nlev1 == 1 || nlev2 == 1 )
        {
          if ( nlev1 == 1 ) { nlev = nlev2; px = p2; }
        }
      else 
        {
          cdoAbort("%s: Number of levels differ (%s[%ld] <-> %s[%ld])",
                   __func__, p1->param.name, nlev1, p2->param.name, nlev2);
        }
    }

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);

  param_meta_copy(&p->param, &px->param);

  int steptype1 = p1->param.steptype;
  int steptype2 = p2->param.steptype;

  if ( p->param.steptype == TIME_CONSTANT )
    {
      if ( steptype1 != TIME_CONSTANT )
        p->param.steptype = steptype1;
      else if ( steptype2 != TIME_CONSTANT )
        p->param.steptype = steptype2;
    }

  p->param.name = p->u.var.nm;
  //printf("%s %s nmiss %ld %ld\n", p->u.var.nm, px->param.name, nmiss1, nmiss2);

  if ( ! init )
    {
      p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));

      for ( size_t k = 0; k < nlev; k++ )
        {
          size_t loff1 = 0, loff2 = 0;
          size_t loff = k*ngp;

          if ( nlev1 > 1 ) loff1 = k*ngp1;
          if ( nlev2 > 1 ) loff2 = k*ngp2;

          const double *restrict idat1 = p1->param.data+loff1;
          const double *restrict idat2 = p2->param.data+loff2;
          double *restrict odat = p->param.data+loff;
          int nmiss = nmiss1 > 0 || nmiss2 > 0;

          if ( ngp1 != ngp2 )
            {
              if ( ngp2 == 1 )
                oper_expr_var_con(oper, nmiss, ngp, missval1, missval2, odat, idat1, idat2[0]);
              else
                oper_expr_con_var(oper, nmiss, ngp, missval1, missval2, odat, idat1[0], idat2);
            }
          else
            {
              oper_expr_var_var(oper, nmiss, ngp, missval1, missval2, odat, idat1, idat2);
            }
          }

      size_t nmiss = 0;
      for ( size_t i = 0; i < ngp*nlev; i++ )
        if ( DBL_IS_EQUAL(p->param.data[i], missval1) ) nmiss++;

      p->param.nmiss = nmiss;
    }
  
  return p;
}

static
void ex_copy_var(int init, nodeType *p2, nodeType *p1)
{
  if ( cdoVerbose ) cdoPrint("\t%s\tcopy\t%s[L%lu][N%lu] = %s[L%lu][N%lu]",
                             ExIn[init], p2->param.name, p2->param.nlev, p2->param.ngp, p1->param.name, p2->param.nlev, p2->param.ngp);
  
  size_t ngp = p1->param.ngp;
  assert(ngp > 0);

  if ( ngp != p2->param.ngp )
    cdoAbort("%s: Number of grid points differ (%s[%d] = %s[%d])",
             __func__, p2->param.name, p2->param.ngp, p1->param.name, ngp);

  size_t nlev = p1->param.nlev;
  assert(nlev > 0);

  if ( nlev != p2->param.nlev )
    cdoAbort("%s: Number of levels differ (%s[%d] = %s[%d])",
             __func__, p2->param.name, p2->param.nlev, p1->param.name, nlev);

  if ( ! init )
    {
      double *restrict odat = p2->param.data;
      const double *restrict idat = p1->param.data;
      for ( size_t i = 0; i < ngp*nlev; ++i ) odat[i] = idat[i];
  
      p2->param.missval = p1->param.missval;
      p2->param.nmiss   = p1->param.nmiss;
    }
}

static
void ex_copy_con(int init, nodeType *p2, nodeType *p1)
{
  double cval = p1->u.con.value;

  if ( cdoVerbose ) cdoPrint("\t%s\tcopy\t%s[L%lu][N%lu] = %g", ExIn[init], p2->param.name, p2->param.nlev, p2->param.ngp, cval);
  
  size_t ngp = p2->param.ngp;
  assert(ngp > 0);

  size_t nlev = p2->param.nlev;
  assert(nlev > 0);

  if ( ! init )
    {
      double *restrict odat = p2->param.data;
      assert(odat != NULL);

      for ( size_t i = 0; i < ngp*nlev; ++i ) odat[i] = cval;
    }
}

static
void ex_copy(int init, nodeType *p2, nodeType *p1)
{
  if ( p1->type == typeCon )
    ex_copy_con(init, p2, p1);
  else
    ex_copy_var(init, p2, p1);
}

static
nodeType *expr(int init, int oper, nodeType *p1, nodeType *p2)
{
  if ( p1 == NULL || p2 == NULL ) return  NULL;
  
  const char *coper = "???";

  if ( cdoVerbose )
    {
      switch ( oper )
        {
        case '+':  coper = "+"; break;
        case '-':  coper = "-"; break;
        case '*':  coper = "*"; break;
        case '/':  coper = "/"; break;
        case '^':  coper = "^"; break;
        case LT:   coper = "<"; break;
        case GT:   coper = ">"; break;
        case LE:   coper = "<="; break;
        case GE:   coper = ">="; break;
        case NE:   coper = "!="; break;
        case EQ:   coper = "=="; break;
        case LEG:  coper = "<=>"; break;
        case AND:  coper = "&&"; break;
        case OR:   coper = "||"; break;
        }
    }

 nodeType *p = NULL;

  if ( p1->type == typeVar && p2->type == typeVar )
    {
      p = expr_var_var(init, oper, p1, p2);
      if ( cdoVerbose )
	cdoPrint("\t%s\tarith\t%s[L%lu][N%lu] = %s %s %s", ExIn[init], p->u.var.nm, p->param.nlev, p->param.ngp, p1->u.var.nm, coper, p2->u.var.nm);
    }
  else if ( p1->type == typeVar && p2->type == typeCon )
    {
      p = expr_var_con(init, oper, p1, p2);
      if ( cdoVerbose )
	cdoPrint("\t%s\tarith\t%s[L%lu][N%lu] = %s %s %g", ExIn[init], p->u.var.nm, p->param.nlev, p->param.ngp, p1->u.var.nm, coper, p2->u.con.value);
    }
  else if ( p1->type == typeCon && p2->type == typeVar )
    {
      p = expr_con_var(init, oper, p1, p2);
      if ( cdoVerbose )
	cdoPrint("\t%s\tarith\t%s[L%lu][N%lu] = %g %s %s", ExIn[init], p->u.var.nm, p->param.nlev, p->param.ngp, p1->u.con.value, coper, p2->u.var.nm);
    }
  else if ( p1->type == typeCon && p2->type == typeCon )
    {
      p = expr_con_con(oper, p1, p2);
      if ( cdoVerbose )
	cdoPrint("\t%s\tarith\t%g = %g %s %g", ExIn[init], p->u.con.value, p1->u.con.value, coper, p2->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  if ( p1->ltmpobj ) node_delete(p1);
  if ( p2->ltmpobj ) node_delete(p2);

  return p;
}

static
nodeType *ex_fun_con(int funcID, nodeType *p1)
{
  int functype = fun_sym_tbl[funcID].type;
  if ( functype != FT_STD ) cdoAbort("Function %s not available for constant values!", fun_sym_tbl[funcID].name);

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type    = typeCon;
  p->ltmpobj = true;

  double (*exprfunc)(double) = (double (*)(double)) fun_sym_tbl[funcID].func;
  p->u.con.value = exprfunc(p1->u.con.value);

  if ( p1->ltmpobj ) node_delete(p1);
  else Free(p1);

  return p;
}

static
nodeType *ex_fun_var(int init, int funcID, nodeType *p1)
{
  const char *funcname = fun_sym_tbl[funcID].name;
  int functype = fun_sym_tbl[funcID].type;
  int funcflag = fun_sym_tbl[funcID].flag;

  size_t ngp  = p1->param.ngp;
  size_t nlev = p1->param.nlev;
  size_t nmiss = p1->param.nmiss;
  double missval = p1->param.missval;

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);

  param_meta_copy(&p->param, &p1->param);

  if ( functype == FT_CONST )
    {
      p->type = typeCon;
      double (*exprfunc)(paramType*) = (double (*)(paramType*)) fun_sym_tbl[funcID].func;
      p->u.con.value = exprfunc(&p1->param);
    }
  else if ( functype == FT_FLD )
    {
      p->param.gridID = pointID;
      p->param.ngp    = 1;
    }
  else if ( functype == FT_VERT )
    {
      p->param.zaxisID = surfaceID;
      p->param.nlev    = 1;
    }

  if ( ! init )
    {
      p->param.data = (double*) Malloc(p->param.ngp*p->param.nlev*sizeof(double));
      double *restrict pdata  = p->param.data;
      double *restrict p1data = p1->param.data;
  
      if ( functype == FT_STD )
        {
          double (*exprfunc)(double) = (double (*)(double)) fun_sym_tbl[funcID].func;
          if ( nmiss > 0 )
            {
              for ( size_t i = 0; i < ngp*nlev; i++ )
                {
                  errno = -1;
                  pdata[i] = DBL_IS_EQUAL(p1data[i], missval) ? missval : exprfunc(p1data[i]);
                  if ( errno == EDOM || errno == ERANGE ) pdata[i] = missval;
                  else if ( isnan(pdata[i]) ) pdata[i] = missval;
                }
            }
          else
            {
              for ( size_t i = 0; i < ngp*nlev; i++ )
                {
                  errno = -1;
                  pdata[i] = exprfunc(p1data[i]);
                  if ( errno == EDOM || errno == ERANGE ) pdata[i] = missval;
                  else if ( isnan(pdata[i]) ) pdata[i] = missval;
                }
            }
        }
      else if ( functype == FT_FLD )
        {
          field_type field;
          double *weights = NULL;
          //if ( funcflag == 1 ) weights = fld_weights(p1->param.gridID, ngp);
          if ( funcflag == 1 )
            {
              weights = p1->param.weight;
              assert(weights!=NULL);
            }
          
          double (*exprfunc)(field_type) = (double (*)(field_type)) fun_sym_tbl[funcID].func;
          for ( size_t k = 0; k < nlev; k++ )
            {
              fld_field_init(&field, nmiss, missval, ngp, p1data+k*ngp, weights);
              pdata[k] = exprfunc(field);
            }
          //if ( weights ) Free(weights);
        }
      else if ( functype == FT_VERT )
        {
          field_type field;
          double *weights = NULL;
          if ( funcflag == 1 ) weights = vert_weights(p1->param.zaxisID, nlev);
          double *array = (double*) Malloc(nlev*sizeof(double));
          double (*exprfunc)(field_type) = (double (*)(field_type)) fun_sym_tbl[funcID].func;
          for ( size_t i = 0; i < ngp; i++ )
            {
              for ( size_t k = 0; k < nlev; k++ ) array[k] = p1data[k*ngp+i];
              fld_field_init(&field, nmiss, missval, nlev, array, weights);
              pdata[i] = exprfunc(field);
            }
          if ( array ) Free(array);
          if ( weights ) Free(weights);
        }
      else if ( functype == FT_CONST )
        {
        }
      else
        cdoAbort("Intermal error, wrong function type (%d) for %s()!", functype, funcname);

      nmiss = 0;
      for ( size_t i = 0; i < p->param.ngp*p->param.nlev; i++ )
        if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

      p->param.nmiss = nmiss;
    }

  if ( p1->ltmpobj ) node_delete(p1);
  else Free(p1);
  
  return p;
}

static
nodeType *ex_fun(int init, int funcID, nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      if ( cdoVerbose ) cdoPrint("\t%s\tfunc\t%s (%s)", ExIn[init], fun_sym_tbl[funcID].name, p1->u.var.nm);
      p = ex_fun_var(init, funcID, p1);
    }
  else if ( p1->type == typeCon )
    {
      if ( cdoVerbose ) cdoPrint("\t%s\tfunc\t%s (%g)", ExIn[init], fun_sym_tbl[funcID].name, p1->u.con.value);
      p = ex_fun_con(funcID, p1);
    }
  else
    cdoAbort("Internal problem!");

  return p;
}

static
size_t get_levidx(size_t nlev, const double *data, double value, const char *funcname)
{
  size_t levidx;
  
  for ( levidx = 0; levidx < nlev; ++levidx )
    if ( IS_EQUAL(data[levidx], value) )
      break;
  if ( levidx == nlev ) cdoAbort("%s(): level %g not found!", funcname, value);

  return levidx;
}

static
nodeType *fun1c(int init, int funcID, nodeType *p1, double value, parse_param_t *parse_arg)
{  
  const char *funcname = fun_sym_tbl[funcID].name;            
  if ( p1->type != typeVar ) cdoAbort("Parameter of function %s() needs to be a variable!", funcname);
  if ( p1->ltmpobj ) cdoAbort("Temorary objects not allowed in function %s()!", funcname);

  size_t ngp   = p1->param.ngp;
  size_t nlev  = p1->param.nlev;
  size_t nmiss = p1->param.nmiss;
  double missval = p1->param.missval;

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);
  param_meta_copy(&p->param, &p1->param);
  p->param.name = p->u.var.nm;

  p->param.nlev = 1;

  int zaxisID = p1->param.zaxisID;
  int coordID = params_get_coordID(parse_arg, 'z', zaxisID);
  size_t levidx = 0;
  
  double *data = NULL;
  if ( init )
    {
      parse_arg->coords[coordID].needed = true;

      data = (double*) Malloc(nlev*sizeof(double));
      cdoZaxisInqLevels(zaxisID, data);
    }
  else
    {
      data = parse_arg->coords[coordID].data;
    }

  if ( strcmp(funcname, "sellevidx") == 0 )
    {
      long ilevidx = lround(value);
      if ( ilevidx < 1 || ilevidx > (long)nlev )
        cdoAbort("%s(): level index %ld out of range (range: 1-%lu)!", funcname, ilevidx, nlev);
      levidx = (size_t) ilevidx - 1;
    }
  else if ( strcmp(funcname, "sellevel") == 0 )
    {
      levidx = get_levidx(nlev, data, value, funcname);
    }
  else
    cdoAbort("Function %s() not implemented!", funcname);

  if ( init )
    {
      double level = data[levidx];
      int zaxisID2 = zaxisCreate(zaxisInqType(zaxisID), 1);
      zaxisDefLevels(zaxisID2, &level);
      p->param.zaxisID = zaxisID2;
    }

  if ( init ) Free(data);

  if ( ! init )
    {
      p->param.data = (double*) Malloc(ngp*sizeof(double));
      double *restrict pdata = p->param.data;
      const double *restrict p1data = p1->param.data+ngp*levidx;

      for ( size_t i = 0; i < ngp; i++ ) pdata[i] = p1data[i];

      if ( nmiss > 0 )
        {
          nmiss = 0;
          for ( size_t i = 0; i < ngp; i++ )
            if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;
        }

      p->param.nmiss = nmiss;
    }

  if ( p1->ltmpobj ) node_delete(p1);

  return p;
}

static
nodeType *coord_fun(int init, int funcID, nodeType *p1, parse_param_t *parse_arg)
{  
  const char *funcname = fun_sym_tbl[funcID].name;            
  if ( p1->type != typeVar ) cdoAbort("Parameter of function %s() needs to be a variable!", funcname);
  if ( p1->ltmpobj ) cdoAbort("Temorary objects not allowed in function %s()!", funcname);
            
  size_t len = 3 + strlen(p1->u.var.nm);
  char *cname = (char*) Calloc(len, 1);
  strcpy(cname, p1->u.var.nm);
            
  if      ( strcmp(funcname, "clon") == 0 ) strcat(cname, ".x");
  else if ( strcmp(funcname, "clat") == 0 ) strcat(cname, ".y");
  else if ( strcmp(funcname, "clev") == 0 ) strcat(cname, ".z");
  else if ( strcmp(funcname, "gridarea")   == 0 ) strcat(cname, ".a");
  else if ( strcmp(funcname, "gridweight") == 0 ) strcat(cname, ".w");
  else cdoAbort("Implementation missing for function %s!", funcname);
  
  Free(p1->u.var.nm);
  p1->u.var.nm = cname;

  nodeType *p = expr_run(p1, parse_arg);
  p->param.lmiss = false;

  if ( ! init )
    {
      /*
      size_t ngp  = p1->param.ngp;
      size_t nlev = p1->param.nlev;
      p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));
      double *restrict pdata  = p->param.data;
      double *restrict p1data = p1->param.data;

      for ( size_t i = 0; i < ngp*nlev; i++ ) pdata[i] = p1data[i];
      */
    }
  /*
  Free(cname);
  Free(p1);
  */
  return p;
}

static
nodeType *ex_uminus_var(int init, nodeType *p1)
{
  size_t ngp   = p1->param.ngp;
  size_t nlev  = p1->param.nlev;
  size_t nmiss = p1->param.nmiss;
  double missval = p1->param.missval;

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);
  param_meta_copy(&p->param, &p1->param);
  p->param.name = p->u.var.nm;

  if ( ! init )
    {
      p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));
      double *restrict pdata = p->param.data;
      const double *restrict p1data = p1->param.data;

      if ( nmiss > 0 )
        {
          for ( size_t i = 0; i < ngp*nlev; ++i )
            pdata[i] = DBL_IS_EQUAL(p1data[i], missval) ? missval : -(p1data[i]);
        }
      else
        {
          for ( size_t i = 0; i < ngp*nlev; ++i )
            pdata[i] = -(p1data[i]);
        }

      p->param.nmiss = nmiss;
    }

  if ( p1->ltmpobj ) node_delete(p1);
  
  return p;
}

static
nodeType *ex_uminus_con(nodeType *p1)
{
  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type    = typeCon;
  p->ltmpobj = true;

  p->u.con.value = -(p1->u.con.value);

  return p;
}

static
nodeType *ex_uminus(int init, nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      if ( cdoVerbose ) cdoPrint("\t%s\tneg\t- (%s)", ExIn[init], p1->u.var.nm);
      p = ex_uminus_var(init, p1);
    }
  else if ( p1->type == typeCon )
    {
      if ( cdoVerbose ) cdoPrint("\t%s\tneg\t- (%g)", ExIn[init], p1->u.con.value);
      p = ex_uminus_con(p1);
    }
  else
    cdoAbort("Internal problem!");

  return p;
}

static
nodeType *ex_ifelse(int init, nodeType *p1, nodeType *p2, nodeType *p3)
{
  if ( p1->type == typeCon ) cdoAbort("expr?expr:expr: First expression is a constant but must be a variable!");

  if ( cdoVerbose )
    {
      fprintf(stderr, "cdo expr:\t%s\tifelse\t%s[L%lu][N%lu] ? ", ExIn[init], p1->u.var.nm, p1->param.nlev, p1->param.ngp);
      if ( p2->type == typeCon )
        fprintf(stderr, "%g : ", p2->u.con.value);
      else
        fprintf(stderr, "%s[L%lu][N%lu] : ", p2->u.var.nm, p2->param.nlev, p2->param.ngp);
      if ( p3->type == typeCon )
        fprintf(stderr, "%g\n", p3->u.con.value);
      else
        fprintf(stderr, "%s[L%lu][N%lu]\n", p3->u.var.nm, p3->param.nlev, p3->param.ngp);
    }

  size_t nmiss1 = p1->param.nmiss;
  size_t ngp1   = p1->param.ngp;
  size_t nlev1  = p1->param.nlev;
  double missval1 = p1->param.missval;

  size_t ngp = ngp1;
  size_t nlev = nlev1;
  nodeType *px = p1;

  double missval2 = missval1;
  double *pdata2 = NULL;
  size_t ngp2 = 1;
  size_t nlev2 = 1;
  
  if ( p2->type == typeCon )
    {
      pdata2 = &p2->u.con.value;
    }
  else
    {
      ngp2 = p2->param.ngp;
      nlev2 = p2->param.nlev;
      missval2 = p2->param.missval;
      pdata2 = p2->param.data;

      if ( ngp2 > 1 && ngp2 != ngp )
        {
          if ( ngp == 1 )
            {
              ngp = ngp2;
              px = p2;
            }
          else
            cdoAbort("expr?expr:expr: Number of grid points differ (ngp1 = %ld, ngp2 = %ld)", ngp1, ngp2);
        }

      if ( nlev2 > 1 && nlev2 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev2;
	      px = p2;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ (nlev = %ld, nlev2 = %ld)", nlev, nlev2);
	}
    }

  double missval3 = missval1;
  double *pdata3 = NULL;
  size_t ngp3 = 1;
  size_t nlev3 = 1;
  
  if ( p3->type == typeCon )
    {
      pdata3 = &p3->u.con.value;
    }
  else
    {
      ngp3 = p3->param.ngp;
      nlev3 = p3->param.nlev;
      missval3 = p3->param.missval;
      pdata3 = p3->param.data;

      if ( ngp3 > 1 && ngp3 != ngp )
        {
          if ( ngp == 1 )
            {
              ngp = ngp3;
              px = p3;
            }
          else
            cdoAbort("expr?expr:expr: Number of grid points differ (ngp1 = %ld, ngp3 = %ld)", ngp1, ngp3);
        }

      if ( nlev3 > 1 && nlev3 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev3;
	      px = p3;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ (nlev = %ld, nlev3 = %ld)", nlev, nlev3);
	}
    }

  nodeType *p = (nodeType*) Calloc(1, sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpobj  = true;
  p->u.var.nm = strdup(tmpvnm);

  param_meta_copy(&p->param, &px->param);
  p->param.name = p->u.var.nm;

  if ( ! init )
    {
      size_t nmiss = 0;
      double *pdata1 = p1->param.data;
      
      p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));

      for ( size_t k = 0; k < nlev; ++k )
        {
          size_t loff1 = (nlev1 == 1) ? 0 : k*ngp1;
          size_t loff  = k*ngp;
          size_t loff2 = (nlev2 == 1) ? 0 : loff;
          size_t loff3 = (nlev3 == 1) ? 0 : loff;

          const double *restrict idat1 = pdata1+loff1;
          const double *restrict idat2 = pdata2+loff2;
          const double *restrict idat3 = pdata3+loff3;
          double *restrict odat = p->param.data+loff;

          double ival1 = idat1[0];
          double ival2 = idat2[0];
          double ival3 = idat3[0];
          for ( size_t i = 0; i < ngp; ++i ) 
            {
              if ( ngp1 > 1 ) ival1 = idat1[i];
              if ( ngp2 > 1 ) ival2 = idat2[i];
              if ( ngp3 > 1 ) ival3 = idat3[i];

              if ( nmiss1 && DBL_IS_EQUAL(ival1, missval1) )
                odat[i] = missval1;
              else if ( IS_NOT_EQUAL(ival1, 0) )
                odat[i] = DBL_IS_EQUAL(ival2, missval2) ? missval1 : ival2;
              else
                odat[i] = DBL_IS_EQUAL(ival3, missval3) ? missval1 : ival3;
            }

          for ( size_t i = 0; i < ngp; i++ )
            if ( DBL_IS_EQUAL(odat[i], missval1) ) nmiss++;
        }

      p->param.nmiss = nmiss;
    }

  if ( p1->ltmpobj ) node_delete(p1);
  if ( p2->ltmpobj ) node_delete(p2);
  if ( p3->ltmpobj ) node_delete(p3);

  return p;
}
/*
static
int exNode(nodeType *p, parse_param_t *parse_arg)
{
  if ( ! p ) return 0;

  // node is leaf
  if ( p->type == typeCon || p->type == typeVar || p->u.opr.nops == 0 )
    {
      return 0;
    }

  // node has children
  for ( int k = 0; k < p->u.opr.nops; k++ )
    {
      exNode(p->u.opr.op[k], parse_arg);
    }

  return 0;
}
*/

static
int param_search_name(int nparam, paramType *params, const char *name)
{
  int varID = -1;

  for ( varID = nparam-1; varID >= 0; --varID )
    {
      if ( strcmp(params[varID].name, name) == 0 ) break;
    }

  return varID;
}


nodeType *expr_run(nodeType *p, parse_param_t *parse_arg)
{
  pointID = parse_arg->pointID;
  surfaceID = parse_arg->surfaceID;
  int init = parse_arg->init;
  paramType *params = parse_arg->params;
  int varID;
  nodeType *rnode = NULL;

  if ( ! p ) return rnode;

  /*  if ( ! init ) { exNode(p, parse_arg); return 0; } */

  switch ( p->type )
    {
    case typeCom:
      {
        const char *cname = p->u.com.cname;
        const char *vname = p->u.com.vname;
        if ( parse_arg->debug ) cdoPrint("\tstatement\t\t%s(%s)", cname, vname);

        varID = param_search_name(parse_arg->nparams, params, vname);
        if ( varID == -1 ) cdoAbort("Variable %s not found, needed for statement %s(%s)!", vname, cname, vname);

        if ( init )
          {
            if ( strcmp(cname, "remove") == 0 )
              {
                params[varID].remove = true;
              }
          }
        else
          {
            if ( strcmp(cname, "print") == 0 )
              {
                size_t maxout = 100;
                int vartsID = parse_arg->tsID;
                int steptype = params[varID].steptype;
                size_t ngp = params[varID].ngp;
                size_t nlev = params[varID].nlev;
                const double *data = params[varID].data;
                long tsID = lround(params[vartsID].data[0]);
                for ( size_t k = 0; k < nlev; ++k )
                  for ( size_t i = 0; i < ngp; ++i )
                    {
                      if ( i < maxout || (ngp > maxout && i >= (ngp-maxout)) )
                        {
                          if ( steptype == TIME_CONSTANT )
                            fprintf(stdout, "   %s[lev=%lu:gp=%lu] = %g\n", vname, k+1, i+1, data[k*ngp+i]);
                          else
                            fprintf(stdout, "   %s[ts=%ld:lev=%lu:gp=%lu] = %g\n", vname, tsID, k+1, i+1, data[k*ngp+i]);
                        }
                      else if ( i == maxout )
                        {
                          fprintf(stdout, "   .......\n");
                        }
                    }
              }
          }
        
        break;
      }
    case typeCon:
      {
        if ( parse_arg->debug ) cdoPrint("\tpush\tconst\t%g", p->u.con.value);

        rnode = p;

        break;
      }
    case typeVar:
      {
        const char *vnm = p->u.var.nm;
        varID = param_search_name(parse_arg->nparams, params, vnm);
        if ( varID == -1 && init )
          {
            size_t len = strlen(vnm);
            int coord = vnm[len-1];
            if ( len > 2 && vnm[len-2] == '.' )
              {
                if ( coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w' )
                  {
                    char *varname = strdup(vnm);
                    varname[len-2] = 0;
                    varID = param_search_name(parse_arg->nparams, params, varname);
                    free(varname);
                    if ( varID == -1 )
                      {
                        cdoAbort("Coordinate %c: variable >%s< not found!", coord, varname);
                      }
                    else
                      {
                        int nvarID = parse_arg->nparams;
                        if ( nvarID >= parse_arg->maxparams )
                          cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);

                        int coordID = params_get_coordID(parse_arg, coord, params[varID].gridID);
                        parse_arg->coords[coordID].needed = true;
                        const char *units = parse_arg->coords[coordID].units;
                        const char *longname = parse_arg->coords[coordID].longname;

                        params[nvarID].coord    = coord;
                        params[nvarID].lmiss    = false;
                        params[nvarID].name     = strdup(vnm);
                        params[nvarID].missval  = params[varID].missval;
                        params[nvarID].gridID   = params[varID].gridID;
                        params[nvarID].zaxisID  = parse_arg->surfaceID;
                        params[nvarID].steptype = TIME_CONSTANT;
                        params[nvarID].ngp      = params[varID].ngp;
                        params[nvarID].nlev     = 1;
                        if ( units ) params[nvarID].units = strdup(units);
                        if ( longname ) params[nvarID].longname = strdup(longname);
                        parse_arg->nparams++;
                        varID = nvarID;
                      }
                  }
                else if ( coord == 'z' )
                  {
                    char *varname = strdup(vnm);
                    varname[len-2] = 0;
                    varID = param_search_name(parse_arg->nparams, params, varname);
                    free(varname);
                    if ( varID == -1 )
                      {
                        cdoAbort("Coordinate %c: variable >%s< not found!", coord, varname);
                      }
                    else
                      {
                        int nvarID = parse_arg->nparams;
                        if ( nvarID >= parse_arg->maxparams )
                          cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);
                          
                        int coordID = params_get_coordID(parse_arg, coord, params[varID].zaxisID);
                        parse_arg->coords[coordID].needed = true;
                        const char *units = parse_arg->coords[coordID].units;
                        const char *longname = parse_arg->coords[coordID].longname;
                                     
                        params[nvarID].coord    = coord;
                        params[nvarID].lmiss    = false;
                        params[nvarID].name     = strdup(vnm);
                        params[nvarID].missval  = params[varID].missval;
                        params[nvarID].gridID   = parse_arg->pointID;
                        params[nvarID].zaxisID  = params[varID].zaxisID;
                        params[nvarID].steptype = TIME_CONSTANT;
                        params[nvarID].ngp      = 1;
                        params[nvarID].nlev     = params[varID].nlev;
                        if ( units ) params[nvarID].units = strdup(units);
                        if ( longname ) params[nvarID].longname = strdup(longname);
                        parse_arg->nparams++;
                        varID = nvarID;
                      }
                  }
              }
          }
        if ( varID == -1 )
          {
            cdoAbort("Variable >%s< not found!", p->u.var.nm);
          }
        else if ( init )
          {
            if ( varID < parse_arg->nvars1 && parse_arg->needed[varID] == false )
              {
                parse_arg->needed[varID] = true;
              }
          }


        param_meta_copy(&p->param, &params[varID]);
        p->param.coord    = params[varID].coord;
        p->param.lmiss    = params[varID].lmiss;
        p->param.name     = params[varID].name;
        p->param.longname = params[varID].longname;
        p->param.units    = params[varID].units;
        p->ltmpobj = false;
        /*
        if ( parse_arg->debug )
          printf("var: u.var.nm=%s name=%s gridID=%d zaxisID=%d ngp=%d nlev=%d  varID=%d\n",
                 p->u.var.nm, p->param.name, p->param.gridID, p->param.zaxisID, p->param.ngp, p->param.nlev, varID);
        */
        if ( ! init )
          {
            p->param.data  = params[varID].data;
            p->param.nmiss = params[varID].nmiss;
          }

        if ( parse_arg->debug ) cdoPrint("\tpush\tvar\t%s[L%lu][N%lu]", vnm, p->param.nlev, p->param.ngp);

        rnode = p;

        break;
      }
    case typeFun1c:
      {
        int funcID = get_funcID(p->u.fun1c.name);
        int functype = fun_sym_tbl[funcID].type;
        
        nodeType *fnode = expr_run(p->u.fun1c.op, parse_arg);
        
        if ( functype == FT_1C )
          {
            double value = p->u.fun1c.value;
            rnode = fun1c(init, funcID, fnode, value, parse_arg);
          }

        break;
      }
    case typeFun:
      {
        int funcID = get_funcID(p->u.fun.name);
        int functype = fun_sym_tbl[funcID].type;

        nodeType *fnode = expr_run(p->u.fun.op, parse_arg);

        if ( functype == FT_COORD )
          {
            rnode = coord_fun(init, funcID, fnode, parse_arg);
          }
        else
          {
            int functype = fun_sym_tbl[funcID].type;
            int funcflag = fun_sym_tbl[funcID].flag;
            if ( functype == FT_FLD && funcflag == 1 )
              {
                int coordID = params_get_coordID(parse_arg, 'w', fnode->param.gridID);
                if ( init )
                  parse_arg->coords[coordID].needed = true;
                else
                  fnode->param.weight = parse_arg->coords[coordID].data;
              }
            rnode = ex_fun(init, funcID, fnode);
            // if ( fnode->ltmpobj ) node_delete(fnode);
            // Free(fnode);
          }
        
        break;
      }
    case typeOpr:
      switch( p->u.opr.oper )
	{
        case '=':
          {
            rnode = expr_run(p->u.opr.op[1], parse_arg);

            const char *varname2 = p->u.opr.op[0]->u.var.nm;

            if ( parse_arg->debug )
              {
                if ( rnode && rnode->type == typeVar)
                  cdoPrint("\tpop\tvar\t%s[L%lu][N%lu]", varname2, rnode->param.nlev, rnode->param.ngp);
                else
                  cdoPrint("\tpop\tconst\t%s", varname2);
              }

            if ( init )
              {
                varID = param_search_name(parse_arg->nparams, params, varname2);
                if ( varID >= 0 )
                  {
                    if ( varID < parse_arg->nvars1 )
                      {
                        params[varID].select = true;
                        parse_arg->needed[varID] = true;
                      }
                    else if ( params[varID].coord )
                      cdoAbort("Coordinate variable %s is read only!", varname2);
                    /*
                      else
                      cdoWarning("Variable %s already defined!", varname2);
                    */
                  }
                else if ( p->u.opr.op[1]->type != typeCon )
                  {
                    varID = parse_arg->nparams;
                    if ( varID >= parse_arg->maxparams )
                      cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);

                    param_meta_copy(&params[varID], &rnode->param);
                    params[varID].coord = 0;
                    params[varID].lmiss = rnode->param.lmiss;
                    params[varID].name  = strdup(varname2);
                    params[varID].nmiss = rnode->param.nmiss;
                    if ( rnode->param.units ) params[varID].units = strdup(rnode->param.units);
                    if ( rnode->param.longname ) params[varID].longname = strdup(rnode->param.longname);
                    parse_arg->nparams++;
                  }
              }
            else
              {
                varID = param_search_name(parse_arg->nparams, params, varname2);
                if ( varID < 0 ) cdoAbort("Variable >%s< not found!", varname2);
                else if ( params[varID].coord ) cdoAbort("Coordinate variable %s is read only!", varname2);
                param_meta_copy(&p->param, &params[varID]);
                p->param.name  = params[varID].name;
                p->param.data  = params[varID].data;
                p->ltmpobj     = false;
                
                ex_copy(init, p, rnode);
                params[varID].nmiss = p->param.nmiss;
              }
            
            if ( rnode && rnode->ltmpobj ) { node_delete(rnode); rnode = NULL; }
            // else Free(rnode);

            break;
          }
        case UMINUS:
          {
            rnode = ex_uminus(init, expr_run(p->u.opr.op[0], parse_arg));

            break;
          }
        case '?':
          {
            rnode = ex_ifelse(init,
                              expr_run(p->u.opr.op[0], parse_arg),
                              expr_run(p->u.opr.op[1], parse_arg),
                              expr_run(p->u.opr.op[2], parse_arg));
            
            break;
          }
        default:
          {
            rnode = expr(init, p->u.opr.oper,
                         expr_run(p->u.opr.op[0], parse_arg),
                         expr_run(p->u.opr.op[1], parse_arg));

            break;
          }
        }
      break;
    }

  return rnode;
}
