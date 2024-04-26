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
#include <stdbool.h>

#ifdef __cplusplus
#ifndef fileno
int fileno(FILE *stream);
#endif
#endif

extern int CDO_parser_errorno;


typedef enum { typeCon, typeVar, typeFun, typeFun1c, typeOpr, typeCom } nodeEnum;

// commands
typedef struct {
  char *cname;                // command name
  char *vname;                // variable name
} comNodeType;

// constants
typedef struct {
  double value;               // value of constant
} conNodeType;

// variables
typedef struct {
  char *nm;                   // variable name
} varNodeType;

// 1c functions
typedef struct {
  char *name;                 // function name
  double value;               // value of constant
  struct nodeTypeTag *op;     // operand
} fun1cNodeType;

// functions
typedef struct {
  char *name;                 // function name
  struct nodeTypeTag *op;     // operand
} funNodeType;

// operators
typedef struct {
  int oper;                   // operator             
  int nops;                   // number of operands   
  struct nodeTypeTag *op[1];  // operands (expandable)
} oprNodeType;

// parameter
typedef struct {
  bool       select;
  bool       remove;
  bool       lmiss;
  int        coord;
  int        gridID;
  int        zaxisID;
  int        steptype;
  size_t     ngp;
  size_t     nlev;
  size_t     nmiss;
  char      *name;
  char      *longname;
  char      *units;
  double     missval;
  double    *data;
  double    *weight;
} paramType;


typedef struct nodeTypeTag {
  bool      ltmpobj;
  paramType param;

  nodeEnum  type;             // type of node

  // union must be last entry in nodeType
  // because operNodeType may dynamically increase
  union {
    comNodeType com;          // commands
    conNodeType con;          // constants
    varNodeType var;          // variables
    funNodeType fun;          // functions
    fun1cNodeType fun1c;      // functions
    oprNodeType opr;          // operators
  } u;
} nodeType;


typedef struct {
  bool       needed;
  int        coord;
  int        cdiID;
  int        size;
  char      *units;
  char      *longname;
  double    *data;
} coordType;


typedef struct {
  bool       init;
  bool       debug;
  bool      *needed;
  int        maxparams;
  int        nparams;
  int        nvars1;
  int        ncoords;
  int        maxcoords;
  int        tsID;
  int        pointID;
  int        surfaceID;
  coordType *coords;
  paramType *params;
} parse_param_t;


typedef union{
  double    cvalue;           // constant value
  char     *varnm;            // variable name 
  char     *fname;            // function name 
  nodeType *nPtr;             // node pointer  
} stype_t;


#define YYSTYPE        stype_t
#define YY_EXTRA_TYPE  parse_param_t *

#define YY_DECL int yylex(YYSTYPE *yylval_param, parse_param_t *parse_arg, void *yyscanner)
YY_DECL;

int  yyparse(parse_param_t *parse_arg, void*);
void yyerror(void *parse_arg, void *scanner, const char *errstr);

int  yylex_init(void **);
int  yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);

nodeType *expr_run(nodeType *p, parse_param_t *parse_arg);
int params_get_coordID(parse_param_t *parse_arg, int coord, int cdiID);

