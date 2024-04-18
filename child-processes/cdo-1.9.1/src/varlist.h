#ifndef VARLIST_H
#define VARLIST_H

typedef struct {
  bool        check_datarange;
  int         gridsize;
  int         datatype;
  double      missval;
  double      addoffset;
  double      scalefactor;
} varlist_t;

#endif
