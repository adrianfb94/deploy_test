#include <errno.h>
#include <limits.h>
#include "cdo_int.h"
#include "cdi.h"


int literal_get_datatype(const char *literal)
{
  if ( literal && *literal )
    {
      char *endptr;
      errno = 0;
      long lval = strtol(literal, &endptr, 10);
      if ( errno == 0 && *endptr == 0 ) return CDI_DATATYPE_INT32;
      else if ( errno == 0 && *(endptr+1) == 0 && (*endptr == 's' || *endptr == 'b') )
        {
          if      ( *endptr == 's' && lval >= SHRT_MIN && lval <= SHRT_MAX )
            return CDI_DATATYPE_INT16;
          else if ( *endptr == 'b' && lval >= SCHAR_MIN && lval <= SCHAR_MAX )
            return CDI_DATATYPE_INT8;
        }
      else
        {
          errno = 0;
          double dval = strtod(literal, &endptr);
          if ( errno == 0 && *endptr == 0 ) return CDI_DATATYPE_FLT64;
          else if ( errno == 0 && *(endptr+1) == 0 )
            {
              if ( *endptr == 'f' && dval >= -FLT_MAX && dval <= FLT_MAX )
                return CDI_DATATYPE_FLT32;
            }
        }
    }

  return -1;
}


int literals_find_datatype(int n, char **literals)
{
  int dtype = -1;

  if ( n )
    {
      dtype = literal_get_datatype(literals[0]);
      if ( dtype != -1 )
        for ( int i = 1; i < n; ++i )
          {
            int xtype = literal_get_datatype(literals[i]);
            if ( dtype != xtype )
              {
                if ( xtype == CDI_DATATYPE_FLT32 || xtype == CDI_DATATYPE_FLT64 )
                  {
                    if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
                      {
                        if ( xtype > dtype ) dtype = xtype;
                      }
                    else dtype = xtype;
                  }
                else
                  {
                    if ( !(dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64) )
                      {
                        if ( xtype > dtype ) dtype = xtype;
                      }
                    else dtype = xtype;
                  }
              }
          }
    }

  return dtype;
}


int literal_to_int(const char *literal)
{
  int ival = INT_MAX;

  if ( literal && *literal )
    {
      char *endptr;
      ival = strtol(literal, &endptr, 10);
    }

  return ival;
}


double literal_to_double(const char *literal)
{
  double dval = DBL_MAX;

  if ( literal && *literal )
    {
      char *endptr;
      dval = strtod(literal, &endptr);
    }

  return dval;
}


#ifdef TEST_LITERAL

int main(void)
{
  const char *literals[] = {"127b", "-32768s", "-2147483647", "-1.e+36f", "1.e+308", "temperature", "surface pressure", "1000."};
  int nliterals = sizeof(literals) / sizeof(literals[0]);

  for ( int i = 0; i < nliterals; ++i )
    {
      int dtype = literal_get_datatype(literals[i]);
      printf("%d %s type = %d", i+1, literals[i], dtype);
      if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
        {
          int ival = literal_to_int(literals[i]);
          printf("  ival = %d", ival);
        }
      else if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
        {
          double dval = literal_to_double(literals[i]);
          printf("  dval = %g", dval);
        }
      else
        {
          printf("  sval = '%s'", literals[i]);
        }
      
      printf("\n");
    }
  
  return 0;
}

#endif
