#include <stdio.h>
#include <string.h>

#include <cdi.h>
#include <cdo.h>
#include <cdo_int.h>
#include "dmemory.h"
#include "field.h"
#include "util.h"


void field_init(field_type *field)
{
  memset(field, 0, sizeof(field_type));
}


field_type **field_allocate(int vlistID, int ptype, int init)
{
  int nvars = vlistNvars(vlistID);

  field_type **field = (field_type **) Malloc(nvars*sizeof(field_type *));

  for ( int varID = 0; varID < nvars; ++varID )
    {
      int nwpv     = vlistInqNWPV(vlistID, varID); // number of words per value; real:1  complex:2
      int gridID   = vlistInqVarGrid(vlistID, varID);
      int gridsize = gridInqSize(gridID);
      int zaxisID  = vlistInqVarZaxis(vlistID, varID);
      int nlevel   = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID, varID);

      field[varID] = (field_type*) Malloc(nlevel*sizeof(field_type));

      for ( int levelID = 0; levelID < nlevel; ++levelID )
	{
	  field_init(&field[varID][levelID]);

	  field[varID][levelID].nwpv     = nwpv;
	  field[varID][levelID].grid     = gridID;
	  field[varID][levelID].size     = gridsize;
	  field[varID][levelID].nsamp    = 0;
	  field[varID][levelID].nmiss    = 0;
	  field[varID][levelID].nmiss2   = 0;
	  if ( ptype & FIELD_FLT ) field[varID][levelID].memtype = MEMTYPE_FLOAT;
	  field[varID][levelID].missval  = missval;
	  field[varID][levelID].ptr      = NULL;
	  field[varID][levelID].ptr2     = NULL;
	  field[varID][levelID].weight   = NULL;

	  if ( ptype & FIELD_PTR )
	    {
              if ( ptype & FIELD_FLT )
                {
                  field[varID][levelID].ptrf = (float*) Malloc(nwpv*gridsize*sizeof(float));
                  if ( init ) memset(field[varID][levelID].ptrf, 0, nwpv*gridsize*sizeof(float));
                }
              else
                {
                  field[varID][levelID].ptr = (double*) Malloc(nwpv*gridsize*sizeof(double));
                  if ( init ) memset(field[varID][levelID].ptr, 0, nwpv*gridsize*sizeof(double));
                }
            }

	  if ( ptype & FIELD_PTR2 )
	    {
              if ( ptype & FIELD_FLT )
                {
                  field[varID][levelID].ptr2 = Malloc(nwpv*gridsize*sizeof(float));
                  if ( init ) memset(field[varID][levelID].ptr2, 0, nwpv*gridsize*sizeof(float));
                }
              else
                {
                  field[varID][levelID].ptr2 = Malloc(nwpv*gridsize*sizeof(double));
                  if ( init ) memset(field[varID][levelID].ptr2, 0, nwpv*gridsize*sizeof(double));
                }
            }

	  if ( ptype & FIELD_WGT )
	    {
	      field[varID][levelID].weight = (double*) Malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].weight, 0, nwpv*gridsize*sizeof(double));
	    }    
	}
    }

  return field;
}


field_type **field_malloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 0);
}


field_type **field_calloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 1);
}


void field_free(field_type **field, int vlistID)
{
  int nvars = vlistNvars(vlistID);
  for ( int varID = 0; varID < nvars; ++varID )
    {
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for ( int levelID = 0; levelID < nlevel; ++levelID )
	{
	  if ( field[varID][levelID].ptr )    Free(field[varID][levelID].ptr);
	  if ( field[varID][levelID].ptrf )   Free(field[varID][levelID].ptrf);
	  if ( field[varID][levelID].ptr2 )   Free(field[varID][levelID].ptr2);
       	  if ( field[varID][levelID].weight ) Free(field[varID][levelID].weight);
	}

      Free(field[varID]);
    }

  Free(field);
}
