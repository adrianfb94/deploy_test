#include <cdi.h>
#include "cdo_int.h"


bool *cdo_read_timestepmask(const char *maskfile, int *n)
{
  *n = 0;

  int streamID = streamOpenRead(maskfile);
  if ( streamID == CDI_UNDEFID ) cdoAbort("Open failed on %s!", maskfile);

  int vlistID = streamInqVlist(streamID);

  int nvars = vlistNvars(vlistID);
  if ( nvars > 1 ) cdoAbort("timestepmask %s contains more than one variable!", maskfile);

  int gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));
  if ( gridsize > 1 ) cdoAbort("timestepmask %s has more than one gridpoint!", maskfile);

  int nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if ( nlev > 1 ) cdoAbort("timestepmask %s has more than one level!", maskfile);

  int nts = vlistNtsteps(vlistID);
  if ( nts == -1 )
    {
      nts = 0;
      while ( streamInqTimestep(streamID, nts) ) nts++;

      if ( cdoVerbose ) cdoPrint("%s: counted %i timeSteps in %s", __func__, nts, maskfile);

      streamClose(streamID);
      streamID = streamOpenRead(maskfile);

    }
  else
    if ( cdoVerbose ) cdoPrint("%s: found %i timeSteps in %s", __func__, nts, maskfile);

  *n = nts;
  bool *imask = (bool*) Malloc(nts*sizeof(bool));

  int nrecs;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      if ( nrecs != 1 ) cdoAbort("Internal error; unexprected number of records!");

      int varID, levelID;
      int nmiss;
      double value;
      streamInqRecord(streamID, &varID, &levelID);
      streamReadRecord(streamID, &value, &nmiss);
      
      imask[tsID] = !(nmiss || IS_EQUAL(value, 0));
      
      tsID++;  
    }
  
  streamClose(streamID);

  return imask;
}


bool *cdo_read_mask(const char *maskfile, int *n)
{
  *n = 0;

  int streamID = streamOpenRead(maskfile);
  if ( streamID == CDI_UNDEFID ) cdoAbort("Open failed on %s!", maskfile);

  int vlistID = streamInqVlist(streamID);

  int nvars = vlistNvars(vlistID);
  if ( nvars > 1 ) cdoAbort("Mask %s contains more than one variable!", maskfile);

  int gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));

  int nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if ( nlev > 1 ) cdoAbort("Mask %s has more than one level!", maskfile);

  *n = gridsize;
  bool *imask = (bool*) Malloc(gridsize*sizeof(bool));
  double *dmask = (double*) Malloc(gridsize*sizeof(double));

  int nrecs = streamInqTimestep(streamID, 0);
  if ( nrecs != 1 ) cdoAbort("Internal error; unexprected number of records!");

  int varID, levelID;
  int nmiss;
  streamInqRecord(streamID, &varID, &levelID);
  streamReadRecord(streamID, dmask, &nmiss);

  for ( int i = 0; i < gridsize; ++i )
    imask[i] = IS_NOT_EQUAL(dmask[i], 0);
      
      Free(dmask);
  
  streamClose(streamID);

  return imask;
}
