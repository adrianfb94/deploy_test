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

#ifndef PSTREAM_WRITE_H
#define PSTREAM_WRITE_H
#include "argument.h"

int     pstreamOpenWrite(const argument_t *argument, int filetype);
void    pstreamClose(int pstreamID);

void    pstreamDefVlist(int pstreamID, int vlistID);

void    pstreamDefTimestep(int pstreamID, int tsID);

void    pstreamDefRecord(int pstreamID, int  varID, int  levelID);

void    pstreamWriteRecord(int pstreamID, double *data, int nmiss);
void    pstreamWriteRecordF(int pstreamID, float *data, int nmiss);

#endif  /* PSTREAM_WRITE_H */
