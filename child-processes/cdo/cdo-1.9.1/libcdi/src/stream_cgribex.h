#ifndef _STREAM_CGRIBEX_H
#define _STREAM_CGRIBEX_H

int cgribexScanTimestep1(stream_t * streamptr);
int cgribexScanTimestep2(stream_t * streamptr);
int cgribexScanTimestep(stream_t * streamptr);

int cgribexDecode(int memtype, void *gribbuffer, int gribsize, void *data, long datasize,
		  int unreduced, int *nmiss, double missval);

size_t cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg,
		     long datasize, const void *data, int nmiss, void *gribbuffer, size_t gribbuffersize);

void *cgribex_handle_new_from_meassage(void *gribbuffer, size_t recsize);
void cgribex_handle_delete(void *gh);

void cgribexChangeParameterIdentification(void *gh, int code, int ltype, int lev);

#endif  /* _STREAM_CGRIBEX_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
