#ifndef _STREAM_GRB_H
#define _STREAM_GRB_H

static inline bool gribbyte_get_bit(int number, int bit) { return (bool)((number >> (8-bit)) & 1); }
static inline void gribbyte_set_bit(int *number, int bit) { *number |= 1 << (8-bit); }
static inline void gribbyte_clear_bit(int *number, int bit) { *number &= ~(1 << (8-bit)); }

int   grbBitsPerValue(int datatype);

int   grbInqContents(stream_t *streamptr);
int   grbInqTimestep(stream_t *streamptr, int tsID);

int   grbInqRecord(stream_t *streamptr, int *varID, int *levelID);
void  grbDefRecord(stream_t *streamptr);
void  grb_read_record(stream_t *streamptr, int memtype, void *data, int *nmiss);
void  grb_write_record(stream_t *streamptr, int memtype, const void *data, int nmiss);
void  grbCopyRecord(stream_t *streamptr2, stream_t *streamptr1);

void  grb_read_var(stream_t *streamptr, int varID, int memtype, void *data, int *nmiss);
void  grb_write_var(stream_t *streamptr, int varID, int memtype, const void *data, int nmiss);

void  grb_read_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, void *data, int *nmiss);
void  grb_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, int nmiss);

int   grib1ltypeToZaxisType(int grib_ltype);
int   grib2ltypeToZaxisType(int grib_ltype);

int   zaxisTypeToGrib1ltype(int zaxistype);
int   zaxisTypeToGrib2ltype(int zaxistype);

struct cdiGribParamChange
{
  int code, ltype, lev;
  bool active;
};

struct cdiGribModeChange
{
  bool mode;
  bool active;
};

struct cdiGribScanModeChange
{
  int value;
  bool active;
};

extern struct cdiGribParamChange cdiGribChangeParameterID;
extern struct cdiGribModeChange cdiGribChangeModeUvRelativeToGrid;
extern struct cdiGribScanModeChange cdiGribDataScanningMode;

// Used in CDO
void streamGrbChangeParameterIdentification(int code, int ltype, int lev);
void streamGrbChangeModeUvRelativeToGrid(int mode);
void streamGrbDefDataScanningMode(int scanmode);
int  streamGrbInqDataScanningMode(void);

#endif  /* _STREAM_GRB_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
