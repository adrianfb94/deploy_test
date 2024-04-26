/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <cdi.h>
#include "cdo_int.h"
#include "pstream.h"
#include "listarray.h"

// NOTE: All operators in this module works only on GRIB edition 1 files!

#ifdef __cplusplus
extern "C" {
#endif
void streamGrbChangeParameterIdentification(int code, int ltype, int lev);
#if defined (__cplusplus)
}
#endif

/*
Supported notations:
======================
Selection provided on commandline:
---------------------------------
cdo selmulti,'(33/34;105;10)'
 - or -
cdo selmulti,'(33/34;105;10);(11/6;109;55)'
cdo selmulti,'(33/34;105;10);(11/6;109;40/55)'
cdo selmulti,'(*;105;10);(11/6;109;40/55);(*;105;2)'
cdo selmulti,'{(33/34;105;10);(11/32,8;109;51/52/53/54/55)}'

NOTE: ' .. ' are mandatory !

Selection provided from a text file:
---------------------------------

cdo selmulti,selection_10m_wind.txt

(*A*) Compact general notation, selection file content:

(1; 103; 0)
(33,34; 105; 10)
(11,17; 105; 2)
(71,73,74,75,61,62,65,117,67,122,121,11,131,66,84,111,112; 105; 0)
# If nothing <'sel(' or 'del('>  is specified then
# the operator -selmulti or -delmulti decides if it will be selection of extraction or delete
# Explicite select or delete is also possible:
#(11; 109; *)
#(*; 105; *)
#del(*; 109; *)
#sel(*; 105; *)
#sel(*; 100; *)

# BUT simple array arithmetics should be also possible ("*" ~= mulc;  "+' ~= addc)
sel(33,34;105,1000,3000):math(*2;)        # not implemented yet
sel(11;105,500,1500,3000):math(+273.15;)  # not implemented yet

(*B*) HIP.X notation (KNMI specific), selection file content:

SELECT, PARAMETER=1, LEVTYPE=103, LEVEL=0
SELECT, PARAMETER=33/34, LEVTYPE=105, LEVEL=10
SELECT, PARAMETER=11/17, LEVTYPE=105, LEVEL=2
SELECT, PARAMETER=71/73/74/75/61/62/65/117/67/122/121/11/131/66/84/111/112, LEVTYPE=105, LEVEL=0
# Explicite delete is also possible:
#DELETE, PARAMETER=128, LEVTYPE=109, LEVEL=*

# BUT simple array arithmetics should be also possible (SCALE ~= mulc;  OFFSET ~= addc)

# The following will convert Pressure from Pa into HPa; Temp from Kelvin to Celsius:
SELECT, PARAMETER=1, LEVTYPE= 103, LEVEL=0, SCALE=0.01
SELECT, PARAMETER=11, LEVTYPE=105, LEVEL=2, OFFSET=273.15
SELECT, PARAMETER=33/34, LEVTYPE=105, LEVEL=10

If SCALE and/or OFFSET are defined, then the data values are scaled as SCALE*(VALUE-OFFSET).
The default value for SCALE is 1.0; the default for OFFSET is 0.0.

**** changemulti ***********

cdo changemulti,'{(134;1;*|1;105;*);{(6;1;*|6;105;*)};{(246;*;*|76;*;*)};{(247;*;*|58;*;*)};{(248;*;*|71;*;*)}' fileIN fileOUT


cdo changemulti,'{(134;1;*|1;105;*)}' fileIN fileOUT
# surface pressure has ECMWF code; change it into Hirlam notation ..
grib_set -w indicatorOfParameter=134,indicatorOfTypeOfLevel=1 -s indicatorOfParameter=1,indicatorOfTypeOfLevel=105 ECMWF_H11_test0.grb ECMWF_H11_test1.grb

cdo changemulti,'{(6;1;*|6;105;*)}' fileIN fileOUT
# orography has wrong level-type, should be 105
grib_set -w indicatorOfParameter=6,indicatorOfTypeOfLevel=1 -s indicatorOfParameter=6,indicatorOfTypeOfLevel=105 ECMWF_H11_test1.grb ECMWF_H11_test2.grb

cdo changemulti,'{(246;*;*|76;*;*)}' fileIN fileOUT
# change code for cloud_water
grib_set -w indicatorOfParameter=246 -s indicatorOfParameter=76 ECMWF_H11_test2.grb ECMWF_H11_test3.grb

cdo changemulti,'{(247;*;*|58;*;*)}' fileIN fileOUT
# change code for cloud_ice
grib_set -w indicatorOfParameter=247 -s indicatorOfParameter=58 ECMWF_H11_test3.grb ECMWF_H11_test4.grb

cdo changemulti,'{(248;*;*|71;*;*)}' fileIN fileOUT
# change code for total_cloud_cover
grib_set -w indicatorOfParameter=248 -s indicatorOfParameter=71 ECMWF_H11_test4.grb ECMWF_H11_test.grb


*/


typedef struct {
  lista_t *codeLST;
  int  ncodes;

  lista_t *levelTypeLST;
  int  nlevelTypes;

  lista_t *levelLST;
  int  nlevels;
  int  sel_or_del_or_change;  // sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete, 3:change
  int  simpleMath;  // 1:  simple array arithmetics ( *,+), 0: do nothing
  float scale;
  float offset;

  int changedCode;     // used only changemulti mode
  int changedLevelType;
  int changedLevel;
} TUPLEREC;



int push_backIntList(int value, lista_t *list, int arraylen);
int checkListContainsInt(int value, lista_t *list, int arraylen);
int getIntFromList(int index, lista_t *list, int arraylen);

#define MAX_TUPLES 1000
static int NUMTUPLES = 0;
static TUPLEREC *SelTUPLEREC[MAX_TUPLES];


TUPLEREC *TUPLERECNew();
void push_backSelTuple(TUPLEREC *tp);
TUPLEREC * getSelTuple(int index);

void printSelectionTuples();
int getNumberOfSelectionTuples();
int getNumberOfDeleteSelectionTuples();

int multiSelectionParser(const char *filenameOrString);


void *Selmulti(void *argument)
{
  int varID, levelID;
  int nlevs, code, zaxisID;
  int ltype = 0;
  int varID2, levelID2;
  int sellevel, selcode, selltype;
  bool lcopy = false;
  int nmiss;
  int simpleMath=0;  // 1:  simple array arithmetics ( *,+), 0: do nothing
  float scale = 1.0;
  float offset = 0.0; // If SCALE and/or OFFSET are defined, then the data values are scaled as SCALE*(VALUE-OFFSET).
                      // The default value for SCALE is 1.0; the default for OFFSET is 0.0.
  double  missval;

  cdoInitialize(argument);

  // clang-format off
  int SELMULTI    = cdoOperatorAdd("selmulti",    0, 0, "filename/string with selection specification ");
  int DELMULTI    = cdoOperatorAdd("delmulti",    0, 0, "filename/string with selection specification ");
  int CHANGEMULTI = cdoOperatorAdd("changemulti", 0, 0, "filename/string with selection specification ");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  char *filenameOrString;

  //operatorCheckArgc(1);
  filenameOrString = operatorArgv()[0];
  if ( cdoDebugExt )
    {
      printf("Given operator arguments (nr=%d): \n", operatorArgc());
      for (int i=0; i < operatorArgc(); i++)
        printf("%s",operatorArgv()[i]);
      printf("\n");
    }
  if (!multiSelectionParser(filenameOrString))
    cdoWarning("Error processing file with selection description!\n%s",filenameOrString);

  if (operatorID == SELMULTI)
    if (getNumberOfSelectionTuples()==0)
      cdoAbort("Error! You must provide at lease ONE selection tuple!\nNotations: 'SELECT,  .. or sel(/;;) or (/;;)'\nCheck the file: %s",filenameOrString);

  if (operatorID == DELMULTI)
    if (getNumberOfDeleteSelectionTuples()==0)
      cdoAbort("Error! You must provide at lease ONE selection tuple!\nNotations: 'DELETE,  .. or del(/;;) or (/;;)'\nCheck the file: %s",filenameOrString);

  if (operatorID == CHANGEMULTI)
    if (getNumberOfSelectionTuples()==0)
      cdoAbort("Error! You must provide at lease ONE selection tuple!\nNotations: 'CHANGE,  .. or (/;;|;;;)'\nCheck the file: %s",filenameOrString);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  int nvars = vlistNvars(vlistID1);

  if ( cdoDebugExt ) cdoPrint(" Total number of variables: %d", nvars);

  for ( varID = 0; varID < nvars; varID++ )
    {
      code    = vlistInqVarCode(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      ltype = zaxis2ltype(zaxisID);
      nlevs   = zaxisInqSize(zaxisID);

      for ( levelID = 0; levelID < nlevs; levelID++ )
        {
          double level = zaxisInqLevel(zaxisID, levelID);

          if (operatorID == DELMULTI) vlistDefFlag(vlistID1, varID, levelID, TRUE); // set initially, override bellow if in selection
          if (operatorID == CHANGEMULTI)
            {
              vlistDefFlag(vlistID1, varID, levelID, TRUE); // change operation copies all fields
              continue;
            }

          for ( int ii=0; ii<NUMTUPLES; ii++ )
            {
              TUPLEREC *tuplerec = getSelTuple(ii);
              //if ( cdoDebugExt ) cdoPrint(" Processing: (code %d, ltype %d, level %d);  nvars=%d, varID=%d", code, ltype, (int)level, nvars, varID);
              // Note: When the list is Empty then function checkListContainsInt() also returns true !
              selcode  = checkListContainsInt(code, tuplerec->codeLST, tuplerec->ncodes);
              selltype = checkListContainsInt(ltype, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
              sellevel = checkListContainsInt((int)level, tuplerec->levelLST, tuplerec->nlevels);
              if ( selcode && selltype && sellevel )
                {
                  if (operatorID == SELMULTI)
                    {
                      switch (tuplerec->sel_or_del_or_change)
                        {
                        case 0:   // operator decides ...
                          vlistDefFlag(vlistID1, varID, levelID, TRUE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        case 1:
                          vlistDefFlag(vlistID1, varID, levelID, TRUE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        case 2:
                          vlistDefFlag(vlistID1, varID, levelID, FALSE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        }
                    }
                  else if (operatorID == DELMULTI)
                    {
                      switch (tuplerec->sel_or_del_or_change)
                        {
                        case 0:   // operator decides ...
                          vlistDefFlag(vlistID1, varID, levelID, FALSE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        case 1:
                          vlistDefFlag(vlistID1, varID, levelID, TRUE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        case 2:
                          vlistDefFlag(vlistID1, varID, levelID, FALSE);
                          if ( cdoDebugExt )
                            {
                              if (!tuplerec->simpleMath)
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", code, ltype, (int)(level), varID, levelID);
                              else
                                cdoPrint(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f", code, ltype, (int)(level), varID, levelID, tuplerec->scale,tuplerec->offset);
                            }
                          break;
                        }
                    }
                  break;
                }
            } //end for ( .. NUMTUPLES
        } //end for ( levelID
    } // end for ( varID

  if ( cdoDebugExt ) cdoPrint(" Writing the selected fields ...");

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);

  nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; ++varID )
    if ( vlistInqVarTimetype(vlistID2, varID) != TIME_CONSTANT ) break;
  if ( varID == nvars ) vlistDefNtsteps(vlistID2, 0);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nrecs = vlistNrecs(vlistID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double *) malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          missval = vlistInqVarMissval(vlistID1, varID);

          if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
            {
              simpleMath=0;  // 1:  simple array arithmetics ( *,+), 0: do nothing
              scale = 1.0;
              offset = 0.0;
              code    = vlistInqVarCode(vlistID1, varID);
              zaxisID = vlistInqVarZaxis(vlistID1, varID);
              double level = zaxisInqLevel(zaxisID, levelID);
              ltype = zaxis2ltype(zaxisID);
              for ( int ii=0; ii<NUMTUPLES; ii++ )
                {
                  TUPLEREC *tuplerec = getSelTuple(ii);
                  // Note: When the list is Empty then function checkListContainsInt() also returns true !
                  selcode  = checkListContainsInt(code, tuplerec->codeLST, tuplerec->ncodes);
                  selltype = checkListContainsInt(ltype, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
                  sellevel = checkListContainsInt((int)level, tuplerec->levelLST, tuplerec->nlevels);
                  lcopy = true;
                  if ( selcode && selltype && sellevel )
                    {
                      if ( operatorID == CHANGEMULTI )
                        {
                          if ( cdoDebugExt )
                            cdoPrint(" Processing: (code %d, ltype %d, level %d);  nvars=%d, varID=%d => (selcode %d, selltype %d, sellevel %d) => change (%d,%d,%d)",
                                     code, ltype, (int)level, nvars, varID, selcode, selltype, sellevel, tuplerec->changedCode, tuplerec->changedLevelType, tuplerec->changedLevel);
                          if ( (tuplerec->changedCode==-1) && (tuplerec->changedLevelType==-1) && (tuplerec->changedLevel==-1) )
                            cdoPrint(" WARNING: Cannot CHANGE identification!");
                          else
                            streamGrbChangeParameterIdentification(tuplerec->changedCode, tuplerec->changedLevelType, tuplerec->changedLevel);
                          // Calling PROXY function streamGrbChangeParameterIdentification()
                          // which results in later calling func. gribapiChangeParameterIdentification(); see stream_gribapi.c
                          // The change happens during the last step of writing a grib-record into a file.
                        }
                      else
                        {
                          if ( cdoDebugExt ) cdoPrint(" Processing: (code %d, ltype %d, level %d);  nvars=%d, varID=%d => (selcode %d, selltype %d, sellevel %d)", code, ltype, (int)level, nvars, varID, selcode, selltype, sellevel);
                        }
                      simpleMath = tuplerec->simpleMath;  // 1:  simple array arithmetics ( *,+), 0: do nothing
                      if (simpleMath)
                        {
                          scale = tuplerec->scale;
                          offset = tuplerec->offset;
                          lcopy = false;
                        }
                      break; // get out of this for loop
                    }
                } // end of for (ii=0; ii<NUMTUPLES ..

              varID2   = vlistFindVar(vlistID2, varID);
              levelID2 = vlistFindLevel(vlistID2, varID, levelID);

              // tijdelijk PATCH M.K.
              if ((varID2==-1) || (levelID2==-1)) {
                cdoPrint(" Warning: Missing varID or levelID with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)] .. #2[varID(%d),levelID(%d)]",code, ltype, (int)(level), varID,levelID, varID2, levelID2);
                continue;
              }
              pstreamDefRecord(streamID2, varID2, levelID2);

              if ( lcopy )
                {
                  if ( cdoDebugExt ) cdoPrint(" Copying record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]",recID, code, ltype, (int)(level), varID,levelID);
                  pstreamCopyRecord(streamID2, streamID1);
                }
              else
                {
                  pstreamReadRecord(streamID1, array, &nmiss);

                  if (!simpleMath)
                    {
                      if ( cdoDebugExt ) cdoPrint(" Writing record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]",recID, code, ltype, (int)(level), varID,levelID);
                      pstreamWriteRecord(streamID2, array, nmiss);
                    }
                  else  // 1:  simple array arithmetics ( *,+)
                    {
                      if ( cdoDebugExt ) cdoPrint(" Writing record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f",recID, code, ltype, (int)(level), varID,levelID,scale,offset);
                      for ( int li = 0; li < gridsize; ++li )
                        if (! DBL_IS_EQUAL(array[li], missval) )
                          {
                            array[li] = scale*(array[li] - offset);
                            // If SCALE and/or OFFSET are defined, then the data values are scaled as SCALE*(VALUE-OFFSET).
                          }
                      pstreamWriteRecord(streamID2, array, nmiss);
                    }
                } // end of else ( lcopy )
            } // end if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
        } // end for ( recID ..

      /*if ( vlistInqFlag(vlistID1, varID, levelID) == FALSE )
        if (operatorID == DELMULTI)
        {
        code    = vlistInqVarCode(vlistID1, varID);
        zaxisID = vlistInqVarZaxis(vlistID1, varID);
        level = zaxisInqLevel(zaxisID, levelID);
        ltype = zaxis2ltype(zaxisID);
        if ( cdoDebugExt ) cdoPrint(" Removing record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]",recID, code, ltype, (int)(level), varID, levelID);
        }
        else*/

      tsID++;
    } // end while ( (nrecs

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  cdoDebugExt = 0;

  return 0;
}


TUPLEREC *TUPLERECNew()
{
  TUPLEREC *tpl = (TUPLEREC *) malloc(sizeof(TUPLEREC));

  tpl->codeLST = lista_new(INT_LISTA);
  tpl->ncodes = 0;

  tpl->levelTypeLST = lista_new(INT_LISTA);
  tpl->nlevelTypes = 0;

  tpl->levelLST = lista_new(INT_LISTA);
  tpl->nlevels = 0;
  tpl->sel_or_del_or_change = 0; // sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete, 3:change

  tpl->simpleMath= 0;  // 1:  simple array arithmetics ( *,+), 0: do nothing
  tpl->scale = 1.0;
  tpl->offset = 0.0;

  tpl->changedCode = -1;    // used only changemulti mode
  tpl->changedLevelType = -1;
  tpl->changedLevel = -1;

  return tpl;
}

void push_backSelTuple(TUPLEREC *tp)
{
 if (NUMTUPLES<MAX_TUPLES) {
  SelTUPLEREC[NUMTUPLES] = tp;
  NUMTUPLES++;
 }
}

TUPLEREC * getSelTuple(int index)
{
 if (index<NUMTUPLES) {
  return SelTUPLEREC[index];
 }
 return NULL;
}

int push_backIntList(int value, lista_t *list, int arraylen)
{
  lista_set_int(list, arraylen, value);
  arraylen +=1;
  return arraylen;
}

int checkListContainsInt(int value, lista_t *list, int arraylen)
// Note: When the list is Empty it also returns true !
{
 int i;
 int *array = (int *) lista_dataptr(list);

 if (arraylen==0) {
        //if ( cdoDebugExt ) cdoPrint(" value=%d: found. (list empty)");
        return 1;
 }
 for (i=0; i<arraylen; i++)
 {
   if (array[i]==-1) {  // this is for '*' selection; can be any code, any level or any level-type
        //if ( cdoDebugExt ) cdoPrint(" value=%d: found.");
        return 1;
   }

   if (array[i]==value) {
        //if ( cdoDebugExt ) cdoPrint(" value=%d: found.");
        return 1;
   }
 }
 //if ( cdoDebugExt ) cdoPrint(" value=%d: NOT found.");
 return 0;
}


int getIntFromList(int index, lista_t *list, int arraylen)
{
  int *array = (int *) lista_dataptr(list);
  if (index <arraylen)
    return array[index];
  else
    return -999;
}





#define MAX_LINE_LEN 65536

static
char *removeSpaces(char *pline)
{
  if (pline==NULL) return NULL;
  while ( isspace((int) *pline) ) pline++;
  return pline;
}


static
char *skipSeparator(char *pline)
{
  if (pline==NULL) return NULL;
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' || *pline == '/' || *pline == ',' ) pline++;
  while ( isspace((int) *pline) ) pline++;
  return pline;
}


static
char *goToNextSeparator(char *pline)
{
  if (pline==NULL) return NULL;
  int separatorFound =0;
  while ( isspace((int) *pline) || !separatorFound)
  {
        if (*pline=='\0') return NULL;
        pline++;
        if (*pline == '|') return NULL;
        if ( *pline == '=' || *pline == ':' || *pline == '/' || *pline == ',' )
                separatorFound = 1;
        else
            if ( *pline == ';')
            { pline++;
              pline = removeSpaces(pline);
              return (pline);
            }
  }
  if (separatorFound) pline++;
  if ( cdoDebugExt>=100 ) cdoPrint("goToNextSeparator():  pline= ('%s') ", pline);
  //while ( isspace((int) *pline) ) pline++;
  pline = removeSpaces(pline);
  return pline;
}


static
char *strContains(char *str, const char *substr)
{
  if (str==NULL) return NULL;
  if (substr==NULL) return NULL;

  str = removeSpaces(str);

  size_t lensub = strlen(substr);
  size_t lenstr = strlen(str);

  if (lensub>lenstr) {
        if ( cdoDebugExt>=100 ) cdoPrint("strContains():  substr('%s') NOT found in str('%s');  lensub(%d)>lenstr(%d) ", substr,str, lensub,lenstr);
        return NULL;
  }
  char *rv = strstr(str, substr);
  if (rv) {
        if ( cdoDebugExt>=100 ) cdoPrint("strContains():  substr('%s') FOUND in str('%s')", substr,str);
        return (rv+lensub);   // points after subStr ..
  }
  else {
        if ( cdoDebugExt>=100 ) cdoPrint("strContains():  substr('%s') NOT found in str('%s')", substr,str);
        return rv;
 }
}

static
char *findParamEnd(char *str)
{
  char *ptr = str;
  char *ptrEnding=NULL;

  if (str==NULL) return NULL;
  // supported endings are: ", " or ";"
  if (ptrEnding==NULL) ptrEnding = strContains(ptr,", "); // HIP notation
  if (ptrEnding==NULL) ptrEnding = strContains(ptr,";");  // compact notation
  if (ptrEnding!=NULL) {
          ptrEnding = removeSpaces(ptrEnding);
          if ( cdoDebugExt>=100 ) cdoPrint(" ptrEnding='%s'", ptrEnding);
          return ptrEnding;
  }
  if ( cdoDebugExt>=100 ) cdoPrint(" ptrEnding=end-of-string");
  size_t lenstr = strlen(str);
  ptrEnding = str+lenstr;
  ptrEnding = removeSpaces(ptrEnding);
  return ptrEnding;
}

static
char *findTupleEnd(char *str)
{
  char *ptr = str;
  char *ptrEnding;

  if (str==NULL) return NULL;
  // supported endings are: ")" or ");"
  ptrEnding = strContains(ptr,")");
  if (ptrEnding==NULL) ptrEnding = strContains(ptr,");");
  if (ptrEnding!=NULL) {
          ptrEnding = removeSpaces(ptrEnding);
          if ( cdoDebugExt>=100 ) cdoPrint(" findTupleEnd='%s'", ptrEnding);
          return ptrEnding;
  }
  if ( cdoDebugExt>=100 ) cdoPrint(" findTupleEnd=end-of-string");
  size_t lenstr = strlen(str);
  ptrEnding = str+lenstr;
  ptrEnding = removeSpaces(ptrEnding);
  return ptrEnding;
}

static
char *readlineForParsing(FILE *gfp, char *strToParsePtr, char *line)
{
  if ( gfp != NULL ) // file is open => we parse text from a file
    {
      int status = readline(gfp, line, MAX_LINE_LEN);
      return (status == 0) ? NULL : line;
    }
  else
    if (strToParsePtr!=NULL)  // we parse a given string
      {
        if ( cdoDebugExt>=30 ) cdoPrint("%s(): Parsing selection string:  %s", __func__, strToParsePtr);
        char *tpEnd = NULL;
        if (strlen(strToParsePtr)>0)
          tpEnd = findTupleEnd(strToParsePtr);
        if (tpEnd==NULL)
          {
            if ( cdoDebugExt>=100 ) cdoPrint("%s(): End of selection string reached.", __func__);
            return NULL;
          }
        else
          {
            tpEnd[0] = 0;
            if (strlen(strToParsePtr)<=MAX_LINE_LEN)
              strcpy(line,strToParsePtr);
            if ( cdoDebugExt>=100) cdoPrint("%s(): Current selection line=%s", __func__, line);
            strToParsePtr = tpEnd +1;
            return strToParsePtr;
          }
      }
    else
      {
        cdoAbort(" Cannot parse selection string:  %s", strToParsePtr);
        return NULL;
      }
}

int multiSelectionParser(const char *filenameOrString)
{
  char line[MAX_LINE_LEN], *pline;
  char *strpos;
  char *parEnd;
  int val;
  float floatval;
  TUPLEREC *tuplerec;
  int selectionRec;
  FILE *gfp = NULL;
  char strToParse[MAX_LINE_LEN];
  char *strToParsePtr = NULL; strToParse[0]=0;
  int convertLevelsHPa2Pa = 0;  // this is HIP - KNMI selection specific
  char first3chars[4];
  strncpy(first3chars,filenameOrString,3);
  
  if ( (filenameOrString[0]=='{') || (filenameOrString[0]=='(') || strContains(first3chars, "del") || strContains(first3chars, "sel") )
  {
      // cdo selmulti,'(33/34;105;10)'
      // - or -
      // cdo selmulti,'{(33/34;105;10);(11/32,8;109;51/52/53/54/55)}'
      // - or -
      // cdo changemulti,'{(134;1;0|1;105;0);(6;1;0|6;105;0)}'
      strncpy(strToParse,filenameOrString,MAX_LINE_LEN);
      strToParsePtr = &strToParse[0];
      if (strToParsePtr[0]=='{')
          strToParsePtr++;
      int strLn = strlen(strToParsePtr);
      if (strToParsePtr[strLn-1]=='}')
        strToParsePtr[strLn-1] = 0;
      if ( cdoDebugExt ) cdoPrint(" Parsing selection string:  %s", strToParsePtr);
  }
  else
  {
      gfp = fopen(filenameOrString, "r");
      if ( cdoDebugExt ) cdoPrint(" Parsing file:  %s", filenameOrString);
      if (gfp==NULL)
      {
        cdoAbort(" Missing file:  %s", filenameOrString);
        return 0;
      }
  }

  while ( (strToParsePtr = readlineForParsing(gfp, strToParsePtr, line)) )
  {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      pline = line;
      if ( cdoDebugExt>=30 ) cdoPrint(": Line: %s", pline);
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      strpos=strContains(pline, "SELECT, ");
      selectionRec = 0; // default is 0;  sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete, 3:change
      if ( strpos != NULL ) selectionRec = 1;
      else {
               strpos=strContains(pline, "DELETE, ");
              if ( strpos != NULL ) selectionRec = 2;
              else {
                      strpos=strContains(pline, "REMOVE, ");
                      if ( strpos != NULL ) selectionRec = 2;
                   }
            }
      if (strpos!= NULL ) // we have SELECT ..
      {
          if ( cdoDebugExt ) cdoPrint(" Parsing notation SELECT: %s", strpos);
          pline = strpos;
          tuplerec = TUPLERECNew();
          tuplerec->sel_or_del_or_change = selectionRec;
          push_backSelTuple(tuplerec);
          // SELECT, PARAMETER=11/17, LEVTYPE=105, LEVEL=2
          while ( (pline!=NULL) && (strlen(pline)!=0) )
          {     if ( cdoDebugExt>=100 ) cdoPrint("pline='%s'", pline);
                strpos=strContains(pline, "PARAMETER=");
                if ( (strpos) != NULL )
                {
                  if ( cdoDebugExt>=100 ) cdoPrint(": PARAMETER=%s", strpos);
                  pline = strpos;
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("code=%d",val);
                        tuplerec->ncodes = push_backIntList(val, tuplerec->codeLST, tuplerec->ncodes);
                        strpos = goToNextSeparator(pline);
                        pline = strpos;
                        if (!strpos) break;
                        pline = skipSeparator(strpos);
                  }
                }
                if ( cdoDebugExt>=100 ) cdoPrint("pline='%s'", pline);
                strpos=strContains(pline, "LEVTYPE=");
                if ( (strpos) != NULL )
                {
                  if ( cdoDebugExt>=100 ) cdoPrint(": LEVTYPE=%s", strpos);
                  pline = strpos;
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd) && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if (val==100)
                        {  convertLevelsHPa2Pa = 1;
                           if ( cdoDebugExt>=1 ) cdoPrint("Detected levelType=100 ! HIP selection specifies HPa but CDO has levels in Pascals.\nSelection HPa's will be converted into Pa's.");
                        } else convertLevelsHPa2Pa = 0;
                        if ( cdoDebugExt>=100 ) cdoPrint("levelType=%d",val);
                        tuplerec->nlevelTypes = push_backIntList(val, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
                        strpos = goToNextSeparator(pline);
                        pline = strpos;
                        if (!strpos) break;
                        pline = skipSeparator(strpos);
                  }
                }
                if ( cdoDebugExt>=100 ) cdoPrint("pline='%s'", pline);
                strpos=strContains(pline, "LEVEL=");
                if ( (strpos) != NULL )
                {
                  if ( cdoDebugExt>=100 ) cdoPrint(": LEVEL=%s", strpos);
                  pline = strpos;
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if ((convertLevelsHPa2Pa) && (val!=-1) )
                        {
                            val*=100;
                        }
                        if ( cdoDebugExt>=100 ) cdoPrint("level=%d",val);
                        tuplerec->nlevels = push_backIntList(val, tuplerec->levelLST, tuplerec->nlevels);
                        strpos = goToNextSeparator(pline);
                        pline = strpos;
                        if (!strpos) break;
                        pline = skipSeparator(strpos);
                  }
                }
                if ( cdoDebugExt>=100 ) cdoPrint("pline='%s' (check SCALE=...)", pline);
                strpos=strContains(pline, "SCALE=");
                if ( (strpos) != NULL )
                {
                  if ( cdoDebugExt>=100 ) cdoPrint(": SCALE= %s", strpos);
                  pline = strpos;
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        floatval = atof(pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("scale=%f",floatval);
                        tuplerec->simpleMath=1; // 1:  simple array arithmetics ( *,+), 0: do nothing
                        tuplerec->scale = floatval;  //tuplerec->offset = 0.0;
                        strpos = goToNextSeparator(pline);
                        pline = strpos;
                        if (!strpos) break;
                        pline = skipSeparator(strpos);
                  }
                }
                if ( cdoDebugExt>=100 ) cdoPrint("pline='%s' (check OFFSET=...)", pline);
                strpos=strContains(pline, "OFFSET=");
                if ( (strpos) != NULL )
                {
                  if ( cdoDebugExt>=100 ) cdoPrint(": OFFSET= %s", strpos);
                  pline = strpos;
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        floatval = atof(pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("offset=%f",floatval);
                        tuplerec->simpleMath=1; // 1:  simple array arithmetics ( *,+), 0: do nothing
                        tuplerec->offset = floatval; // tuplerec->scale = 1.0;
                        strpos = goToNextSeparator(pline);
                        pline = strpos;
                        if (!strpos) break;
                        pline = skipSeparator(strpos);
                  }
                }
          } // end while pline
          continue;
      } // end if SELECT,

      // Here comes the short notation
      selectionRec = 0; // default is 0;  sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete, 3:change
      strpos=strContains(pline, "sel(");
      if ( strpos != NULL ) selectionRec = 1;
      else {
              strpos=strContains(pline, "del(");
              if ( strpos != NULL ) selectionRec = 2;
              else {
                  strpos=strContains(pline, "|");
                  if ( strpos != NULL ) selectionRec = 3;
                  else {
                          strpos=strContains(pline, "(");
                          if ( strpos != NULL ) selectionRec = 0;
                      }
              }
           }
      if ( strpos != NULL )
        {
          if ( cdoDebugExt ) cdoPrint(" Parsing notation (code(s),..; levelType(s),..; level(s),..) : %s; [selectionRec =%d]", strpos, selectionRec );
          if (selectionRec == 3)
              strpos = strContains(pline, "(");
          pline = strpos;
          tuplerec = TUPLERECNew();
          tuplerec->sel_or_del_or_change = selectionRec;
          push_backSelTuple(tuplerec);
          // (33/34; 105; 10)
          while ( (pline!=NULL) && (strlen(pline)!=0) && (pline[0]!=')') )
          {     if ( cdoDebugExt>=100 ) cdoPrint("[1]: pline='%s'", pline);
                // 1st is code
                {
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("code=%d",val);
                        tuplerec->ncodes = push_backIntList(val, tuplerec->codeLST, tuplerec->ncodes);
                        strpos = goToNextSeparator(pline);
                        if (!strpos) { pline=parEnd; break;} else pline = strpos;
                        pline = skipSeparator(strpos);
                  }
                }
                // 2nd is level type
                if ( cdoDebugExt>=100 ) cdoPrint("[2]: pline='%s'", pline);
                {
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if (val==100)
                        {  convertLevelsHPa2Pa = 1;
                           if ( cdoDebugExt>=1 ) cdoPrint("Detected levelType=100 ! Selection specifies HPa but CDO has levels in Pascals.\nSelection HPa's will be converted into Pa's.");
                        } else convertLevelsHPa2Pa = 0;

                        if ( cdoDebugExt>=100 ) cdoPrint("levelType=%d",val);
                        tuplerec->nlevelTypes  =push_backIntList(val, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
                        strpos = goToNextSeparator(pline);
                        if (!strpos) { pline=parEnd; break;} else pline = strpos;
                        pline = skipSeparator(strpos);
                  }
                }
                // 3rd is level
                if ( cdoDebugExt>=100 ) cdoPrint("[3]: pline='%s'", pline);
                {
                  pline = removeSpaces(pline);
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if ((convertLevelsHPa2Pa) && (val!=-1) )
                        {
                            val*=100;
                        }
                        if ( cdoDebugExt>=100 ) cdoPrint("level=%d",val);
                        tuplerec->nlevels = push_backIntList(val, tuplerec->levelLST, tuplerec->nlevels);
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        strpos = goToNextSeparator(pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if (!strpos) {
                            strpos = strContains(pline,"|");  // compact notation for  changemulti
                            if (strpos)  pline=strpos-1; // strContains returns character after...
                                   else  pline=parEnd;
                            if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                            break;
                            }
                        else pline = strpos;
                        pline = skipSeparator(strpos);
                  }
                }

                // OPTIONAL:
                // cdo changemulti,'{(134;1;*|1;105;*);{(6;1;*|6;105;*)};{(246;*;*|76;*;*)};{(247;*;*|58;*;*)};{(248;*;*|71;*;*)}' fileIN fileOUT
                if ( cdoDebugExt>=100 ) cdoPrint("[OPT]: pline='%s'", pline);
                {
                  pline = removeSpaces(pline);
                  // pline points to: "=1;105;*);....."
                  if (pline[0]=='|')
                  {     // changemulti specification
                        tuplerec->sel_or_del_or_change = 3;
                        pline =  &pline[1];
                        // Get changedCode:
                        parEnd = findParamEnd(pline);
                        if ( (!parEnd) || (pline[0]==0))
                                cdoAbort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        tuplerec->changedCode = val;
                        if ( cdoDebugExt>=100 ) cdoPrint("changedCode=%d",val);
                        strpos = goToNextSeparator(pline);
                        if (!strpos)
                            cdoAbort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        pline = skipSeparator(strpos);
                        // Get changedLevelType:
                        parEnd = findParamEnd(pline);
                        if ( (!parEnd) || (pline[0]==0))
                                cdoAbort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        tuplerec->changedLevelType = val;
                        if ( cdoDebugExt>=100 ) cdoPrint("changedLevelType=%d",val);
                        strpos = goToNextSeparator(pline);
                        if (!strpos)
                            cdoAbort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        pline = skipSeparator(strpos);
                        // Get changedLevel:
                        parEnd = findParamEnd(pline);
                        if ( (!parEnd) || (pline[0]==0))
                                cdoAbort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        tuplerec->changedLevel = val;
                        if ( cdoDebugExt>=100 ) cdoPrint("changedLevel=%d",val);
                        pline=parEnd;
                  }// changemulti specification
                  parEnd = findParamEnd(pline);
                  if ((!parEnd) && (!pline))
                        if ( cdoDebugExt>=100 ) cdoPrint("strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd,pline);
                  while ( (pline!=parEnd)  && (strlen(pline)>0))
                  {
                        pline = removeSpaces(pline);
                        if (pline[0]=='*') val=-1; else val = atoi(pline);
                        if ((convertLevelsHPa2Pa) && (val!=-1) )
                        {
                            val*=100;
                        }
                        if ( cdoDebugExt>=100 ) cdoPrint("level=%d",val);
                        tuplerec->nlevels = push_backIntList(val, tuplerec->levelLST, tuplerec->nlevels);
                        strpos = goToNextSeparator(pline);
                        if (!strpos) { pline=parEnd; break;} else pline = strpos;
                        pline = skipSeparator(strpos);
                  }
                }


          } // end while pline
       } // end if "("
  } // end while ( readline(gfp, line, MAX_LINE_LEN) )
  if (gfp!=NULL)
      fclose(gfp);

  printSelectionTuples();
  return 1;

}


void printSelectionTuples()
{
  if ( cdoDebugExt ) cdoPrint(" Printing selection tuples:");

  int ii, ri;
  for (ii=0; ii<NUMTUPLES; ii++)
  {
        TUPLEREC *tuplerec = getSelTuple(ii);
        char strval[1000]; strcpy(strval,"(");
        char bff[200]; bff[0]='\0';
        if ( cdoDebugExt ) cdoPrint(" Selection tuple [%d]: ncodes=%d, nlevelTypes=%d, nlevels=%d",ii,tuplerec->ncodes,tuplerec->nlevelTypes,tuplerec->nlevels);
        for (ri=0; ri< tuplerec->ncodes; ri++)
        {
                 sprintf(bff,"%d", lista_get_int(tuplerec->codeLST, ri) );
                 strcat(strval,bff);
                 if ( (ri+1)< tuplerec->ncodes )
                        strcat(strval,"/");
                 else    strcat(strval,";");
        }
        for (ri=0; ri< tuplerec->nlevelTypes; ri++)
        {
                 sprintf(bff,"%d", lista_get_int(tuplerec->levelTypeLST, ri) );
                 strcat(strval,bff);
                 if ( (ri+1)< tuplerec->nlevelTypes )
                        strcat(strval,"/");
                 else    strcat(strval,";");
        }
        for (ri=0; ri< tuplerec->nlevels; ri++)
        {
                 sprintf(bff,"%d", lista_get_int(tuplerec->levelLST, ri) );
                 strcat(strval,bff);
                 if ( (ri+1)< tuplerec->nlevels )
                        strcat(strval,"/");
                 else    strcat(strval,")");
        }


        if (tuplerec->simpleMath)
        {
            sprintf(bff," {scale = %f; offset = %f}",tuplerec->scale, tuplerec->offset);
            strcat(strval,bff);
        }

        // sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete, 3:change
        if (tuplerec->sel_or_del_or_change==1)
        {
            if ( cdoDebugExt ) cdoPrint(" Selection tuple [%d] = %s (select)",ii,strval);
        }
        else
        if (tuplerec->sel_or_del_or_change==2)
        {
            if ( cdoDebugExt ) cdoPrint(" Selection tuple [%d] = %s (delete)",ii,strval);
        }
        if (tuplerec->sel_or_del_or_change==3)
        {
            if ( cdoDebugExt ) cdoPrint(" Selection tuple [%d] = %s (change) => (%d; %d;%d;)",ii,strval,tuplerec->changedCode,tuplerec->changedLevelType,tuplerec->changedLevel);
        }
        else
            cdoPrint(" Selection tuple [%d] = %s (select/delete)",ii,strval);
  }
}


int getNumberOfSelectionTuples()
{
  int ii;
  int nn=0;
  for (ii=0; ii<NUMTUPLES; ii++)
  {
        TUPLEREC *tuplerec = getSelTuple(ii);
        if (tuplerec->sel_or_del_or_change!=2) nn++;
  }
  return nn;
}

int getNumberOfDeleteSelectionTuples()
{
  int ii;
  int nn=0;
  cdoPrint("NUMTUPLES = %d",NUMTUPLES);
  printSelectionTuples();
  for (ii=0; ii<NUMTUPLES; ii++)
  {
        TUPLEREC *tuplerec = getSelTuple(ii);
        cdoPrint("getNumberOfDeleteSelectionTuples() [%d]=%d\n", ii, tuplerec->sel_or_del_or_change);
        if (tuplerec->sel_or_del_or_change!=1) nn++;
  }
  return nn;
}
