#include <cdi.h>
#include <signal.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include <unistd.h>
#include "uthash.h"
#include "util.h"

#ifdef __cplusplus
  extern "C" {
#endif
#include "cmor.h"
#ifdef __cplusplus
  }
#endif

#include "netcdf.h"
#include "pmlist.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

/* */
/* Read Mapping Table */
/* */
int stringToParam(const char *paramstr);

static char *kv_get_a_val(list_t *kvl, const char *key, const char *replacer)
{
  keyValues_t *kv = kvlist_search(kvl, key);
  if ( kv )
    return kv->values[0];
  if ( replacer )
    return (char *)replacer;
  else
    return NULL;
}

list_t *maptab_search_miptab(list_t *pmlist, const char *cmorname, const char *miptab, const char *key)
{
  if ( pmlist && cmorname && miptab )
    {
      listNode_t *node = pmlist->head;
      list_t *listlatest = NULL;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvlist = *(list_t **)node->data;
              keyValues_t *kvcn = kvlist_search(kvlist, key);
              if ( kvcn && kvcn->nvalues > 0 && *(kvcn->values[0]) == *cmorname && strcmp(kvcn->values[0], cmorname) == 0 )
                {
                  keyValues_t *kvmt = kvlist_search(kvlist, "mip_table");
                  if ( ( kvmt && kvmt->nvalues > 0 && *(kvmt->values[0]) == *miptab && strcmp(kvmt->values[0], miptab) == 0 ) || !kvmt )
                    return kvlist;
                  else
                    listlatest = kvlist;
                }
            }
          node = node->next;
        }
      if ( listlatest )
        {
          if ( cdoVerbose )
            cdoPrint("No attribute 'mip_table' found in mapping table line for cmorname '%s'.\n          The latest line of the mapping table is used.", cmorname);
          return listlatest;
        }
    }

  return NULL;
}

static
char *readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  int ichar;
  size_t ipos = 0;
  while ( *buffersize )
    {
      ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if ( ichar == '\r' )
        {
          if ( *buffersize )
            {
              ichar = *buffer;
              if ( ichar == '\n' )
                {
                  (*buffersize)--;
                  buffer++;
                }
            }
          break;
        }
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %ld)!\n", len);
          break;
        }
    }
  line[ipos] = 0;
  if ( *buffersize == 0 && ipos == 0 ) buffer = NULL;
  return buffer;
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return pline;
}

static void *handleError(list_t *kvl, int errnum, char *argument)
{
  char *filename = NULL;
  if ( kvl )
    filename = kv_get_a_val(kvl, "workfile4err", NULL);
  else 
    cdoAbort("In parsing the command line:\n          More than 100 values for a key are not supported.");
  if ( !filename )
    cdoAbort("In parsing file:\n          Cannot resolve recent working file.");
  switch ( errnum )
    {
    case ( 1 ): cdoAbort("In parsing file '%s':\n          Unexpected blank in line:\n          '%s'\n          Check syntax.", filename, argument);
    case ( 2 ): cdoAbort("In parsing file '%s':\n          Unexpected separator sign ',' in a key of line:\n          '%s'\n          Check syntax.", filename, argument);
    case ( 3 ): cdoAbort("In parsing file '%s':\n          More than 100 values for a key are not supported.", filename);
    case ( 4 ): cdoAbort("In parsing file '%s':\n          No values found for a keyword in line:\n          '%s'.\n          Check syntax.", filename, argument);
    case ( 5 ): cdoAbort("In parsing file '%s':\n          A value for a keyword begins with ',' in line:\n          '%s'.\n          Check syntax.", filename, argument);
    case ( 6 ): cdoAbort("In parsing file '%s':\n          A Value for a keyword has a start quote sign but no end quote sign in line:\n          '%s'.\n          Check syntax.",filename,  argument);
    case ( 7 ): cdoAbort("In parsing file '%s':\n          Unexpected separator sign '=' or ':' is found in values in line:\n          '%s'.", filename, argument);
    case ( 8 ): cdoAbort("In parsing file '%s':\n          Connected lines for one keyvalue contain more than the allowed 4096 characters.",filename);
    case ( 9 ): cdoAbort("In parsing file '%s':\n          A ',' is found at end of a line without information in next line.",filename);
    case ( 10): cdoAbort("In parsing file '%s':\n          A ',' is found at end of file.",filename);
    }
}

static
char *getElementName(char *pline, char *name, int *errh)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  size_t pos = 0;
  while ( pos < len && *(pline+pos) != '=' )
    {
      if ( isspace((int) *(pline+pos)) )
        { *errh=1; return pline; }
      if ( *(pline+pos) == ',' )
        { *errh=2; return pline; }
      if ( *(pline+pos) == ':' )
        cdoWarning("In parsing file:\n          Separator sign ':' is not supported. Use '=' instead.\n          ...'%s'...", pline);
      name[pos] = tolower(*(pline+pos));
      pos++;
    }
  name[pos] = 0;

  pline += pos;
  return pline;
}

static int copy_value(char *value, char **values, int *nvalues)
{
  if ( *nvalues > 100 )
    return 3;
  values[*nvalues] = strdup(value);
  (*nvalues)++;
  values[*nvalues] = NULL;
  return 0;
}


static void free_array(char **tofree)
{
  int i = 0;
  while ( tofree[i] )
    {
      free(tofree[i]);
      i++;
    }
  Free(tofree);
}

static void quote_replace(char **values, int nvalues, int i)
{
  char *source = values[nvalues];
  char *useful;
  *source++;
  source[i-2] = 0;
  useful = strdup(source);
  free(values[nvalues]);
  values[nvalues] = strdup(useful);
  free(useful);
}

static
char *getElementValues(char *pline, char **values, int *nvalues, char *keyword, int *errh)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  while ( isspace((int) *(pline+(int)len)) ) len--;
  *(pline+len) = 0;
  if ( (int)len == 0 )
    {*errh=4; return keyword;}
  *nvalues = 0;
  int i = 0;
  while ( i < len && len )
    {
      if ( *(pline+i) == ',')
        {
          if ( i == 0 )
            {
              *errh=5;
              return pline;
            }
          *errh = copy_value(pline, values, nvalues);
          if ( *errh > 0 )
            return pline;
          if ( *(values[*nvalues-1]+i-1) == '"' || *(values[*nvalues-1]+i-1) == '\'' )
            quote_replace(values, *nvalues-1,i);
          else
            *(values[*nvalues-1]+i) = 0;

          i++;
          pline+=i;
          len-=i;
          i=0;
        }
      else if ( *(pline+i) == '"' )
        {
          i++;
          while ( *(pline+i) != '"' )
            {
              i++;
              if ( *(pline+i) == 0 )
                {*errh=6; return pline;}
            }
          i++;
        }
      else if ( isspace((int) *(pline+i)) )
        break;          
      else if ( *(pline+i) == '=' || *(pline+i) == ':' )
        {
          *errh=7;
          return pline;
        }
      else
        i++;
    }
  *errh = copy_value(pline, values, nvalues);
  if ( *errh > 0 )
    return pline;
  if ( *(values[*nvalues-1]+i-1) == '"' )
    quote_replace(values, *nvalues-1, i);
  else
    *(values[*nvalues-1]+i) = 0;
  pline+=i;
  return pline;
}

static int parse_line_to_list(list_t *list, char *pline, const char *kvlname, int checkpml, int lowprior)
{
  int errh = 0;
  char name[256];
  int i = 0, nvalues;
  list_t *kvl = NULL;
  if ( checkpml )
    {
      kvl = list_new(sizeof(keyValues_t *), free_keyval, kvlname);
      list_append(list, &kvl);
    }
  while ( *pline != 0 )
    {
      char **values = (char **) Malloc( 100 * sizeof(char *) );
      
      pline = getElementName(pline, name, &errh);
      if ( strcmp(name, "conventions") == 0 )
        strcpy(name, "Conventions");
      if ( strcmp(name, "cordex_domain") == 0 )
        strcpy(name, "CORDEX_domain");
      if ( errh > 0 )
        return errh;
      if ( *pline == 0 )
        {
          cdoPrint("In parsing a line of a file:\n          Could not find values for key: '%s'. Use correct 'key=value' syntax.", name);
          break;
        }
      pline = skipSeparator(pline);
      if ( *pline == 0 )
        {
          cdoPrint("In parsing a line of a file:\n          Could not find values for: '%s'. Use correct 'key=value' syntax.", name);
          break;
        }
      pline = getElementValues(pline, values, &nvalues, name, &errh);
      if ( errh > 0 )
        return errh;

      if ( checkpml )
        kvlist_append(kvl, name, (const char **)values, nvalues);
      else
        {
          if ( lowprior )
            {
              keyValues_t *already = kvlist_search(list, name);
              if ( already )
                continue;
            }
          kvlist_append(list, name, (const char **)values, nvalues);
        }
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '/' )
        *pline = 0;
      free_array(values);
    }
  return 0;
}

static void remove_space_and_comms(char **pline, char *line)
{
  while ( isspace((int) *(*pline)) ) (*pline)++;
  char *tester = *pline;
  int i = 0;
  while ( *tester != 0 )
    {
      if ( *tester == '#' || *tester == '!' )
        {
          line[i] = '\0';
          break;
        }
      i++; tester++;
    }
}

static int add_lines_tester(char *line)
{
  char *tester = line;
  while ( *tester != 0 ) tester++;
  tester--;
  while ( isspace((int) *tester) ) tester--;
  if ( *tester == ',' )
    return 1;
  return 0;
}

static int add_lines(char *line, char **buffer, size_t *buffersize)
{
  int len = strlen(line);
  char nextline[4096];
  if ( (*buffer = readLineFromBuffer(*buffer, buffersize, nextline, sizeof(nextline))) )
    {
      char *nexttester = nextline;
      remove_space_and_comms(&nexttester, nextline);
      if ( *nexttester != '\0' && *nexttester != '&' )
        {
          if ( strlen(nexttester) + len > 4096 )
            return 8;
          strcat(line, nextline);
        }
      else
        return 9;
    }
  else
    return 10;
}  

void parse_buffer_to_list(list_t *list, size_t buffersize, char *buffer, int checkpml, int lowprior, list_t *kvl)
{
  char line[4096];
  char name[256];
  char *pline;
  const char *listkeys[] = {"axis_entry:", "variable_entry:", "&parameter", NULL};
  int linenumber = 0;
  int listtype = 0;
  int errh = 0;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      remove_space_and_comms(&pline, line);
      if ( *pline == '\0' ) continue;
      while ( add_lines_tester(line) )
        {
          errh = add_lines(line, &buffer, &buffersize);
          if ( errh > 0 )
            handleError(kvl, errh, NULL);
        }
      //  len = (int) strlen(pline);
      if ( listtype == 0 && *pline == '&' ) listtype = 1;
/* MAXNVALUES*/
      int i = 0;
      while ( listkeys[i] )
        {
          if ( strlen(pline) > strlen(listkeys[i]) )
            if ( strncmp(pline, listkeys[i], strlen(listkeys[i])) == 0 )
              {
	        pline += strlen(listkeys[i]);
 	        listtype = 2;
                break;
	      }
          i++;
        }
      if ( listtype )
        errh = parse_line_to_list(list, pline, listkeys[i], checkpml, lowprior);
      else
        errh = parse_line_to_list(list, pline, "keyvals", checkpml, lowprior);
      if ( errh > 0 )
        handleError(kvl, errh, pline);
    }
}


static void kv_insert_a_val(list_t *kvl, const char *key, char *value, int replace)
{
  if ( key )
    {
      keyValues_t *kv = kvlist_search(kvl, key);
      if ( !kv )
        {
          const char *apvalue[] = {value};
          kvlist_append(kvl, key, apvalue, 1);
        }
      else if ( replace )
        {
          Free(kv->values[0]);
          kv->values[0] = strdup((const char *) value);
        }
    }
}

static char **kv_get_vals(list_t *kvl, const char *key, int *numvals)
{
  keyValues_t *kv = kvlist_search(kvl, key);
  if ( kv )
    { *numvals = kv->nvalues; return kv->values; }
  return NULL;
}

list_t *cdo_parse_cmor_file(const char *filename, list_t *kvl)
{
  assert(filename != NULL);

  size_t filesize = fileSize(filename);
  if ( filesize == 0 )
    {
      fprintf(stderr, "In reading the mapping table:\n          Empty table file: %s.", filename);
      return NULL;
    }

  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "In reading the mapping table:\n          Open failed on %s: %s.", filename, strerror(errno));
      return NULL;
    }

  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "In reading the mapping table:\n          Read failed on %s.", filename);
      return NULL;
    }
 
  list_t *pml = list_new(sizeof(list_t *), free_kvlist, filename);

/*  if ( buffer[0] == '{' )
    parse_json_buffer_to_pml(pml, filesize, buffer);
  else */
  parse_buffer_to_list(pml, filesize, buffer, 1, 0, kvl);

  if ( buffer ) Free(buffer);

  return pml;
}

static void get_stringcode(int vlistID, int varID, char *varcodestring)
{
  int varcode;
  varcode = vlistInqVarCode(vlistID, varID);
  sprintf(varcodestring, "%03d", varcode);
}

static void get_ifilevalue(char *ifilevalue, const char *key, int vlistID, int varID)
{
  if ( strcmp(key, "name") == 0 )
    vlistInqVarName(vlistID, varID, ifilevalue);
  else if ( strcmp(key, "code") == 0 )
    {
      char varcodestring[CDI_MAX_NAME];
      get_stringcode(vlistID, varID, varcodestring);
      strcpy(ifilevalue, varcodestring);
    }
}

static int getVarIDToMap(int vlistID, int nvars, const char *key, const char *value)
{
  for ( int varID = 0; varID < nvars; varID++ )
    {
      char ifilevalue[CDI_MAX_NAME];
      get_ifilevalue(ifilevalue, key, vlistID, varID);
      if ( strcmp(ifilevalue, value) == 0 )
        return varID;
    }
  return CDI_UNDEFID;
}


static const char *check_short_key(char *key)
{
  const char *short_keys[]={"cn", "n", "c", "u", "cm", "vc", "p", "szc", "i", "ca", "za", "gi", "rtu", "mt", "om", "ms", "dr", "d", "lc", "dj", NULL};
  const char *long_keys[]={"cmor_name", "name", "code", "units", "cell_methods", "variable_comment", "positive", "scalar_z_coordinate", "info", "character_axis", "z_axis",  "grid_info", "required_time_units", "mapping_table", "output_mode", "max_size", "drs_root", "drs", "last_chunk", "dataset_json", NULL};

  for ( int i = 0; short_keys[i]; i++ )
    if ( strcmp(key, short_keys[i]) == 0 || strcmp(key, long_keys[i]) == 0 )
      return short_keys[i];
/*  if ( strcmp(key, "cmor_name") == 0 ) short_key = strdup("cn");
  else if ( strcmp(key, "name") == 0 ) short_key = strdup("n");
  else if ( strcmp(key, "code") == 0 ) short_key = strdup("c");
  else if ( strcmp(key, "units") == 0 ) short_key = strdup("u");
  else if ( strcmp(key, "cell_methods") == 0 ) short_key = strdup("cm");
  else if ( strcmp(key, "comment") == 0 ) short_key = strdup("k");
  else if ( strcmp(key, "positive") == 0 ) short_key = strdup("p");
  else if ( strcmp(key, "scalar_z_coordinate") == 0 ) short_key = strdup("szc");
  else if ( strcmp(key, "info") == 0 ) short_key = strdup("i");
  else if ( strcmp(key, "character_axis") == 0 ) short_key = strdup("ca");
  else if ( strcmp(key, "z_axis") == 0 ) short_key = strdup("za");
  else if ( strcmp(key, "grid_info") == 0 ) short_key = strdup("gi");
  else if ( strcmp(key, "required_time_units") == 0 ) short_key = strdup("rtu");
  else if ( strcmp(key, "mapping_table") == 0 ) short_key = strdup("mt");
  else if ( strcmp(key, "calendar") == 0 ) short_key = strdup("l");
  else if ( strcmp(key, "output_mode") == 0 ) short_key = strdup("om");
  else if ( strcmp(key, "max_size") == 0 ) short_key = strdup("ms");
  else if ( strcmp(key, "drs_root") == 0 ) short_key = strdup("dr");
  else if ( strcmp(key, "drs") == 0 ) short_key = strdup("d");
  else if ( strcmp(key, "last_chunk") == 0 ) short_key = strdup("lc"); 
  else if ( strcmp(key, "dataset_json") == 0 ) short_key = strdup("dj"); */
  return NULL;
}

static void map_it(list_t *kvl, int vlistID, int varID, char *var2map)
{
  for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      const char *key = ( check_short_key((char *)kv->key ) ) ? check_short_key((char *)kv->key) : NULL;
      if ( !key )
        {
          if ( strcmp(kv->key, "mip_table") != 0 )
            cdoWarning("In variable mapping:\n           you try to assign '%s' to variable '%s'.\n          This mapping table keyword is skipped. Check allowed mapping table keywords.", kv->key, var2map);
          continue;
         }
      const char *value = kv->values[0];
/*      printf("'%s' = '%s'\n", key, value); */
      if ( !value ) continue;
/* Not necessary because cmor_name is what we use for renaming :
      else if ( STR_IS_EQ(key, "name")          ) vlistDefVarName(vlistID, varID, parameter2word(value));
*/
      else if ( STR_IS_EQ(key, "cn") )
        {
          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, varID, name);
          if ( name[0] != 0 )
            cdiDefAttTxt(vlistID, varID, "original_name", (int) strlen(name), parameter2word(name));
          vlistDefVarName(vlistID, varID, parameter2word(value));
         }
      else if ( STR_IS_EQ(key, "u")        ) vlistDefVarUnits(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "cm")  ) cdiDefAttTxt(vlistID, varID, "cell_methods", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "ca")  ) cdiDefAttTxt(vlistID, varID, "character_axis", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "za")  ) cdiDefAttTxt(vlistID, varID, "z_axis", (int) strlen(value), value);
/*      else if ( STR_IS_EQ(key, "factor")        ) {}
      else if ( STR_IS_EQ(key, "delete")        ) {}
      else if ( STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "code")          ) {}
              // else if ( STR_IS_EQ(key, "code")          ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              // else if ( STR_IS_EQ(key, "out_code")      ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
*/
      else if ( STR_IS_EQ(key, "vc")       ) cdiDefAttTxt(vlistID, varID, "comment", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "p")       )
        {
          if ( !isspace(value[0]) )
            cdiDefAttTxt(vlistID, varID, "positive", (int) strlen(value), value);
        }
/*      else if ( STR_IS_EQ(key, "cell_measures") ) cdiDefAttTxt(vlistID, varID, "cell_measures", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "convert")       ) {} 
      else if ( STR_IS_EQ(key, "missval")       )   {}  
      else if ( STR_IS_EQ(key, "valid_min")     ){}
      else if ( STR_IS_EQ(key, "valid_max")     ){}
      else if ( STR_IS_EQ(key, "ok_min_mean_abs") ) {}
      else if ( STR_IS_EQ(key, "ok_max_mean_abs") ){}
      else if ( STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type") ) {} */
      else
        {
          if ( cdoVerbose ) cdoPrint("In applying the mapping table:\n          Key: '%s' is ignored.", key);
        }
    }
}


static int change_name_via_name(int vlistID, char *map_name, char *cmor_name)
{
  char name[CDI_MAX_NAME];
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      vlistInqVarName(vlistID, varID, name);
      if ( strcmp(name, map_name) == 0 )
        {
          vlistDefVarName(vlistID, varID, parameter2word((const char *) cmor_name));
          return 1;
        }
    }
  return 0;
}

static int change_name_via_code(int vlistID, char *map_code, char *cmor_name)
{
  int code;
  char codestring[4];
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      code = vlistInqVarCode(vlistID, varID);
      sprintf(codestring, "%03d", code);
      if ( strcmp(codestring, map_code) == 0 )
        {
          vlistDefVarName(vlistID, varID, parameter2word((const char *) cmor_name));
          return 1;
        }
    }
  return 0;
}


static int maptab_via_key(list_t *pml, int vlistID, int varID, int nventry, const char **ventry, const char *key, char *miptabfreq)
{
  char ifilevalue[CDI_MAX_NAME];
  if ( strcmp(key, "cmor_name") != 0 )
    get_ifilevalue(ifilevalue, key, vlistID, varID);
  else
    get_ifilevalue(ifilevalue, "name", vlistID, varID);

  if ( ifilevalue[0] )
    {
      list_t *kvl = maptab_search_miptab(pml, ifilevalue, miptabfreq, key);
      if ( kvl )
        {
          if ( cdoVerbose )
            cdoPrint("Start to map via '%s'.", key);
          map_it(kvl, vlistID, varID, ifilevalue);
          return 1;
        }
      if ( cdoVerbose )
        cdoPrint("In variable mapping:\n          Variable named '%s' with varID '%d' could not be mapped via '%s' because no corresponding key '%s' was found in mapping table file.", ifilevalue, varID, key, key);
      return 0;
    }
  else
    {
      if ( cdoVerbose )
        cdoPrint("In variable mapping:\n          Variable with varID '%d' could not be mapped via '%s' because it does not possess a '%s' in infile.", varID, key, key);
      return 0;
    }
}

static int maptab_via_cn_and_key(list_t *kvl_oname, int vlistID, int nvars, const char *key)
{
  keyValues_t *kv = kvlist_search(kvl_oname, key);
  keyValues_t *kvcn = kvlist_search(kvl_oname, "cmor_name");
  if ( kv )
    {
      int varID = ( strcmp(key, "cmor_name") == 0 ) ? getVarIDToMap(vlistID, nvars, "name", kv->values[0]) : getVarIDToMap(vlistID, nvars, key, kv->values[0]);
      if ( varID != CDI_UNDEFID )
        {
          if ( cdoVerbose )
            cdoPrint("Started mapping of variable via '%s'.", key);
          map_it(kvl_oname, vlistID, varID, kv->values[0]);
          return 1;
        }
      cdoPrint("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped via key '%s' because no infile variable '%s' equals '%s'.", kvcn->values[0], key, key, kv->values[0]);
    }
  else
    if ( cdoVerbose )
      cdoPrint("In variable mapping:\n          Variable '%s' configured via cmor_name\n          could not be mapped via key '%s' because it possesses no corresponding key '%s' in mapping file.", kvcn->values[0], key, key);
  return 0;
}

static void maptab_via_cmd(list_t *pml, const char *origValue, int vlistID, int nvars, const char *key, char *cmorName, char *miptabfreq)
{
  int varIDToMap = getVarIDToMap(vlistID, nvars, key, origValue);
  if ( varIDToMap == CDI_UNDEFID )
    cdoAbort("In variable mapping:\n          Variable with '%s': '%s' configured via cmdline could not be found in infile '%s'.", key, origValue, cdoStreamName(0)->args);
  list_t *kvl_maptab = maptab_search_miptab(pml, cmorName, miptabfreq, "cmor_name");
  if ( !kvl_maptab )
    {
      cdoWarning("In variable mapping:\n          The registered cmor_name '%s' via cmdline could not be found in mapping table.\n          No mapping table is applied.", cmorName);
      vlistDefVarName(vlistID, varIDToMap, parameter2word((const char *) cmorName));
    }
  else
    {
      if ( cdoVerbose )
        cdoPrint("Started mapping of variable via '%s'.", key);
      map_it(kvl_maptab, vlistID, varIDToMap, cmorName);
    }
}

static void maptab_via_cn(list_t *pml, char **request, int vlistID, int nvars, int numvals, char *miptabfreq, int filetype)
{
  for ( int j = 0; j<numvals; j++)
    {
      list_t *kvl_oname = maptab_search_miptab(pml, request[j], miptabfreq, "cmor_name");
      if ( kvl_oname )
        {
          if ( filetype == FILETYPE_GRB || filetype ==  FILETYPE_GRB2 )
            {
              if ( maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "code") )
                {
                  if ( cdoVerbose )
                    cdoPrint("Successfully mapped variable via code to cmor_name '%s'.", request[j]); 
                  continue;
                }
            }
          if ( maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "name") )
            {
              if ( cdoVerbose )
                cdoPrint("Successfully mapped variable via name to cmor_name: '%s'.", request[j]); 
              continue;
            }
          else if ( ( filetype != FILETYPE_GRB && filetype != FILETYPE_GRB2 ) && maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "code") )
            {
              if ( cdoVerbose )
                cdoPrint("Successfully mapped variable via code to cmor_name '%s'.", request[j]); 
              continue;
            }
          else
            {
              if ( cdoVerbose )
                cdoPrint("In variable mapping:\n          Try to use cmor_name for selecting the infile variable.", request[j]);
              if ( maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "cmor_name") )
                {
                  cdoPrint("Successfully mapped variable via cmor_name to cmor_name '%s'.", request[j]); 
                  continue;
                }
              cdoWarning("In variable mapping:\n          Mapping table line of cmor_name '%s' could neither be mapped via 'name', 'code' nor 'cmor_name'.\n          No mapping for cmor_name: '%s'.", request[j],request[j]);
              continue;
            }
        }
      else
        {
          cdoWarning("In variable mapping:\n          Requested cmor_name: '%s' is not found in mapping table.\n          No mapping for cmor_name: '%s'", request[j],request[j]);
          continue;
        }
    }
}

/* */
/*... until here */
/* */

struct mapping
{
  int help_var;
  int cdi_varID;
  int cmor_varID;
  int zfactor_id;
  int charvars;
  char datatype;
  void *data;
};

static struct mapping *construct_var_mapping(int streamID)
{
  int nvars_max = vlistNvars(pstreamInqVlist(streamID));
  struct mapping *vars =
    (struct mapping *) Malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  vars[0].data = NULL;
  vars[0].charvars = 0;
  return vars;
}

static void destruct_var_mapping(struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    Free(vars[i].data);
  Free(vars);
}

static struct mapping *map_var(int cdi_varID, struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    if ( cdi_varID == vars[i].cdi_varID )
      return &vars[i];
  return NULL;
}

static struct mapping *new_var_mapping(struct mapping vars[])
{
  int i;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  vars[i + 1].cdi_varID = CDI_UNDEFID;
  vars[i + 1].data = NULL;
  vars[i + 1].charvars = 0;
  return &vars[i];
}

static int *new_axis_id(int *axis_ids)
{
  int i;
  for ( i = 0; axis_ids[i] != CMOR_UNDEFID; i++ );
  axis_ids[i + 1] = CMOR_UNDEFID;
  return &axis_ids[i];
}

static int count_axis_ids(int *axis_ids)
{
  int i;
  for ( i = 0; axis_ids[i] != CMOR_UNDEFID; i++ );
  return i;
}

static void addcharvar(keyValues_t *charvars, int vlistID, const char *key, struct mapping vars[])
{
  if ( cdoVerbose )
    cdoPrint("Start to merge variables to one character coordinate.");
  int varIDs[charvars->nvalues];
  int nvars = vlistNvars(vlistID);
  for ( int i = 0; i < charvars->nvalues; i++)
    {
      varIDs[i] = getVarIDToMap(vlistID, nvars, key, charvars->values[i]);
      if ( varIDs[i] == CDI_UNDEFID )
        cdoAbort("In merging variables to a variable with a character coordinate:\n          Could not find '%s' in infile '%s' to build a variable with character coordinate.", charvars->values[i], cdoStreamName(0)->args);
    }

  int gridID = vlistInqVarGrid(vlistID, varIDs[0]);
  int subgridID;
  int zaxisID = vlistInqVarZaxis(vlistID, varIDs[0]);
  int subzaxisID;
  int ntsteps = vlistNtsteps(vlistID);

  if ( cdoStreamName(0)->args[0] == '-' )
    cdoAbort("No variables can be merged to one character axis since you piped several cdo operators.");

  int streamID2 = pstreamOpenRead(cdoStreamName(0));
  if ( ntsteps == -1 )
    {
      ntsteps = 0;
      int dummy;
      while ( ( dummy = pstreamInqTimestep(streamID2, ntsteps++) ) );
    }

  int axissize[3];
  double *xvals, *yvals, *zvals, *subsvals;

  subsvals = (double *) Malloc(charvars->nvalues * sizeof(double));
  for ( int i = 0; i < charvars->nvalues; i++ )
    subsvals[i] = i+1;

  axissize[0] = gridInqXsize(gridID);
  axissize[1] = gridInqYsize(gridID);
  axissize[2] = zaxisInqSize(zaxisID);


  if ( axissize[0] != 1 && axissize[1] != 1 && axissize[2] != 1 )
    cdoAbort("In merging variables to a variable with a character coordinate:\n          No axis found to merge. One axis may not be allocated with more than one value.");

  int oldgridsize = axissize[0] * axissize[1];
  double *buffer_old = (double *)Malloc(oldgridsize * sizeof(double));

  for ( int i = 1; i < charvars->nvalues; i++)
    {
      gridID = vlistInqVarGrid(vlistID, varIDs[i]);  
      zaxisID = vlistInqVarZaxis(vlistID, varIDs[i]);
      if ( axissize[0] != gridInqXsize(gridID) )
        cdoAbort("In merging variables to a variable with a character coordinate:\n          Size of x-axis: '%d' of variable '%s'\n          differ from x-axis size of variable '%s': '%d'.", gridInqXsize(gridID), charvars->values[i], charvars->values[0], axissize[0]);
      if ( axissize[1] != gridInqYsize(gridID) )
        cdoAbort("In merging variables to a variable with a character coordinate:\n          Size of y-axis: '%d' of variable '%s'\n          differ from y-axis size of variable '%s': '%d'.", gridInqYsize(gridID), charvars->values[i], charvars->values[0], axissize[1]);
      if ( axissize[2] != zaxisInqSize(zaxisID) )
        cdoAbort("In merging variables to a variable with a character coordinate:\n          Size of z-axis: '%d' of variable '%s'\n          differ from z-axis size of variable '%s': '%d'.", zaxisInqSize(zaxisID), charvars->values[i], charvars->values[0], axissize[2]);
    }

  if ( axissize[0] == 1 )
    {
      xvals = subsvals;
      yvals = (double *) Malloc(axissize[1] * sizeof(double));
      zvals = (double *) Malloc(axissize[2] * sizeof(double));
      gridInqYvals(gridID, yvals);
      zaxisInqLevels(zaxisID, zvals); 
      axissize[0] = charvars->nvalues;
    }
  else if ( axissize[1] == 1 )
    {
      xvals = (double *) Malloc(axissize[0] * sizeof(double));
      yvals = subsvals;
      zvals = (double *) Malloc(axissize[2] * sizeof(double));
      gridInqXvals(gridID, xvals);
      zaxisInqLevels(zaxisID, zvals); 
      axissize[1] = charvars->nvalues;
    }
  else if ( axissize[2] == 1 )
    {
      xvals = (double *) Malloc(axissize[0] * sizeof(double));
      yvals = (double *) Malloc(axissize[1] * sizeof(double));
      zvals = subsvals;
      gridInqXvals(gridID, xvals);
      gridInqYvals(gridID, yvals);
      axissize[2] = charvars->nvalues;
    }

  subgridID = gridCreate(GRID_GENERIC, axissize[0]*axissize[1]);
  subzaxisID = zaxisCreate(zaxisInqType(zaxisID), axissize[2]);

  gridDefXsize(subgridID, axissize[0]); 
  gridDefYsize(subgridID, axissize[1]); 
  gridDefXvals(subgridID, xvals); 
  gridDefYvals(subgridID, yvals); 
  zaxisDefLevels(subzaxisID, zvals); 

  struct mapping *var = new_var_mapping(vars);
  var->cdi_varID = vlistDefVar(vlistID, subgridID, subzaxisID,  TIME_VARYING);
  vlistDefVarName(vlistID, getVarIDToMap(vlistID, nvars+1, key, charvars->values[0]), "ChangedForMap");
  vlistDefVarName(vlistID, var->cdi_varID, charvars->values[0]);
  vlistDefVarDatatype(vlistID, var->cdi_varID,  DATATYPE_FLT64);
  vlistDefVarMissval(vlistID, var->cdi_varID, vlistInqVarMissval(vlistID, varIDs[0]));
  var->datatype = 'd';
  var->data = Malloc(ntsteps*axissize[0]*axissize[1]*axissize[2]*sizeof(double));

  int testzaxisID = vlistInqVarZaxis(vlistID, var->cdi_varID);

  int tsID = 0, nrecs = 0;

  while ( (nrecs = pstreamInqTimestep(streamID2, tsID)) )
    {
      while ( nrecs-- )
        {
          int varIDrw, levelIDrw, nmiss;
          pstreamInqRecord(streamID2, &varIDrw, &levelIDrw);
          for ( int i = 0; i < charvars->nvalues; i++ )
            if ( varIDrw == varIDs[i] )
              {
                pstreamReadRecord(streamID2, buffer_old, &nmiss);
                int newIndex;
                for ( int j = 0; j < oldgridsize; j++ )
                  {
/* (lev x lat, basin ) 
            newIndex = j * levdim + levelID; */
                    if ( oldgridsize == axissize[0]*axissize[1] )
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+i*axissize[0]*axissize[1] + j;
                    else if ( axissize[0] == charvars->nvalues )
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+i*axissize[1]*axissize[2] + j*axissize[2] + levelIDrw;
                    else
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+levelIDrw*axissize[0]*axissize[1] + i*axissize[0]*axissize[1]/oldgridsize + j;
                    ((double *)var->data)[newIndex] = (double) buffer_old[j];
                  }
              }
        }
      tsID++;
    }
  var->charvars = 1;

  pstreamClose(streamID2);
  Free(buffer_old);

  if ( cdoVerbose )
    cdoPrint("Successfully merged variables into one character axis. The final variable is called '%s' and has the ID: '%d'", charvars->values[0], var->cdi_varID);
}

static char *trim(char *s)
{
  if (s == NULL) return s;
  while ( *s != '\0' && (isspace(*s) || *s == '"') )
    s++;
  int n = strlen(s);
  while ( n > 0 && (isspace(s[n - 1]) || s[n - 1] == '"') )
    n--;
  s[n] = '\0';
  return s;
}

static int file_exist(const char *tfilename, int force, const char *fileart)
{
  assert(tfilename != NULL);
  size_t filesize = fileSize(tfilename);
  if ( filesize == 0 && force)
    cdoAbort("In checking wether a %s file exist:\n          Empty file: '%s'.", fileart, tfilename);
  else if ( filesize == 0 && !force )
    {
      cdoPrint("In checking wether a %s file exist:\n          Empty file: '%s'.", fileart, tfilename);
      return 0;
    }
  if ( strstr(tfilename, ".nc") || strstr(tfilename, ".grb") )
    return 1;
  FILE *fp = fopen(tfilename, "r");
  if ( fp == NULL && force )
    cdoAbort("In checking wether a %s file exist:\n          Open failed on: '%s'.", fileart, tfilename);
  else if ( fp == NULL && !force )
    {
      cdoPrint("In checking wether a %s file exist::\n          Open failed on: '%s'.", fileart, tfilename);
      return 0;
    }

  fclose(fp);
  return 1;
}
  
static int parse_kv_file(list_t *kvl, const char *filename)
{
  file_exist(filename, 1, "Configuration");

  FILE *fp = fopen(filename, "r");
  size_t filesize = fileSize(filename);
  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);
  fclose(fp);

  kv_insert_a_val(kvl, "workfile4err", (char *)filename, 1);
  parse_buffer_to_list(kvl, filesize, buffer, 0, 1, kvl);

  if ( buffer ) Free(buffer);
  return 1;
}

static void check_compare_set(char **finalset, char *attribute, const char *attname, const char *defaultstr)
{
  if ( !(*finalset) )
    {
      if ( !attribute )
        {
          if ( defaultstr )
            *finalset = strdup(defaultstr);
          else
            cdoAbort("In comparison of configuration attribute and infile attribute:\n          Required value for attribute '%s' is neither found in input file nor in the configuration.", attname);
        }
      else 
        *finalset = strdup(attribute);
    }
  else if ( attribute )
    {
      if ( strcmp(attribute, *finalset) != 0 )
        {
          cdoPrint("In comparison of configuration attribute and infile attribute:\n          '%s' of variable in input file: '%s' does not agree with configuration attribute %s: '%s'.\n          Cmor libary is called with attribute unit '%s'.", attname, *finalset, attname, attribute, attribute);
          Free(*finalset);
          *finalset = strdup(attribute);
        }
    }
}

static void add_globalhybrids(list_t *kvl)
{
  const char *longAtt[] = {"required_time_units", "grid_info", "mapping_table", NULL};
  const char *shortAtt[] = {"rtu", "gi", "mt", NULL};

  int i = 0;
  while ( longAtt[i] != NULL )
    {
      keyValues_t *kv_latt = kvlist_search(kvl, longAtt[i]);      
      keyValues_t *kv_satt = kvlist_search(kvl, shortAtt[i]);      
      if ( kv_latt && !kv_satt )
        kv_insert_a_val(kvl, shortAtt[i], kv_latt->values[0], 1);
      else if ( !kv_latt && kv_satt )
        kv_insert_a_val(kvl, longAtt[i], kv_satt->values[0], 1);      
      i++;
    }
}

static int check_attarray(list_t *kvl, const char **reqAtt)
{
  int i = 0;
  while ( reqAtt[i] != NULL )
    {
      keyValues_t *kv_reqatt = kvlist_search(kvl, reqAtt[i]);
      if ( !kv_reqatt || strcmp(kv_reqatt->values[0], "notSet") == 0 )
        return i;
      if ( cdoVerbose )
        cdoPrint("Attribute '%s' is '%s'. ", reqAtt[i], kv_reqatt->values[0]);
      i++;
    }
  return -1;
}

static void attErr(const char **reqAtt, int errnum)
{
  char errStr[CMOR_MAX_STRING];
  int i = 1;

  sprintf(errStr, "Attribute '%s' is required. Either it is missing or notSet.\n          Make sure that you have configured all following attributes:\n          %s", reqAtt[errnum], reqAtt[0]);
  while ( reqAtt[i] )
    {
      sprintf(errStr, "%s, %s", errStr, reqAtt[i]);
      i++;
    }
  cdoAbort(errStr);
}

static int check_attr(list_t *kvl, char *project_id)
{
/* Project id moved to main void fct */
  const char *reqAtt[] = {"institution", "contact", "model_id", "source",
            "experiment_id", "required_time_units", NULL};
  const char *reqAttCMIP5[] = {"institute_id", "product", "member", NULL};
  const char *reqAttCMIP5CMOR3[] = {"modeling_realm", NULL};
  const char *reqAttCORDEX[] = {"institute_id", "product", "member", "CORDEX_domain", "driving_model_id", NULL};
  const char *reqAttCMIP6CMOR3[] = {"outpath", "output_path_template", "output_file_template", "tracking_prefix", NULL};
  const char *reqAttCMIP6[] = {"Conventions", "activity_id", "experiment", "further_info_url", "grid", "grid_label", "institution_id", "license", "mip_era", "nominal_resolution", "product", "source_id", "source_type", "variant_label", NULL};
  const char *expdepAttCMIP6[] = {"parent_experiment_id", "parent_activity_id", "parent_mip_era", "parent_source_id", "parent_variant_label", "parent_time_units", "sub_experiment", "sub_experiment_id", "branch_method", NULL};
/* In all Projects needed Attributes are tested first */

  int errnum = 0;
  if ( ( errnum = check_attarray(kvl, reqAtt) ) != -1 )
    attErr(reqAtt, errnum);
/* Set default attributes */
  keyValues_t *kv = kvlist_search(kvl, "references");

  if ( !kv || strcmp(kv->key, "notSet") == 0 )
    {
      keyValues_t *kv_model_id = kvlist_search(kvl, "model_id");
      char *references = (char *) Malloc(strlen(kv_model_id->values[0]) + 28);
      strcpy(references, "No references available for ");
      strcat(references, kv_model_id->values[0]);
      cdoPrint("Attribute 'references' is set to '%s' ", references);
      kv_insert_a_val(kvl, "references", references, 1);
      Free(references);
    }

/* Special check for CMIP or CORDEX projects */
  if ( strcmp(project_id, "CMIP5") == 0 )
    {
      if ( cdoVerbose )
        cdoPrint("Since the project id is CMIP5 further attributes are tested. ");
      if ( ( errnum = check_attarray(kvl, reqAttCMIP5) ) != -1 )
        attErr(reqAttCMIP5, errnum);
#if ( CMOR_VERSION_MAJOR == 3 )
/**************/
/* Add additional attributes for CMIP5 */
/* allthough using CMOR 3 */
/**************/
      if ( ( errnum = check_attarray(kvl, reqAttCMIP5CMOR3) ) != -1 )
        attErr(reqAttCMIP5CMOR3, errnum);
#endif
    }
  else if (strcmp(project_id, "CORDEX") == 0 )
    {
      if ( cdoVerbose )
        cdoPrint("Since the project id is CORDEX further attributes are tested.");
      if ( ( errnum = check_attarray(kvl, reqAttCORDEX) ) != -1 )
        attErr(reqAttCORDEX, errnum);
    }
  else if (strcmp(project_id, "CMIP6") == 0 )
    {
      kv_insert_a_val(kvl, "_cmip6_option", (char *)"CMIP6", 1);
      if ( cdoVerbose )
        cdoPrint("Since the project_id is CMIP6 further attributes are tested.");
      if ( ( errnum = check_attarray(kvl, reqAttCMIP6) ) != -1 )
        attErr(reqAttCMIP6, errnum);
      int j = 0;
      while ( expdepAttCMIP6[j] != NULL )
        {
          keyValues_t *kv_reqatt = kvlist_search(kvl, expdepAttCMIP6[j]);
          if ( !kv_reqatt || strcmp(kv_reqatt->values[0], "notSet") == 0 )
            if ( cdoVerbose )
              cdoPrint("Depending on the experiment, attribute '%s' may be required. Either it is missing or notSet", expdepAttCMIP6[j]);
          else if ( cdoVerbose )
            cdoPrint("Attribute '%s' is '%s'. ", expdepAttCMIP6[j], kv_reqatt->values[0]);
          j++;
        }
    }
  return 1;
} 

static int check_mem(list_t *kvl, char *project_id)
{
  const char *ripcharCMIP5[] = {"realization", "initialization_method", "physics_version"};
  const char *ripcharCMIP6[] = {"realization_index","initialization_index","physics_index","forcing_index"};
  char ripchar[5];
  strcpy(ripchar, "ripf");

  char *kv_member = kv_get_a_val(kvl, "member", "");
  int ripcharlen = 3;
  if ( strcmp(project_id, "CMIP6") == 0 )
    {
      kv_member= kv_get_a_val(kvl, "variant_label", "");
      kv_insert_a_val(kvl, (char *)"member", kv_member, 1);
      ripcharlen = 4;
    }
  int memberlen = strlen(kv_member);
  char ripvaluechar[memberlen][CMOR_MAX_STRING];
  for ( int k = 0; k < memberlen; k++ )
    strcpy(ripvaluechar[k], kv_member);

  int firstNum, tonull, j = 0;
  for ( int i = 0; i < ripcharlen; i++ )
    {
      while ( kv_member[j] != ripchar[i] && j < memberlen )
        j++;
      if ( j == memberlen )
        {
          j = -1;
          if ( ripcharlen == 3)
            cdoPrint("Attribute 'member' has no RIP format! Default setting is used: \n          realization=-1 \n          initialization=-1 \n          physics=-1. \n          forcing=-1.");
          else
            cdoPrint("Attribute 'variant_label' has no RIPF format! Default setting is used: \n          realization=-1 \n          initialization=-1 \n          physics=-1. \n          forcing=-1.");
          break;
        }
      firstNum = j+1;
      j = 0;
      if ( i < (ripcharlen - 1) )
        while ( kv_member[j] != ripchar[i+1] && j < memberlen )
          j++;
      if ( j == memberlen )
        {
          j = -1;
          if ( ripcharlen == 3)
            cdoPrint("Attribute 'member' has no RIP format! Default setting is used: \n          realization=-1 \n          initialization=-1 \n          physics=-1. \n          forcing=-1.");
          else
            cdoPrint("Attribute 'variant_label' has no RIPF format! Default setting is used: \n          realization=-1 \n          initialization=-1 \n          physics=-1. \n          forcing=-1.");
          break;
        }
      tonull = j;
      j = 0;
      char *temp = ripvaluechar[i];
      temp += firstNum;
      strcpy(ripvaluechar[i], temp);

      if ( i < ripcharlen - 1 )
        ripvaluechar[i][tonull-firstNum] = '\0';
    }
  if ( j != -1 )
    {
      if ( ripcharlen == 3 )
        for ( int i = 0; i < ripcharlen; i++ )
          kv_insert_a_val(kvl, (char *)ripcharCMIP5[i], ripvaluechar[i], 1);
      else
        for ( int i = 0; i < ripcharlen; i++ )
          kv_insert_a_val(kvl, (char *)ripcharCMIP6[i], ripvaluechar[i], 1);
    }
  else
    {
      if ( ripcharlen == 3 )
        for ( int i = 0; i < ripcharlen; i++ )
          kv_insert_a_val(kvl, (char *)ripcharCMIP5[i], (char *)"-1", 1);
      else
        for ( int i = 0; i < ripcharlen; i++ )
          kv_insert_a_val(kvl, (char *)ripcharCMIP6[i], (char *)"-1", 1);
    }
  return 0;

/*
  char workchar[CMOR_MAX_STRING]; 
  int realization, initialization_method, physics_version, forcing;
  int ipos=0, ppos=0;

/* Test for the right member, else abort or warn */
/*
  
  if ( strlen(kv_member) >= 6 && kv_member[0] == 'r' )
    {
      strcpy(crealiz, &kv_member[1]);
      if ( strtok_r(crealiz, "i", &cinitial) )
        {
          strtok_r(cinitial, "p", &cphysics); 
          realization = strtol(crealiz, NULL, 10);
          initialization_method = strtol(cinitial, NULL, 10);
          physics_version = strtol(cphysics, NULL, 10);
        }
      else cphysics=NULL;
    }
  else {crealiz[0] = '\0'; cinitial[0] = '\0'; cphysics[0] = '\0';};
  if ( realization && initialization_method && physics_version)
    {
      char *ripvaluechar[] = {crealiz, cinitial, cphysics};
      for ( int i = 0; i < 3; i++ )
        kv_insert_a_val(kvl, ripchar[i], ripvaluechar[i], 1);
      return 1;
    }
  else if ( strcmp(kv_member, "notSet") == 0 )
    {
      cdoPrint("The member has no RIP format! Default setting is used: \n          realization=-1 \n           initialization_method=-1 \n           physics_version=-1. ");
      for ( int i = 0; i < 3; i++ )   
        kv_insert_a_val(kvl, ripchar[i], (char *)"-1", 1);
    }
/* Now abort or warn */
/* 
  if (strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CORDEX") == 0)
    cdoAbort("Attribute member has no RIP format (at least 6 characters and in RIP order)! Found for \n          member: %s. This is interpreted as \n           Realization: %s \n           Initialization: %s \n           Physics: %s \n             But three Integers are needed", kv_member, crealiz, cinitial, cphysics);
*/
} 


/*
static void dump_global_attributes(list_t *pml, int streamID)
{
  int natts;
  int vlistID = pstreamInqVlist(streamID);
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char *value = NULL;
      char buffer[8];
      int type, len;
      cdiInqAtt(vlistID, CDI_GLOBAL, i, name, &type, &len);
      switch ( type )
        {
        case CDI_DATATYPE_TXT:
          value = Malloc(len + 1);
          cdiInqAttTxt(vlistID, CDI_GLOBAL, name, len, value);
          value[len] = '\0';
          break;
        case CDI_DATATYPE_INT32:
          value = Malloc(CDI_MAX_NAME);
          cdiInqAttInt(vlistID, CDI_GLOBAL, name, len, (int *)buffer);
          snprintf(value, CDI_MAX_NAME, "%i", *(int *)buffer);
          break;
        case CDI_DATATYPE_FLT64:
          value = Malloc(CDI_MAX_NAME);
          cdiInqAttFlt(vlistID, CDI_GLOBAL, name, len, (double *)buffer);
          snprintf(value, CDI_MAX_NAME, "%e", *(double *)buffer);
          break;
        default:
          cdoWarning("Unsupported type %i name %s\n", type, name);
        }
      hinsert(ht, name, value);
      if ( value ) Free(value);
    }
}
*/

static void dump_special_attributes(list_t *kvl, int streamID)
{
  int vlistID = pstreamInqVlist(streamID);
  int fileID = pstreamFileID(streamID);
  size_t old_historysize;
  char *new_history = kv_get_a_val(kvl, "history", NULL);
  size_t historysize;
  int natts;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  if ( natts > 0 )
    old_historysize = (size_t) streamInqHistorySize(fileID);
  else
    old_historysize = 0;

  if ( old_historysize )
    {
      historysize = old_historysize;
      if ( new_history )
        historysize += strlen(new_history) + 1;
    }
  else if ( new_history )
    {
      historysize = strlen(new_history);
    }

  if ( historysize )
    {
      char *history = (char *)Malloc(historysize + 1);
      memset(history, 0, historysize + 1);
      if ( old_historysize )
        {
          streamInqHistoryString(fileID, history);
          if ( new_history )
            {
              strcat(history, " ");
              strcat(history, new_history);
            }
        }
      else
        {
          strcpy(history, new_history);
        }
      kv_insert_a_val(kvl, "history", history, 1);
      Free(history);
    }
}

static void read_config_files(list_t *kvl)
{
  if ( cdoVerbose )
    cdoPrint("1. Start to read configuration files.");
  /* Files from info key in command line. */
  keyValues_t *info = kvlist_search(kvl, "i");
  int i = 0;
  if ( info )
    while ( i < info->nvalues )
      {
        if ( cdoVerbose )
          cdoPrint("1.1. Try to parse file: '%s' configured with key 'info'.", info->values[i]);
        if ( parse_kv_file(kvl, info->values[i]) == 0 )
          cdoAbort("File '%s' does not exist.", info->values[i]);
        if ( cdoVerbose )
          cdoPrint("1.1. Successfully parsed file: '%s' configured with key 'info'.", info->values[i]);
        i++;
      }

  /* Config file in user's $cwd directory. */
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  const char *dotconfig = ".cdocmorinfo";
  char *workfile = (char *)Malloc((strlen(cwd) + strlen(dotconfig) + 2 ) * sizeof(char));
  sprintf(workfile, "%s/%s", cwd, dotconfig);
  if ( cdoVerbose )
    cdoPrint("1.2. Try to parse default file: '%s'.", workfile);
  if ( parse_kv_file(kvl, workfile) == 0 )
    cdoWarning("Default file for keyword 'info': '%s' does not exist.", workfile);
  else if ( cdoVerbose )
    cdoPrint("1.2. Successfully parsed default file: '%s'.", workfile);
  Free(workfile);
  
  if ( i == 0 )
    {
      keyValues_t *info2 = kvlist_search(kvl, "i");
      if ( info2 )
        while ( i < info2->nvalues )
          {
            if ( cdoVerbose )
              cdoPrint("1.3. Try to parse file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
            if ( parse_kv_file(kvl, info2->values[i]) == 0 )
              cdoAbort("File '%s' does not exist.", info2->values[i]);
            if ( cdoVerbose )
              cdoPrint("1.3. Successfully parsed file: '%s' configured with key 'info' in file '.cdocmorinfo'.", info2->values[i]);
            i++;
          }
    }
  if ( cdoVerbose )
    cdoPrint("1. Successfully read configuration files.");
}

static int in_list(char **list, const char *needle, int num)
{
  for ( int i = 0; i < num; i++ )
    if ( strcmp(list[i], needle) == 0 )
      return 1;
  return 0;
}

static int get_netcdf_file_action(list_t *kvl)
{
  char *proj = kv_get_a_val(kvl, "project_id", NULL);
  char *chunk = kv_get_a_val(kvl, "om", "r");
  if ( strcmp(proj, "CORDEX") == 0 )
    {
      if ( chunk[0] == 'r' )
        return CMOR_REPLACE_4;
      else if ( chunk[0] == 'a')
        return CMOR_APPEND_4;
      else if ( chunk[0] == 'p')
        return CMOR_PRESERVE_4;
      else
        {
          cdoWarning("You set output_mode = '%s', but valid are 'a' for append ,'r' for replace or 'p' for preserve.\n          CMOR output mode is set to: replace.", chunk);
          return CMOR_REPLACE;
        }
    }
  else
    {
      if ( chunk[0] == 'r' )
        return CMOR_REPLACE;
      else if ( chunk[0] == 'a')
        return CMOR_APPEND;
      else if ( chunk[0] == 'p')
        return CMOR_PRESERVE;
      else
        {
          cdoWarning("You set output_mode = '%s', but valid are 'a' for append ,'r' for replace or 'p' for preserve.\n          CMOR output mode is set to: replace.", chunk);
          return CMOR_REPLACE;
        }
    }
}

static int get_cmor_verbosity(list_t *kvl)
{
  char *verbos = kv_get_a_val(kvl, "set_verbosity", NULL);
  if ( !verbos )
    return CMOR_NORMAL;
  if ( strcmp(verbos, "CMOR_QUIET") == 0 )
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int get_cmor_exit_control(list_t *kvl)
{
  char *exit = kv_get_a_val(kvl, "exit_control", NULL);
  if ( !exit )
    return CMOR_NORMAL;
  if ( strcasecmp(exit,"CMOR_EXIT_ON_MAJOR") == 0 )
    return CMOR_EXIT_ON_MAJOR;
  else if ( strcasecmp(exit,"CMOR_EXIT_ON_WARNING")  == 0 )
    return CMOR_EXIT_ON_WARNING;
  else
    return CMOR_NORMAL;
}

static char *get_calendar_ptr(int calendar)
{
  char *calendar_ptr = (char *)Malloc(CMOR_MAX_STRING * sizeof(char));
  switch ( calendar )
    {
    case CALENDAR_STANDARD:
    case CALENDAR_GREGORIAN:
      strcpy(calendar_ptr, "gregorian"); break;
    case CALENDAR_PROLEPTIC:
      strcpy(calendar_ptr, "proleptic_gregorian"); break;
    case CALENDAR_360DAYS:
      strcpy(calendar_ptr, "360_day"); break;
    case CALENDAR_365DAYS:
      strcpy(calendar_ptr, "noleap"); break;
    case CALENDAR_366DAYS:
      strcpy(calendar_ptr, "all_leap"); break;
    default:
      Free(calendar_ptr); return NULL;
    }
  return calendar_ptr;
}

static int get_calendar_int(char *calendar)
{
  if ( !calendar )
    return 0;
  if ( strcmp(calendar, "gregorian") == 0 )
    return CALENDAR_STANDARD;
  else if ( strcmp(calendar, "proleptic_gregorian") == 0 )
    return CALENDAR_PROLEPTIC;
  else if ( strcmp(calendar, "360_day") == 0 )
    return CALENDAR_360DAYS;
  else if ( strcmp(calendar, "noleap") == 0 )
    return  CALENDAR_365DAYS;
  else if ( strcmp(calendar, "all_leap") == 0 )
    return  CALENDAR_366DAYS;
  else
    {
      cdoWarning("You set calendar type = '%s' which is not supported by CMOR.", calendar);
      return 0;
    }
}

static char *get_txtatt(int vlistID, int varID, const char *key)
{
  int natts;
  cdiInqNatts(vlistID, varID, &natts);
  char *txtatt = NULL;
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char buffer[8];
      int type, len;
      cdiInqAtt(vlistID, varID, i, name, &type, &len);
      if ( strcmp(name, key) == 0 )
        {
          txtatt = (char *)Malloc(CMOR_MAX_STRING * sizeof(char));
          cdiInqAttTxt(vlistID, varID, name, len, txtatt);
          txtatt[len] = '\0';
          return txtatt;
        }
    }
  return txtatt;
}

/***********************************************/
/*Time related functions************************/
/***********************************************/

static char *get_time_units(int taxisID)
{
  char *units = (char *)Malloc ( CMOR_MAX_STRING * sizeof(char) );
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second;
  cdiDecodeDate(taxisInqRdate(taxisID), &year, &month, &day);
  cdiDecodeTime(taxisInqRtime(taxisID), &hour, &minute, &second);
  if ( timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES )
    timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_3HOURS || timeunit == TUNIT_6HOURS ||
       timeunit == TUNIT_12HOURS )
    timeunit = TUNIT_HOUR;

  sprintf(units, "%s since %d-%d-%d %02d:%02d:%02d\0", tunitNamePtr(timeunit),
          year, month, day, hour, minute, second);
  return units;
}

static int get_time_step_int(char *time_step)
{
  if ( strcmp(time_step, "hours") == 0  )
    return TUNIT_HOUR;
  else if ( strcmp(time_step, "days") == 0 ) 
    return TUNIT_DAY;
  else if ( strcmp(time_step, "months") ==  0 )
    return TUNIT_MONTH;
  else if ( strcmp(time_step, "years") == 0  )
    return TUNIT_YEAR;
  else
    {
      cdoWarning("You set required_time_units = '%s since...'.\n          This time step is not yet implemented in cmor.", time_step);
      return 0;
    }
}

static int check_time_units(char *time_units)
{
/* Required attribute in check_att */
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char time_step[CMOR_MAX_STRING];
  if ( sscanf(time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s",
                  time_step, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond) != 7)
    {
      cdoWarning("You set required_time_units = '%s'\n          but it requires the form 'timestep since year-month-day hour:minute:second.\n          Could not read all 7 required time unit values.", time_units);
      return 0;
    }
  if ( !get_time_step_int(time_step) )
    return 0;
  return 1;
}

static void get_time_method(list_t *kvl, int vlistID, int varID, char *cmor_time_name, char *project_id, int miptab_freq, int *time_axis)
{
  if ( ( strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CMIP6") == 0 ) && miptab_freq )
    switch ( miptab_freq )
      {
      case 1: strcpy(cmor_time_name, "time2"); *time_axis=2; break;
      case 2: strcpy(cmor_time_name, "time"); *time_axis=0; break;
      case 4: strcpy(cmor_time_name, "time"); *time_axis=0; break;
      case 5: strcpy(cmor_time_name, "time1"); *time_axis=1; break;
      case 6: strcpy(cmor_time_name, "time1"); *time_axis=1; break;
      case 7: strcpy(cmor_time_name, "time3"); *time_axis=3; break;
      }
  if ( cmor_time_name[0] != 't' )
    {
      char *time_method = get_txtatt(vlistID, varID, "cell_methods");
      char *att_time_method = kv_get_a_val(kvl, "cm", NULL);
      check_compare_set(&time_method, att_time_method, "cell_methods", " ");
      if ( time_method[0] == 'm' )      { strcpy(cmor_time_name, "time \0"); *time_axis=0; }
      else if ( time_method[0] == 'p' ) { strcpy(cmor_time_name, "time1\0"); *time_axis=1; }
      else if ( time_method[0] == 'c' ) { strcpy(cmor_time_name, "time2\0"); *time_axis=2; }
      else if ( time_method[0] == 'd' ) { strcpy(cmor_time_name, "time3\0"); *time_axis=3; }
      else if ( time_method[0] == 'n' ) { strcpy(cmor_time_name, "none\0"); *time_axis=4; }
      else
        {
          cdoWarning("You set cell method = '%s' which is not valid. Check CF-conventions for allowed time cell methods.\n          Time cell method is set to 'mean'. ", time_method);
          strcpy(cmor_time_name, "time \0");
        }
      Free(time_method);
    }
}

static void get_taxis(char *required_time_units, int *sdate, int *stime, int *timeunit)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char atttimeunit[CMOR_MAX_STRING];

  sscanf(required_time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s",
                  atttimeunit, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond);
  *sdate = cdiEncodeDate(attyear, attmonth, attday);
  *stime = cdiEncodeTime(atthour, attminute, attsecond);
  *timeunit = get_time_step_int(atttimeunit);
}

static double *get_branch_times(list_t *kvl, int calendar, char *time_units)
{
  if ( cdoVerbose )
    cdoPrint("6.1.3. Start to compute attribute 'branch_time'.");
  double *branch_time = (double *)Malloc(2 * sizeof(double));
  branch_time[0] = 0.0;
  branch_time[1] = 0.0;

  int numdates;
  char **branch_times_p = kv_get_vals(kvl, "branch_times", &numdates);

  if ( numdates == 2 && kv_get_a_val(kvl, "parent_experiment_id", NULL) )
    {
      int parentdates[2], parentyears[2], parentmonths[2], parentdays[2], parentvdates[2], parentvtimes[2];
      juldate_t parentstartdate, parentbranchdate, childstartdate ;
      for ( int i = 0; i < 2; i++ )
        {
          parentdates[i]  = atol(branch_times_p[i]);
          parentyears[i]  =  parentdates[i]/100/100;
          parentmonths[i] = (parentdates[i] - parentyears[i]*100*100)/100;
          parentdays[i]   =  parentdates[i] - parentyears[i]*100*100 - parentmonths[i]*100;
          parentvdates[i] = cdiEncodeDate(parentyears[i], parentmonths[i], parentdays[i]);
          parentvtimes[i] = 0;
        }
      parentstartdate   = juldate_encode(calendar, parentvdates[0], parentvtimes[0]);
      parentbranchdate  = juldate_encode(calendar, parentvdates[1], parentvtimes[1]);

      int childsdate, childstime, childtimeunit;
      get_taxis(time_units, &childsdate, &childstime, &childtimeunit);
      childstartdate    = juldate_encode(calendar, childsdate, childstime);

/* If time unit is always "days since.." */
      branch_time[0] = juldate_to_seconds(juldate_sub(parentbranchdate, parentstartdate)) / 86400;
      branch_time[1] = juldate_to_seconds(juldate_sub(parentbranchdate, childstartdate)) / 86400;
    }
  if ( cdoVerbose )
    cdoPrint("6.1.3. Successfully computed 'branch_time': '%f', '%f'.", branch_time[0], branch_time[1]);
  return branch_time;
}

static char *check_required_time_units(list_t *kvl, int taxisID)
{
  if ( cdoVerbose )
    cdoPrint("6.1.2. Start to check attribute 'required_time_units'.");
  char *time_units = get_time_units(taxisID);
  char *required_time_units = kv_get_a_val(kvl, "required_time_units", NULL);
  if ( check_time_units(required_time_units) )
    check_compare_set(&time_units, required_time_units, "time_units", NULL);
  else 
    cdoAbort("Required Attribute 'required_time_units' from configuration is invalid!");
  kv_insert_a_val(kvl, "required_time_units", time_units, 1);
  if ( cdoVerbose )
    cdoPrint("6.1.2. Successfully checked attribute 'required_time_units'.");
  return time_units;
}

static char *check_calendar(list_t *kvl, int taxisID, int *calendar)
{ 
  if ( cdoVerbose )
    cdoPrint("6.1.1. Start to check attribute 'calendar'.");
  char *attcalendar = kv_get_a_val(kvl, "calendar", NULL);
  char *calendarptr = get_calendar_ptr(taxisInqCalendar(taxisID));
  if ( *calendar = get_calendar_int(attcalendar) )
    check_compare_set(&calendarptr, attcalendar, "calendar", NULL);
  else 
    {
      if ( cdoVerbose )
        cdoPrint("Try to use infile calendar.");
      if ( !get_calendar_int(calendarptr) )
        cdoAbort("In validating calendar:\n          No valid configuration calendar and no valid infile calendar found.");
      else
        *calendar = get_calendar_int(calendarptr);
    }
  if ( cdoVerbose )
    cdoPrint("6.1.1. Successfully checked attribute 'calendar'.");
  return calendarptr;
}

/*********/
/* main: */
/*********/

static void setup_dataset(list_t *kvl, int streamID, int *calendar)
{
  if ( cdoVerbose )
    cdoPrint("6. Start to process cmor_setup and cmor_dataset.");
  int netcdf_file_action = get_netcdf_file_action(kvl);
  int set_verbosity = get_cmor_verbosity(kvl);
  int exit_control = get_cmor_exit_control(kvl);
  int creat_subs = 1;
  char *drs = kv_get_a_val(kvl, "d", "y");
  if ( drs[0] == 'n' )
    creat_subs = 0;
  else if ( drs[0] != 'y' )
    {
      cdoWarning("In preparing cmor_setup:\n          You set 'd' = '%s' which is not valid.\n          Allowed are: 'n' or 'y'. d is set to 'y'.", drs);
      kv_insert_a_val(kvl, "d", (char *)"y", 1);
    }

  int vlistID = pstreamInqVlist(streamID);

  int cmf = cmor_setup(kv_get_a_val(kvl, "inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             kv_get_a_val(kvl, "logfile", NULL),
             &creat_subs);
  if ( cmf != 0 )
    cdoAbort("Function cmor_setup failed!");

  int taxisID = vlistInqTaxis(vlistID);

/*
  char *attcomment = kv_get_a_val(kvl, "comment", NULL);
  char *comment = get_txtatt(vlistID, CDI_GLOBAL, "comment");
*/

/* First compare file calendar and config calendar and retrieve pointer and integer
   Then check the required time units from config and retrieve
   Then compute branch_time_in_parent and branch_time_in_child */

  if ( cdoVerbose )
    cdoPrint("6.1. Start to check model calendar as well as 'required_time_units' and 'branch_times' attributes.");
  char *calendarptr = check_calendar(kvl, taxisID, calendar);  
  char *time_units = check_required_time_units(kvl, taxisID);
  double *branch_times = get_branch_times(kvl, *calendar, time_units); 
  if ( cdoVerbose )
    cdoPrint("6.1. Successfully found valid calendar, 'required_time_units' and 'branch_times'.");
#if defined(CMOR_VERSION_MAJOR)
#if ( CMOR_VERSION_MAJOR == 2 )
    {
      cmf = cmor_dataset(kv_get_a_val(kvl, "dr", "./"),
               kv_get_a_val(kvl, "experiment_id", ""),
               kv_get_a_val(kvl, "institution", ""),
               kv_get_a_val(kvl, "source", ""),
               calendarptr,
               atoi(kv_get_a_val(kvl, "realization", "")),
               kv_get_a_val(kvl, "contact", ""),
               kv_get_a_val(kvl, "history", ""),
               kv_get_a_val(kvl, "comment", ""),
               kv_get_a_val(kvl, "references", ""),
               atoi(kv_get_a_val(kvl, "leap_year", "")),
               atoi(kv_get_a_val(kvl, "leap_month", "")),
               NULL,
               kv_get_a_val(kvl, "model_id", ""),
               kv_get_a_val(kvl, "forcing", ""),
               atoi(kv_get_a_val(kvl, "initialization_method", "")),
               atoi(kv_get_a_val(kvl, "physics_version", "")),
               kv_get_a_val(kvl, "institute_id", ""),
               kv_get_a_val(kvl, "parent_experiment_id", ""),
               &(branch_times[0]),
               kv_get_a_val(kvl, "parent_experiment_rip", ""));
      if ( cmf != 0 )
        cdoAbort("Function cmor_dataset failed!");
    }
  const char *allneeded2[] = {"cordex_domain",  "driving_experiment", "driving_model_id", "driving_model_ensemble_member", "driving_experiment_name", "rcm_version_id", NULL};
  int ind = 0;
  if ( strcmp(kv_get_a_val(kvl, "project_id", NULL),"CORDEX") == 0 )
    while ( allneeded2[ind] )
      {
        char *tmp = kv_get_a_val(kvl, allneeded2[ind], NULL );
        if ( tmp )
          cmf = cmor_set_cur_dataset_attribute((char *)allneeded2[ind], tmp, 1);
        if ( cmf != 0 )
          cdoAbort("Function cmor_set_cur_dataset_attribute failed!");
        ind++;
      }
#elif ( CMOR_VERSION_MAJOR == 3 )
    {
/***/
/* Could not give CMOR all attributes separately because some are required to be in a json file (outpath,...). /
/* Better collect them in this file. */
/* todo this **/
/* If a Json file is denoted, read this file and check attributes */
/***/

/*
      char *filename = kv_get_a_val(kvl, "dj", NULL); */

      char cwd[1024];
      getcwd(cwd, sizeof(cwd));
      char *dataset_path = (char *) Malloc( (strlen(cwd) + 1 + strlen("dataset.json") + 3) * sizeof(char));;
      FILE *dataset_json;
      int procID = getpid();

      sprintf(dataset_path, "%s/dataset%d.json\0", cwd,procID);

      dataset_json = fopen(dataset_path, "w+");
      if ( !dataset_json )
        cdoAbort("In preparing cmor_dataset:\n          Could not open a dataset file '%s' for cmor_dataset.", dataset_path);
      fputs("{\n", dataset_json);

      const char *allneeded[] = /*CMIP5*/{"project_id", "experiment_id", "institution", "source", "realization", "contact", "history", "comment", "references", "leap_year", "leap_month", "source_id", "model_id", "forcing", "initialization_method", "modeling_realm", "physics_version", "institute_id", "parent_experiment_rip", 
/*CORDEX */
  "CORDEX_domain",  "driving_experiment", "driving_model_id", "driving_model_ensemble_member", "driving_experiment_name", "rcm_version_id",
/* CMIP6: */
  /* Glob Atts */
"_cmip6_option", "Conventions", "activity_id", "branch_method", "experiment", "experiment_id", "forcing_index", "further_info_url", "grid", "grid_label", "initialization_index", "institution", "institution_id", "license", "mip_era", "nominal_resolution", "physics_index", "product", "realization_index", "source", "source_id", "source_type", "sub_experiment", "sub_experiment_id", "table_id", "variant_label", "parent_experiment_id", "parent_activity_id", "parent_mip_era", "parent_source_id", "parent_variant_label",
"parent_time_units", NULL};
      int i = 0;

      while ( allneeded[i] )
        {
          char *tmp = kv_get_a_val(kvl, allneeded[i], "notSet");
          if ( strncmp(tmp, "notSet", 6) != 0 )
            {
                  int linelen = strlen(allneeded[i]) + strlen(tmp) + 10;
                  char line[linelen];
                  sprintf(line, "\"%s\" : \"%s\",\n", allneeded[i], tmp);
                  fputs((const char *) line, dataset_json);
            }
          i++;
        }

/*      char *miptabfreq = kv_get_a_val(kvl, "miptab_freq", NULL);
      char *freq = (char *) Malloc(11 * sizeof(char));
      if ( strstr(miptabfreq, "yr") || strstr(miptabfreq, "Yr") )
        strcpy(freq, "yr");
      else if ( strstr(miptabfreq, "mon") || strstr(miptabfreq, "Mon") )
        {
          strcpy(freq, "mon");
          if ( strstr(miptabfreq, "ClimMon") )
            strcpy(freq, "ClimMon");
        }
      else if ( strstr(miptabfreq, "day") || strstr(miptabfreq, "Day") )
        strcpy(freq, "day");
      else if ( strstr(miptabfreq, "6hr") )
        strcpy(freq, "6hr");
      else if ( strstr(miptabfreq, "3hr") )
        strcpy(freq, "3hr");
      else if ( strstr(miptabfreq, "1hr") )
        {
          strcpy(freq, "1hr");
          if ( strstr(miptabfreq, "1hrClimMon") )
            strcpy(freq, "1hrClimMon");
        }
      else if ( strstr(miptabfreq, "fx") )
        strcpy(freq, "fx"); */

      char *branch_time_in_parent = (char *) Malloc(2*sizeof(double));
      char *branch_time_in_child = (char *) Malloc(2*sizeof(double));
      snprintf(branch_time_in_parent, sizeof(double), "%.12f", branch_times[0]);
      snprintf(branch_time_in_child, sizeof(double), "%.12f", branch_times[1]);

  /* CMOR internal */ 
          fputs("\"outpath\" : \"", dataset_json);
          fputs(kv_get_a_val(kvl, "dr", "./"), dataset_json);
          fputs("\",\n", dataset_json); 
          fputs("\"output_path_template\" : \"", dataset_json);
          fputs(kv_get_a_val(kvl, "output_path_template", "<mip_era><activity_id><institution_id><source_id><experiment_id><variant_label><table><variable_id><grid_label><version>"), dataset_json);
          fputs("\",\n", dataset_json); 
          fputs("\"output_file_template\" : \"", dataset_json);
          fputs(kv_get_a_val(kvl, "output_file_template", "<variable_id><table><source_id><experiment_id><variant_label><grid_label>"), dataset_json);
          fputs("\",\n", dataset_json); 
          fputs("\"tracking_prefix\" : \"", dataset_json);
          fputs(kv_get_a_val(kvl, "tracking_prefix", "hdl:21.14100"), dataset_json);
          fputs("\",\n", dataset_json); 
/*          fputs("\"frequency\" : \"", dataset_json);
          fputs(freq, dataset_json);
          fputs("\",\n", dataset_json);  */

/* cdo cmor preprocessed: */
          fputs("\"calendar\" : \"", dataset_json);
          fputs(calendarptr, dataset_json);
          fputs("\",\n", dataset_json); 
          fputs("\"branch_time_in_parent\" : \"", dataset_json);
          fputs(branch_time_in_parent, dataset_json);
          fputs("\",\n", dataset_json);
          fputs("\"branch_time_in_child\" : \"", dataset_json);
          fputs(branch_time_in_child, dataset_json);
          fputs("\",\n", dataset_json); 

          fputs("}\n", dataset_json);
          fclose(dataset_json);
          cmf = cmor_dataset_json(dataset_path);
          if ( cmf != 0 )
            cdoAbort("Function cmor_dataset_json failed!");
          Free(dataset_path);


      Free(branch_time_in_parent);
      Free(branch_time_in_child);
/*      Free(freq); */
    }
#else
    cdoAbort("Cmor version %d not yet enabled!", (int) CMOR_VERSION_MAJOR);
#endif
#else
    cdoAbort("It is not clear which CMOR version is installed since\n          Makros CMOR_VERSION_MAJOR and CMOR_VERSION_MINOR are not available.");
#endif
  Free(calendarptr);
  Free(branch_times);
  Free(time_units);
  if ( cdoVerbose )
    cdoPrint("6. Successfully finished cmor_setup and cmor_dataset.");
}


static void gen_bounds(int n, double *vals, double *bounds)
{
  for ( int i = 0; i < n-1; ++i )
    {
      bounds[2*i+1]   = 0.5*(vals[i] + vals[i+1]);
      bounds[2*(i+1)] = 0.5*(vals[i] + vals[i+1]);
    }

  bounds[0]     = 2*vals[0] - bounds[1];
  bounds[2*n-1] = 2*vals[n-1] - bounds[2*(n-1)];
}

static void get_zcell_bounds(int zaxisID, double *zcell_bounds, double *levels, int zsize)
{
  double *lbounds;
  lbounds = (double *)Malloc(zsize * sizeof(double));
  zaxisInqLbounds(zaxisID, lbounds);
  double *ubounds;
  ubounds = (double *)Malloc(zsize * sizeof(double));
  zaxisInqUbounds(zaxisID, ubounds);
  if ( !lbounds || !ubounds || pow((ubounds[1] - ubounds[0]),2) < 0.001 || pow((lbounds[1] - lbounds[0]), 2) < 0.001 )
    gen_bounds(zsize, levels, zcell_bounds);
  else
    {
      if ( lbounds )
        zcell_bounds[0] = lbounds[0];
      else
        zcell_bounds[0] = 0; 
      for ( int i = 0; i < zsize-1; ++i )
        {
          zcell_bounds[2*i+1]   = ubounds[i];
          zcell_bounds[2*(i+1)] = lbounds[i+1];
        }
      if ( ubounds )
        zcell_bounds[2*zsize-1] = ubounds[zsize-1];
      else
        zcell_bounds[2*zsize-1] = levels[zsize-1] + ( levels[zsize-1] - zcell_bounds[2*zsize-2] );
    }
  Free(lbounds);
  Free(ubounds);
}

static void get_zhybrid(int zaxisID, double *p0, double *alev_val, double *alev_bnds, double *b_val, double *b_bnds, double *ap_val, double *ap_bnds)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  double *vct = (double *)Malloc(vctsize * sizeof(double) );
  zaxisInqVct(zaxisID, vct);
  for ( int i = 0; i<(zsize+1); i++)
    {
      ap_bnds[i] = vct[i];
      b_bnds[i] = vct[zsize+1+i];
    }
  for ( int i = 0; i<zsize; i++)
    {
      ap_val[i] = (ap_bnds[i]+ ap_bnds[i+1]) / 2.0;
      b_val[i] = (b_bnds[i]+ b_bnds[i+1]) / 2.0;
      alev_val[i] = ap_val[i]/p0[0] + b_val[i];
      alev_bnds[i] = ap_bnds[i]/p0[0] + b_bnds[i];
    }
  alev_bnds[zsize] = ap_bnds[zsize]/p0[0] + b_bnds[zsize];
  Free(vct);
}

static int get_strmaxlen(char **array, int len)
{
  int result = 0, i;
  for (i = 0; i < len; i++)
    if ( result < strlen(array[i]) )
      result = strlen(array[i]);
  return result;     
}

static void register_z_axis(list_t *kvl, int vlistID, int varID, int zaxisID, char *varname, int *axis_ids, int *zfactor_id, char *project_id, int miptab_freq)
{

  *zfactor_id = 0;
  int cmf = 0;
  int zsize = zaxisInqSize(zaxisID);
  double *levels;

  char *chardimatt = kv_get_a_val(kvl, "ca", NULL);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");
  if ( strcmp(chardim, "vegtype") == 0 || strcmp(chardim, "oline") == 0  )
    {
      if ( zsize )
        cdoPrint("In Z-axis registration:\n          You configured a character coordinate '%s' but a zaxis is found with '%d' numerical values. The zaxis attributes are ignored and the '%d' levels are interpreted as the character coordinates in the order they are given for '%s'.", chardim, zsize, zsize, varname);
      int numchar = 0;
      char *charvalstring = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(charvalstring, "char_axis_%s", chardim);
      char **charvals = kv_get_vals(kvl, charvalstring, &numchar);
      Free(charvalstring);
      if ( charvals )
        {
          int maxlen = get_strmaxlen(charvals, numchar);
          void *charcmor = (void *) Malloc ( (numchar * maxlen + 1) * sizeof(char));
          sprintf((char *)charcmor, "%s", charvals[0]);
          char blanks[maxlen];
          for ( int i = 0; i < maxlen; i++)
            blanks[i] = ' ';
          sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[0]), blanks);         
          for ( int i = 1; i < numchar; i++ )
            {
              sprintf((char *)charcmor, "%s%s", (char *)charcmor, charvals[i]);
              sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[i]), blanks);         
            }
          if ( numchar == zsize )
            cmf = cmor_axis(new_axis_id(axis_ids), chardim, (char *)"", numchar, (void *)charcmor, 'c',  NULL, maxlen, NULL); 
          else
            cdoAbort("In registration of a character coordinate as a substitution for a vertical coordinate:\n          The number of registered character coordinates '%d' differ from the number of axis levels '%d'.", numchar, zsize);
          if ( cmf != 0 )
            cdoAbort("Function cmor_axis failed!");
          Free(charcmor);
        }
      else
        cdoAbort("In registration of a character coordinate as a substitution for a vertical coordinate:\n          You configured a character coordinate '%s' but no values are found! Configure values via attribute 'char_dim_vals'.", chardim);
    }
  else
  {
  if ( zsize > 1)
    {
      levels = (double *)Malloc(zsize * sizeof(double));
      zaxisInqLevels(zaxisID, levels);
      double *zcell_bounds;
      zcell_bounds = (double *)Malloc( 2*zsize * sizeof(double) );
      get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
      if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
        {
          char *zaxis = get_txtatt(vlistID, varID, "z_axis");
          char *attzaxis = kv_get_a_val(kvl, "za", NULL);
          check_compare_set(&zaxis, attzaxis, "z_axis", "notSet");

          if ( strcmp(zaxis, "notSet") != 0 )
            cmf = cmor_axis(new_axis_id(axis_ids),
                        zaxis,
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
          else if ( strcmp(project_id, "CMIP5") != 0 && strcmp(project_id, "CMIP6") != 0 )
            cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *)"plevs",
                        (char *)"Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
          else
            {
              if ( strcmp(project_id, "CMIP6") == 0 )
                {
                  switch ( miptab_freq )
                    {
                    case 12: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plev19",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    case 4: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plev8",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    default: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plevs",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    }
                }
              else
                {
                  switch ( miptab_freq )
                    {
                    case 3: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plev7",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    case 4: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plev8",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    case 5: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plev3",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    default: cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "plevs",
                        (char *) "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                    }         
                }       
            }
          if ( zaxis ) Free(zaxis);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
        {
          double *alev_val = (double *) Malloc(zsize * sizeof(double));
          double *alev_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
          double *ap_val = (double *) Malloc(zsize * sizeof(double));
          double *ap_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
          double *b_val = (double *) Malloc(zsize * sizeof(double));
          double *b_bnds = (double *) Malloc((zsize + 1) * sizeof(double));
          double *p0 = (double *) Malloc(sizeof(double));
          p0[0] = 101325.0;

          char *mtproof = kv_get_a_val(kvl, "mtproof", NULL);
          if ( mtproof )
            {
              if ( cdoVerbose )
                cdoPrint("Mapping table: '%s' is applied for ps.", mtproof);
              kv_insert_a_val(kvl, "workfile4err", mtproof, 1);
              list_t *pml = cdo_parse_cmor_file(mtproof, kvl);
              const char *tempo[] = {"ps"};
              maptab_via_cn(pml, (char **)tempo, vlistID, vlistNvars(vlistID), 1, kv_get_a_val(kvl, "miptab_freq", NULL),  FILETYPE_NC ); 
              list_destroy(pml);
            }
          else
            {
              cdoPrint("Ps needs to be one infile variable name.");
            }
          int psID = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
          if ( psID == CDI_UNDEFID )
            cdoAbort("In registration of a vertical axis:\n          Could not find a surface pressure variable in infile. Cannot register a hybrid zaxis without surface pressure.");

          get_zhybrid(zaxisID, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
/*cmor_zfactor (int *zfactor_id,int zaxis_id, char *zfactor_name, char *units, int ndims, int axis_ids[], char type, void *zfactor_values, void *zfactor_bounds)*/
          cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "alternate_hybrid_sigma",
                        (char *) "",
                        zsize,
                        (void *)alev_val,
                        'd', alev_bnds,  1, NULL);
          int lev_id = axis_ids[count_axis_ids(axis_ids)-1];
          int lev_id_array[2];
          lev_id_array[0] = lev_id;
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *)"p0", (char *)"Pa", 0, 0, 'd', (void *)p0, NULL);
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *)"b",  (char *)"", 1, &lev_id_array[0], 'd', (void *)b_val, (void *)b_bnds);
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *)"ap", (char *)"Pa", 1, &lev_id_array[0], 'd', (void *)ap_val, (void *)ap_bnds);
          cmf = cmor_zfactor(zfactor_id, lev_id, (char *)"ps", (char *)"Pa", count_axis_ids(axis_ids)-1, axis_ids, 'd', NULL, NULL);  
          Free(alev_val);  
          Free(alev_bnds);  
          Free(ap_val);  
          Free(ap_bnds);  
          Free(b_val);  
          Free(b_bnds);  
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA )
        {
          zcell_bounds[0] = (double) 0;
          cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "depth_coord",
                        (char *) "m",
                        zsize,
                        (void *)levels,
                        'd', zcell_bounds,  2, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_LAND )
        {
          zcell_bounds[0] = (double) 0;
          cmf = cmor_axis(new_axis_id(axis_ids),
                        (char *) "sdepth",
                        (char *) "cm",
                        zsize,
                        (void *)levels,
                        'd', zcell_bounds, 2, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_GENERIC || zaxisInqType(zaxisID) == ZAXIS_HEIGHT)
        {
          char *zaxisname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
          zaxisInqName(zaxisID, zaxisname);
          if ( strcmp(zaxisname, "rho") == 0 )
            {
              char *zaxisunits = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
              zaxisInqUnits(zaxisID, zaxisunits);
              if ( strcmp(zaxisunits, "kg m-3") != 0 )
                {
                  cdoAbort("For zaxis with name 'rho' the units must be kg m-3 but are: '%s'", zaxisunits);
                }
              else
                {
                  levels = (double *) Malloc(zsize * sizeof(double));
                  zaxisInqLevels(zaxisID, levels);
                  double *zcell_bounds;
                  zcell_bounds = (double *) Malloc( 2*zsize * sizeof(double) );
                  get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
                  cmf = cmor_axis(new_axis_id(axis_ids),
                      (char *) "rho",
                      (char *) "kg m-3",
                      zsize,
                      (void *) levels,
                      'd', zcell_bounds, 2, NULL);
                }
              Free(zaxisunits);
            }
          else
            cdoAbort("In registration of a vertical axis:\n          Z-axis type %d with name '%s' not yet enabled.", zaxisInqType(zaxisID), zaxisname);
          Free(zaxisname);
        }
      else
        cdoAbort("In registration of a vertical axis:\n          Invalid Z-axis type %d . ", zaxisInqType(zaxisID));
      Free(zcell_bounds);
      Free(levels);
    }
  char *szc_name = kv_get_a_val(kvl, "szc", NULL);
  if ( zsize == 1 &&  szc_name )
    {
      char *szc_value = NULL;
      strtok_r(szc_name, "_", &szc_value);
      if ( !szc_value || !szc_value[0] )
        cdoAbort("Could not find an underscore '_' in szc value '%s' to seperate axis name from axis value", szc_name);
      levels = (double *)Malloc(sizeof(double));
      levels[0] = (double) atof(szc_value);
      if ( cdoVerbose )
        cdoPrint("Attribute szc is found.\n          Scalar z coordinate name is: '%s'\n          Scalar z coordinate value is: '%f'\n          ", szc_name, levels[0]);
      cmf = cmor_axis(new_axis_id(axis_ids),
                      szc_name,
                      (char *) "m",
                      zsize,
                      (void *) levels,
                      'd', NULL, 0, NULL);
      Free(levels);
    }
  }
  Free(chardim);
  if ( cmf != 0 )
    cdoAbort("Function cmor_axis failed!");
}

/*
static void register_character_dimension(int *axis_ids, char *filename)
{
  printf("The grid type is generic and a dimension 'basin' is found.\nTherefore, it is tried to read the character dimension.\n");
  int nc_file_id, nfiledims, nvars, ngatts, unlimdimid;
  nc_type xtypep;
  int varndims, varnattsp;
  int *vardimids;

  char *varname = Malloc(36 * sizeof(char));
  char *dimname = Malloc(36 * sizeof(char));

  size_t dimlength, dimstrlength;

  nc_open(filename, NC_NOWRITE, &nc_file_id);
  nc_inq(nc_file_id, &nfiledims, &nvars, &ngatts, &unlimdimid);
  vardimids = Malloc(nfiledims * sizeof(int));
  void *final_chardim;
  for ( int i = 0; i < nvars; i++ )
    {
      nc_inq_var(nc_file_id, i, varname, &xtypep, &varndims, vardimids, &varnattsp);
      if ( strcmp(varname, "region") == 0 )
        {
          nc_inq_dim(nc_file_id, vardimids[1], dimname, &dimstrlength);
          nc_inq_dim(nc_file_id, vardimids[0], dimname, &dimlength);

          final_chardim = (void *)Malloc(dimstrlength * dimlength *sizeof(char));
          nc_get_var(nc_file_id, i, final_chardim);
        }
    }
  nc_close(nc_file_id);
  cmor_axis(new_axis_id(axis_ids), dimname, "", dimlength, final_chardim, 'c',  NULL, dimstrlength, NULL); 
  Free(varname);
  Free(dimname);
  Free(vardimids);
}
*/

static void change_grid(char *grid_file, int gridID, int vlistID)
{
  if ( cdoVerbose )
    cdoPrint("You configured a grid_info file: '%s'. It is tested for a valid use as substitution.\n");
  argument_t *fileargument = file_argument_new(grid_file);
  int streamID2 = pstreamOpenRead(fileargument); 
  int vlistID2 = pstreamInqVlist(streamID2);
  int gridID2 = vlistInqVarGrid(vlistID2, 0); 

  if ( !gridID2 )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because of internal problems.", grid_file);

  int a,b;
  a = gridInqSize(gridID);
  b = gridInqSize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because total size of $IFILE: '%d' is not identical to total size of ginfo file: '%d'.", grid_file, a, b);

  a = gridInqYsize(gridID);
  b = gridInqYsize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because ysize of $IFILE: '%d' is not identical to ysize of ginfo file: '%d'.", grid_file, a, b);

  a = gridInqXsize(gridID);
  b = gridInqXsize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n          because xsize of $IFILE: '%d' is not identical to xsize of ginfo file: '%d'.", grid_file, a, b);

  vlistChangeGrid(vlistID, gridID, gridID2);
  cdoPrint("Successfully substituted grid.");

  pstreamClose(streamID2);
}

static void move_lons(double *xcoord_vals, double *xcell_bounds, int xsize, int xboundsize, int xnbounds)
{  
  int testbool = 0;
  for ( int i = 0; i < xsize; i++)
    if ( xcoord_vals[i] < 0.0 )
      {
        testbool = 1;
        break;
      }
  if ( testbool > 0 )
    for ( int i = 0; i < xsize; i++ )
      if ( xcoord_vals[i] < 0 )
        xcoord_vals[i] += 360.0;
  if ( xnbounds > 1 && testbool > 0 )
    for ( int j = 0; j < xboundsize; j++ )
      if ( xcell_bounds[j] < 0 )
        xcell_bounds[j] += 360.0;
}

static void inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals, double *xcell_bounds, double *ycell_bounds)
{
  gridInqYvals(gridID, ycoord_vals);
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
}

static void get_cmor_table(list_t *kvl, char *project_id)
{
  int gridtable_id;
  int cmf = 0;
  char gridtable[CMOR_MAX_STRING];
  char *mip_table_dir = kv_get_a_val(kvl, "mip_table_dir", NULL);
  if ( mip_table_dir && project_id )
    {
#if ( CMOR_VERSION_MAJOR == 2 )
      sprintf(gridtable, "%s/%s_grids\0", mip_table_dir, project_id);
#elif ( CMOR_VERSION_MAJOR == 3 )
      sprintf(gridtable, "%s/%s_grids.json\0", mip_table_dir, project_id);
#endif
      if ( file_exist(gridtable, 0, "Cmor-grid_table") )  
        {
          cmf = cmor_load_table(gridtable, &gridtable_id);
          cmf = cmor_set_table(gridtable_id);
        }
      else
        cdoAbort("In grid registration:\n          A project grid table is required for this type of grid but not found in the mip table directory '%s'.", mip_table_dir);
    }
  else
    {
      cdoAbort("In grid registration:\n          A project grid table is required for this type of grid but not found in the mip table directory. Check attributes 'mip_table_dir' and 'project_id' !");
    } 
  if ( cmf != 0 )
    cdoAbort("Function cmor_load_table or cmor_set_table failed!"); 
}

static void check_and_gen_bounds(int gridID, int nbounds, int length, double *coord_vals, double *cell_bounds, int x)
{
  if ( nbounds != 2 * length )
    {
      gen_bounds(length, coord_vals, cell_bounds);
      if ( x )
        {
          gridDefNvertex(gridID, 2);
          gridDefXbounds(gridID, cell_bounds);
        }
      else
        gridDefYbounds(gridID, cell_bounds);
    }
}

static double lonbnds_mids_trans_check(double value1, double value2)
{
  if ( abs(value1 - value2) < 180.0 )
    return (value1 + value2) * 0.5;
  else 
    {
      if ( value1 + value2 < 360.0 )
        return (value1 + value2 + 360.0) * 0.5;
      else
        return (value1 + value2 + 360.0) * 0.5 - 360.0;
    }
}

static double lonbnds_bnds_trans_check(double value1, double value2)
{
  if ( abs(value1 - value2) < 180 )
    {
      if ( 2*value1 < value2 )
        return (2*value1 - value2 + 360.0);
      else if ( 2*value1 > value2 + 360.0 )
        return (2*value1 - value2 - 360.0);
      else
        return (2*value1 - value2);
    }
  else if ( value1 - value2 > 180  )
    return (2*value1 - value2 - 360.0);
  else
    return (2*value1 - value2 + 360.0);
}

static void check_and_gen_bounds_curv(int gridID, int totalsize, int xnbounds, int xlength, double *xcoord_vals, double *xcell_bounds, int ynbounds, int ylength, double *ycoord_vals, double *ycell_bounds)
{ 
  if ( xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (xcell_bounds[1] == 0.00 && xcell_bounds[2] == 0.00) || (ycell_bounds[1] == 0.00 && ycell_bounds[2] == 0.00) )
    {
      double halflons[xlength+1][ylength];
      double halflats[xlength][ylength+1];
      double halflonsOnhalflats[xlength+1][ylength+1];
      double halflatsOnhalflons[xlength+1][ylength+1];

/**/
/*************Half-lons with 360-0 transmission check**************/
/**/
      for ( int j = 0; j < ylength; j++ )
        {
          for ( int i = 1; i < xlength; i++ )
            halflons[i][j] = lonbnds_mids_trans_check(xcoord_vals[i-1+j*xlength], xcoord_vals[i+j*xlength]);
/*left and right boundary: */
          halflons[0][j]       = lonbnds_bnds_trans_check(xcoord_vals[j*xlength], halflons[1][j]);
          halflons[xlength][j] = lonbnds_bnds_trans_check(xcoord_vals[j*xlength-1], halflons[xlength-1][j]);
        }
/**/
/*************Half-lats **************/
/**/
      for ( int i = 0; i < xlength; i++ )
        {
          for ( int j = 1; j < ylength; j++ )
            halflats[i][j] = (ycoord_vals[i+(j-1)*xlength] + ycoord_vals[i+j*xlength]) * 0.5;
/*upper and lower boundary: */
          halflats[i][0]       = 2*ycoord_vals[i] - halflats[i][1];
          halflats[i][ylength] = 2*ycoord_vals[i+(ylength-1)*xlength] - halflats[i][ylength-1];
        }
/**/
/****************Half-lons-on-half-lats with 0-360 transmission check**********/
/****************Half-lats-on-half-lons                              **********/
/**/

      for ( int i = 1; i < xlength; i++ )
        {
          for ( int j = 1; j < ylength; j++ )
            {
              halflonsOnhalflats[i][j] = lonbnds_mids_trans_check(halflons[i][j-1], halflons[i][j]);
              halflatsOnhalflons[i][j] = ( halflats[i-1][j] + halflats[i][j] ) * 0.5;
            }
/*upper and lower boundary: */
          halflonsOnhalflats[i][0]       = lonbnds_bnds_trans_check(halflons[i][0], halflonsOnhalflats[i][1]);
          halflonsOnhalflats[i][ylength] = lonbnds_bnds_trans_check(halflons[i][ylength-1], halflonsOnhalflats[i][ylength-1]);
          halflatsOnhalflons[i][0]       = ( halflats[i-1][0] + halflats[i][0] ) * 0.5;
          halflatsOnhalflons[i][ylength] = ( halflats[i-1][ylength] + halflats[i][ylength] ) * 0.5;
        }      

/*left and right boundary: */
      for ( int j = 1; j < ylength; j++ )
        {
          halflonsOnhalflats[0][j]       = lonbnds_mids_trans_check(halflons[0][j-1], halflons[0][j]);
          halflonsOnhalflats[xlength][j] = lonbnds_mids_trans_check(halflons[xlength][j-1], halflons[xlength][j]);

          halflatsOnhalflons[0][j]       = 2*halflats[0][j] - halflatsOnhalflons[1][j];
          halflatsOnhalflons[xlength][j] = 2*halflats[xlength-1][j] - halflatsOnhalflons[xlength-1][j];
        }
      halflatsOnhalflons[0][0]             = 2*halflats[0][0] - halflatsOnhalflons[1][0];
      halflatsOnhalflons[0][ylength]       = 2*halflats[0][ylength] - halflatsOnhalflons[1][ylength];
      halflatsOnhalflons[xlength][0]       = 2*halflats[xlength-1][0] - halflatsOnhalflons[xlength-1][0];
      halflatsOnhalflons[xlength][ylength] = 2*halflats[xlength-1][ylength] - halflatsOnhalflons[xlength-1][ylength];

      halflonsOnhalflats[0][0]             = lonbnds_bnds_trans_check(halflons[0][0], halflonsOnhalflats[0][1]);
      halflonsOnhalflats[0][ylength]       = lonbnds_bnds_trans_check(halflons[0][ylength-1], halflonsOnhalflats[0][ylength-1]);
      halflonsOnhalflats[xlength][0]       = lonbnds_bnds_trans_check(halflons[xlength][0], halflonsOnhalflats[xlength][1]);
      halflonsOnhalflats[xlength][ylength] = lonbnds_bnds_trans_check(halflons[xlength][ylength-1], halflonsOnhalflats[xlength-1][ylength]);

      for ( int i = 0; i < xlength; i++ )
        for ( int j = 0; j < ylength; j++ )
          {
            xcell_bounds[4*(j*xlength+i)]   = halflonsOnhalflats[i][j+1];
            xcell_bounds[4*(j*xlength+i)+1] = halflonsOnhalflats[i][j];
            xcell_bounds[4*(j*xlength+i)+2] = halflonsOnhalflats[i+1][j];
            xcell_bounds[4*(j*xlength+i)+3] = halflonsOnhalflats[i+1][j+1];
            ycell_bounds[4*(j*xlength+i)]   = halflatsOnhalflons[i][j+1];
            ycell_bounds[4*(j*xlength+i)+1] = halflatsOnhalflons[i][j];
            ycell_bounds[4*(j*xlength+i)+2] = halflatsOnhalflons[i+1][j];
            ycell_bounds[4*(j*xlength+i)+3] = halflatsOnhalflons[i+1][j+1];
          }
      gridDefNvertex(gridID, 4);
      gridDefXbounds(gridID, xcell_bounds);
      gridDefYbounds(gridID, ycell_bounds);
    }
}
/*

  if ( xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (xcell_bounds[1] == 0.00 && xcell_bounds[2] == 0.00) || (ycell_bounds[1] == 0.00 && ycell_bounds[2] == 0.00) )
    {
      for ( int j = 1; j < ylength-1; j++ )
        for ( int i = 1; i < xlength-1; i++ )
          {
            double *star[9] =
 { &xcoord_vals[(j-1)*xlength+(i-1)], &xcoord_vals[j*xlength+(i-1)], &xcoord_vals[(j+1)*xlength+(i-1)],
   &xcoord_vals[(j-1)*xlength+i],     &xcoord_vals[j*xlength+i],     &xcoord_vals[(j+1)*xlength+i],
   &xcoord_vals[(j-1)*xlength+(i+1)], &xcoord_vals[j*xlength+(i+1)], &xcoord_vals[(j+1)*xlength+i+1] };
            double max = 0, min = 0;
            for ( int k = 0; k < 9; k++ )
              {
                max = (max < *star[k]) ? *star[k] : max;
                min = (min > *star[k]) ? *star[k] : min;
              }
            if ( ( max - min ) > 270 )
              {
                if ( *star[4] < 90 )
                  for ( int l = 0; l < 9; l++)
                    *star[l] = (*star[l] > 270 ) ? *star[l] - 360.0 : *star[l];
                else if ( *star[4] > 270 )
                  for ( int l = 0; l < 9; l++)
                    *star[l] = (*star[l] < 90 ) ? *star[l] + 360.0 : *star[l];
              }

            ycell_bounds[4*(j*xlength+i)]   = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+1] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+2] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+3] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)]   = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+1] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+2] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+3] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
            for ( int m = 0; m < 4; m++ )
              {
                if ( xcell_bounds[4*(j*xlength+i)+m] > 360 ) 
                  xcell_bounds[4*(j*xlength+i)+m] = xcell_bounds[4*(j*xlength+i)+m] - 360.0;
                else if ( xcell_bounds[4*(j*xlength+i)+m] < 0 ) 
                  xcell_bounds[4*(j*xlength+i)+m] = xcell_bounds[4*(j*xlength+i)+m] + 360.0;
              }
            for ( int l = 0; l < 9; l++)
              {
                *star[l] = (*star[l] > 360 ) ? *star[l] - 360.0 : *star[l];
                *star[l] = (*star[l] < 0 )   ? *star[l] + 360.0 : *star[l];
              }
          }

      ycell_bounds[0] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      ycell_bounds[1] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[2] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[3] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      double xcyclicKorr;
      if ( xcyclicKorr =  xcoord_vals[0] - xcoord_vals[xlength] > 270 )
        if ( xcoord_vals[0] < 90 )
           xcyclicKorr -= 
      xcell_bounds[0] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[1] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[2] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[3] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;


      ycell_bounds[4*(xlength-1)]   = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+1] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+2] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+3] = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      xcell_bounds[4*(xlength-1)]   = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+1] = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+2] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+3] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;


      ycell_bounds[4*(totalsize-xlength)]   =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+1] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+2] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+3] =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)]   =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+1] =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+2] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+3] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;


      ycell_bounds[4*totalsize-4] =    ycoord_vals[totalsize-1] + ( ycoord_vals[totalsize-1] - ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-3] = (  ycoord_vals[totalsize-1] + ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-2] = (  ycoord_vals[totalsize-1] + ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-1] =    ycoord_vals[totalsize-1] + ( ycoord_vals[totalsize-1] - ycoord_vals[totalsize-1-xlength] ) * 0.5;
      xcell_bounds[4*totalsize-4] = (  xcoord_vals[totalsize-1] + xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-3] = (  xcoord_vals[totalsize-1] + xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-2] =    xcoord_vals[totalsize-1] + ( xcoord_vals[totalsize-1] - xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-1] =    xcoord_vals[totalsize-1] + ( xcoord_vals[totalsize-1] - xcoord_vals[totalsize-2] ) * 0.5;

      for ( int i = 1; i < xlength-1; i++)
        {

          ycell_bounds[4*i]   = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i-1] + ycoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*i+1] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+2] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+3] = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i+1] + ycoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i]   = ( ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i+1] =   ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5;
          xcell_bounds[4*i+2] =   ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5;
          xcell_bounds[4*i+3] = ( ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;


          ycell_bounds[4*(totalsize-xlength+i)]   =     ycoord_vals[totalsize-xlength+i] + ( ycoord_vals[totalsize-xlength+i] - ycoord_vals[totalsize-2*xlength+i] ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+1] = ( ( ycoord_vals[totalsize-xlength+i] + ycoord_vals[totalsize-2*xlength+i] ) * 0.5 + ( ycoord_vals[totalsize-xlength+i-1] + ycoord_vals[totalsize-2*xlength+i-1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+2] = ( ( ycoord_vals[totalsize-xlength+i] + ycoord_vals[totalsize-2*xlength+i] ) * 0.5 + ( ycoord_vals[totalsize-xlength+i+1] + ycoord_vals[totalsize-2*xlength+i+1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+3] =     ycoord_vals[totalsize-xlength+i] + ( ycoord_vals[totalsize-xlength+i] - ycoord_vals[totalsize-2*xlength+i] ) * 0.5;

          xcell_bounds[4*(totalsize-xlength+i)]   = ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i-1] ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+1] = ( ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i-1] ) * 0.5 + ( xcoord_vals[totalsize-2*xlength+i] + xcoord_vals[totalsize-2*xlength+i-1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+2] = ( ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i+1] ) * 0.5 + ( xcoord_vals[totalsize-2*xlength+i] + xcoord_vals[totalsize-2*xlength+i+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+3] = ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i+1] ) * 0.5;
        }

     for ( int j = 1; j < ylength-1; j++)
        {

          ycell_bounds[4*j*xlength]   = (   ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+1] = (   ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+2] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*j*xlength+3] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength]   =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5;
          xcell_bounds[4*j*xlength+1] =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5; 
          xcell_bounds[4*j*xlength+2] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j-1)*xlength] + xcoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength+3] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j+1)*xlength] + xcoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;


          ycell_bounds[4*(j+1)*xlength-4] =  ( ( ycoord_vals[(j+1)*xlength-1] + ycoord_vals[(j+2)*xlength-1] ) * 0.5 + ( ycoord_vals[(j+1)*xlength-2] + ycoord_vals[(j+2)*xlength-2] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-3] =  ( ( ycoord_vals[(j+1)*xlength-1] + ycoord_vals[j*xlength-1] )     * 0.5 + ( ycoord_vals[(j+1)*xlength-2] + ycoord_vals[j*xlength-2] )     * 0.5 ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-2] =  (   ycoord_vals[(j+1)*xlength-1] + ycoord_vals[j*xlength-1] ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-1] =  (   ycoord_vals[(j+1)*xlength-1] + ycoord_vals[(j+2)*xlength-1] ) * 0.5;

          xcell_bounds[4*(j+1)*xlength-4] =  ( ( xcoord_vals[(j+1)*xlength-1] + xcoord_vals[(j+1)*xlength-2] ) * 0.5 + ( xcoord_vals[(j+2)*xlength-1] + xcoord_vals[(j+2)*xlength-2] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-3] =  ( ( xcoord_vals[(j+1)*xlength-1] + xcoord_vals[(j+1)*xlength-2] ) * 0.5 + ( xcoord_vals[j*xlength-1] + xcoord_vals[j*xlength-2] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-2] =      xcoord_vals[(j+1)*xlength-1] + ( xcoord_vals[(j+1)*xlength-1] - xcoord_vals[(j+1)*xlength-2] ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-1] =      xcoord_vals[(j+1)*xlength-1] + ( xcoord_vals[(j+1)*xlength-1] - xcoord_vals[(j+1)*xlength-2] ) * 0.5;
        }
      gridDefNvertex(gridID, 4);
      gridDefXbounds(gridID, xcell_bounds);
      gridDefYbounds(gridID, ycell_bounds);
    } 
} */

/*
static void select_and_register_character_dimension(char *grid_file, int *axis_ids)
{
  char *ifile = cdoStreamName(0)->args;
  if ( ifile[0] == '-' )
    cdoAbort("Cdo cmor cannot register a character dimension when several cdo operators are piped.");
  if ( strcmp(grid_file, "") == 0 )  
    register_character_dimension(axis_ids, ifile);
  else
    register_character_dimension(axis_ids, grid_file);
}
*/
static void register_lon_axis(int gridID, int xlength, int *axis_ids)
{
  double *xcoord_vals = (double *) Malloc(xlength * sizeof(double));
  if ( gridInqXvals(gridID, xcoord_vals) == 0 )
    Free(xcoord_vals);
  else
    {
      double *xcell_bounds = (double *) Malloc(2 * xlength * sizeof(double));
      int xnbounds = gridInqXbounds(gridID, xcell_bounds);
      check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
      int cmf = cmor_axis(new_axis_id(axis_ids),    (char *) "longitude",    (char *) "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);
      if ( cmf != 0 )
        cdoAbort("Function cmor_axis failed!");
      if ( xcell_bounds ) Free(xcell_bounds);
      if ( xcoord_vals ) Free(xcoord_vals);
    }
}

static void register_lat_axis(int gridID, int ylength, int *axis_ids)
{
  double *ycoord_vals = (double *) Malloc(ylength * sizeof(double));
  if ( gridInqYvals(gridID, ycoord_vals) == 0 )
    Free(ycoord_vals);
  else
    {
      double *ycell_bounds = (double *) Malloc(2 * ylength * sizeof(double));
      int ynbounds = gridInqYbounds(gridID, ycell_bounds);
      check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
      int cmf = cmor_axis(new_axis_id(axis_ids),    (char *) "latitude",    (char *) "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
      if ( cmf != 0 )
        cdoAbort("Function cmor_axis failed!");
      if ( ycell_bounds ) Free(ycell_bounds);
      if ( ycoord_vals ) Free(ycoord_vals);  
    }
}

static void register_char_axis(int numchar, char **charvals, int *axis_ids, char *chardim)
{
  int maxlen = get_strmaxlen(charvals, numchar);
  void *charcmor = (void *) Malloc ( numchar * maxlen * sizeof(char));
  sprintf((char *)charcmor, "%.*s", strlen(charvals[0]), charvals[0]);
  char blanks[maxlen];
  for ( int i = 0; i < maxlen; i++)
    blanks[i] = ' ';
  sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[0]), blanks);
  for ( int i = 1; i < numchar; i++ )
    {
      sprintf((char *)charcmor, "%s%s", (char *)charcmor, charvals[i]);
      sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[i]), blanks);   
    }
  int cmf = cmor_axis(new_axis_id(axis_ids), chardim, (char *) "", numchar, (void *)charcmor, 'c',  NULL, maxlen, NULL); 
  if ( cmf != 0 )
    cdoAbort("Function cmor_axis failed!");
  Free(charcmor);
}

static void register_projection(int *grid_ids, int projID, double *ycoord_vals, double *xcoord_vals, double *ycell_bounds, double *xcell_bounds, int xlength, int ylength)
{
              int cmf = 0;
              int pxnbounds;
              int pynbounds;
              int pylength = gridInqYsize(projID);
              int pxlength = gridInqXsize(projID);
              double *pxcoord_vals = (double *) Malloc(pxlength * sizeof(double));
              double *pycoord_vals = (double *) Malloc(pylength * sizeof(double));
              double *pxcell_bounds = (double *) Malloc(2 * pxlength * sizeof(double));
              double *pycell_bounds = (double *) Malloc(2 * pylength * sizeof(double));
              inquire_vals_and_bounds(projID, &pxnbounds, &pynbounds, pxcoord_vals, pycoord_vals, pxcell_bounds, pycell_bounds);
              check_and_gen_bounds(projID, pxnbounds, pxlength, pxcoord_vals, pxcell_bounds, 1);
              check_and_gen_bounds(projID, pynbounds, pylength, pycoord_vals, pycell_bounds, 0);

              int projtype = gridInqProjType(projID);

              char p_rll_cmor[CMOR_MAX_STRING];
              int l_p_rll = strlen("grid_north_pole_longitude")+1;
              memcpy(p_rll_cmor, "grid_north_pole_latitude\0 grid_north_pole_longitude\0north_pole_grid_longitude\0", 3*l_p_rll);

              char u_rll_cmor[CMOR_MAX_STRING];
              int l_u_rll = strlen("degrees_north")+1;
              memcpy(u_rll_cmor, "degrees_north\0degrees_east\0 degrees_east\0 ", 3*l_u_rll);


              char p_lcc_cmor[CMOR_MAX_STRING];
              int l_p_lcc = strlen("longitude_of_central_meridian")+1;
              memcpy(p_lcc_cmor, "standard_parallel1\0           longitude_of_central_meridian\0latitude_of_projection_origin\0standard_parallel2\0           ", 4*l_p_lcc);


              char u_lcc_cmor[CMOR_MAX_STRING];
              int l_u_lcc = 6;
              memcpy(u_lcc_cmor, "      \0      \0      \0      \0", 4*l_u_lcc);

              const char *p_rll[] = {"grid_north_pole_latitude",
                               "grid_north_pole_longitude",
                               "north_pole_grid_longitude", NULL};

              const char *p_lcc[] = {"standard_parallel1",
                               "longitude_of_central_meridian",
                               "latitude_of_projection_origin",
                               "standard_parallel2", NULL};

              double *parameter_values = NULL;

              char mapping[CDI_MAX_NAME]; mapping[0] = 0;
              cdiGridInqKeyStr(projID, CDI_KEY_MAPPING, CDI_MAX_NAME, mapping);

              int atttype, attlen;
              char attname[CDI_MAX_NAME];

              int natts;
              cdiInqNatts(projID, CDI_GLOBAL, &natts);

              int p_len;
              switch ( projtype )
                {
                case CDI_PROJ_RLL: 
                 p_len = (sizeof(p_rll) / sizeof(p_rll[0]))-1; break;
                case CDI_PROJ_LAEA: cdoAbort("In grid registration:\n          This grid projection is not yet enabled."); break;
                case CDI_PROJ_LCC:
                  p_len = (sizeof(p_lcc) / sizeof(p_lcc[0]))-1; break;
                case CDI_PROJ_SINU: cdoAbort("In grid registration:\n          This grid projection is not yet enabled."); break;
                }
              if ( natts != p_len )
                cdoWarning("In grid registration:\n          Number of required grid mapping attributes '%d' differs from the number of given grid mapping attributes '%d'.\n          Note that all required mapping attributes are set to 0.0 by default in case they are not given.", p_len, natts);
 
              parameter_values = (double *) Malloc(p_len * sizeof(double));
              for ( int i = 0; i < p_len; i++ )
                parameter_values[i] = 0.0;

              for ( int iatt = 0; iatt < natts; ++iatt )
                {
                  cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
                  if ( atttype == DATATYPE_FLT32 || atttype == DATATYPE_FLT64 )
                    {
                      if ( attlen > 1 )
                        cdoAbort("In grid registration:\n          Dont know what to do with grid mapping attribute '%s'.", attname);
                      double attflt[attlen];
                      cdiInqAttFlt(projID, CDI_GLOBAL, attname, attlen, attflt);
                      int i = 0;
                      for ( i = 0; i < p_len; i++ )
                        {
                          if ( projtype == CDI_PROJ_RLL )
                            if ( strcmp(attname, p_rll[i]) == 0 )
                                {
                                  parameter_values[i] = attflt[0];
                                  break;
                                }
                          else if ( projtype == CDI_PROJ_LCC ) 
                            if ( strcmp(attname, p_lcc[i]) == 0 )
                                {
                                  parameter_values[i] = attflt[0];
                                  break;
                                }

                        }
                      if ( i == p_len )
                        cdoWarning("In grid registration:\n          grid mapping attribute '%s' is neglected.", attname);
                    }
                  else if ( atttype  == DATATYPE_TXT )
                    {
                      char atttxt[attlen];
                      cdiInqAttTxt(projID, CDI_GLOBAL, attname, attlen, atttxt);
                    }
                }

              int grid_axis[2];
              if ( projtype == CDI_PROJ_RLL )
                {
                  cmf = cmor_axis(&grid_axis[0],    (char *) "grid_latitude",   (char *) "degrees_north",    pylength,    (void *)pycoord_vals,    'd',    0, 0,   NULL);
                  cmf = cmor_axis(&grid_axis[1],    (char *)"grid_longitude",   (char *) "degrees_east",    pxlength,    (void *)pxcoord_vals,    'd',    0, 0,   NULL);
                  cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    4,     (void *)ycell_bounds,    (void *)xcell_bounds);
#if ( CMOR_VERSION_MAJOR == 2 )
                  cmf = cmor_set_grid_mapping(grid_ids[0], "rotated_latitude_longitude", p_len, (char **) p_rll_cmor, l_p_rll, parameter_values, (char **)u_rll_cmor,  l_u_rll);
#elif ( CMOR_VERSION_MAJOR == 3 )
                  cmf = cmor_set_grid_mapping(grid_ids[0], (char *)"rotated_latitude_longitude", p_len, p_rll_cmor, l_p_rll, parameter_values, u_rll_cmor,  l_u_rll);
#endif
                }
              else if ( projtype == CDI_PROJ_LCC )
                {
                  double *xii = (double *) Malloc(xlength * sizeof(double));
                  double *yii = (double *) Malloc(ylength * sizeof(double));
                  for ( int i = 0; i < xlength; i++ )
                    xii[i] = (double) i;
                  for ( int i = 0; i < ylength; i++ )
                    yii[i] = (double) i;
                  cmf = cmor_axis(&grid_axis[0],  (char *)  "x",  (char *)  "m",    ylength,    (void *)yii,    'd',    0, 0,   NULL);
                  cmf = cmor_axis(&grid_axis[1],  (char *)  "y",  (char *)  "m",    xlength,    (void *)xii,    'd',    0, 0,   NULL);
                  cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    4,     (void *)ycell_bounds,    (void *)xcell_bounds);
#if ( CMOR_VERSION_MAJOR == 2 )
                  cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len,(char **)p_lcc_cmor, l_p_lcc, parameter_values, (char **)u_lcc_cmor,  l_u_lcc);
#elif ( CMOR_VERSION_MAJOR == 3 )
                  cmf = cmor_set_grid_mapping(grid_ids[0], mapping, p_len, p_lcc_cmor, l_p_lcc, parameter_values, u_lcc_cmor,  l_u_lcc);
#endif
                  Free(xii); Free(yii);
                }
              Free(parameter_values);
              Free(pxcell_bounds);
              Free(pycell_bounds);
              Free(pxcoord_vals);
              Free(pycoord_vals);
              if ( cmf != 0 )
                cdoAbort("Function cmor_axis or cmor_set_grid_mapping failed!");
}

static void register_grid(list_t *kvl, int vlistID, int varID, int *axis_ids, int *grid_ids, char *project_id)
{
  int cmf = 0;
  int gridID = vlistInqVarGrid(vlistID, varID);

  char *grid_file = kv_get_a_val(kvl, "gi", NULL);

  char *chardimatt = kv_get_a_val(kvl, "ca", NULL);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");

  if ( grid_file )
    {
      change_grid(grid_file, gridID, vlistID);
      gridID = vlistInqVarGrid(vlistID, varID);
    }

  int type = gridInqType(gridID);
  int projID = gridInqProj(gridID);
  int ylength = gridInqYsize(gridID);
  int xlength = gridInqXsize(gridID);
  int totalsize = gridInqSize(gridID);

  double *xcoord_vals;
  double *ycoord_vals;
  double *xcell_bounds;
  double *ycell_bounds;
  double *x2cell_bounds;
  double *y2cell_bounds;
  int xnbounds;
  int ynbounds;

  if ( totalsize > 1 )
    {
      if ( type == GRID_GAUSSIAN || type == GRID_LONLAT )
        {
          grid_ids[0] = 0;
          xcoord_vals = (double *) Malloc(xlength * sizeof(double));
          ycoord_vals = (double *) Malloc(ylength * sizeof(double));
          xcell_bounds = (double *) Malloc(2 * xlength * sizeof(double));
          ycell_bounds = (double *) Malloc(2 * ylength * sizeof(double));
          inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);

          check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
          check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);

          cmf = cmor_axis(new_axis_id(axis_ids),  (char *)  "latitude",   (char *) "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
          cmf = cmor_axis(new_axis_id(axis_ids),  (char *)  "longitude",  (char *) "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);

          Free(xcell_bounds);
          Free(ycell_bounds);
          Free(xcoord_vals);
          Free(ycoord_vals);
        }
      else if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED )
            {
              xcoord_vals = (double *) Malloc(totalsize * sizeof(double));
              ycoord_vals = (double *) Malloc(totalsize * sizeof(double));
        /* maximal 4 gridbounds per gridcell permitted */
              xcell_bounds = (double *) Malloc(4 * totalsize * sizeof(double));
              ycell_bounds = (double *) Malloc(4 * totalsize * sizeof(double));
              inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
       /* In a projection, this is done by setting mapping parameter */
              move_lons(xcoord_vals, xcell_bounds, totalsize, 4 * totalsize, xnbounds);   
              get_cmor_table(kvl, project_id);
              int grid_axis[2];
              check_and_gen_bounds_curv(gridID, totalsize, xnbounds, xlength, xcoord_vals, xcell_bounds, ynbounds, ylength, ycoord_vals, ycell_bounds);
              if ( type == GRID_CURVILINEAR && projID == CDI_UNDEFID )
                {
                  double *xncoord_vals;
                  double *yncoord_vals;
                  xncoord_vals = (double *) Malloc(xlength * sizeof(double));
                  yncoord_vals = (double *) Malloc(ylength * sizeof(double)); 
                  for (int j=0; j<ylength; j++) 
                    yncoord_vals[j]= (double) j;
                  for (int j=0; j<xlength; j++)
                    xncoord_vals[j]= (double) j;
                  cmf = cmor_axis(&grid_axis[0],(char *) "j_index",   (char *) "1",    ylength,    (void *)yncoord_vals,
            'd', 0, 0, NULL);
                  cmf = cmor_axis(&grid_axis[1],(char *) "i_index",   (char *) "1",    xlength,    (void *)xncoord_vals,    'd', 0, 0, NULL);
                  cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    4,     (void *)ycell_bounds,    (void *)xcell_bounds);
                  Free(xncoord_vals);
                  Free(yncoord_vals);
                  Free(xcoord_vals);
                  Free(ycoord_vals);
                  Free(xcell_bounds);
                  Free(ycell_bounds);
                }
              /*else
                { 
                  cmf = cmor_axis(&grid_axis[0],    "grid_longitude",   "degrees",    xlength,    (void *)xcoord_vals,    'd', 0, 0, NULL);
                  cmf = cmor_axis(&grid_axis[1],    "grid_latitude",    "degrees",    ylength,    (void *)ycoord_vals,    'd', 0, 0, NULL);
                  cmf = cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    2,     (void *)ycell_bounds,    (void *)xcell_bounds); 
                }*/
            }
      else if ( type == GRID_GENERIC && ( strcmp(chardim, "oline") == 0 || strcmp(chardim, "basin") == 0 ))
            {
              if ( cdoVerbose )
                cdoPrint("Start to define a character axis '%s' instead of a grid axis'.", chardim);
              grid_ids[0] = 0;
              int numchar = 0;
              char *charvalstring = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
              sprintf(charvalstring, "char_axis_%s", chardim);
              char **charvals = kv_get_vals(kvl, charvalstring, &numchar);
              Free(charvalstring);
              if ( ( xlength > 0 && xlength != numchar ) && ( ylength > 0 && ylength != numchar ) )
                cdoAbort("In registration of a character coordinate as substitution for a horizontal axis:\n          You configured a character coordinate '%s' with '%d' string values but you also registered a grid with '%d' numerical values on X axis and '%d' numerical values on Y axis. One axis must match the number of string values.", chardim, numchar, xlength, ylength);
              if ( !charvals )
                cdoAbort("In registration of a character coordinate as substitution for a horizontal axis:\n          You configured a character coordinate '%s' but no values are found! Configure values via attribute 'char_dim_vals'!", chardim);
              if ( charvals && ( xlength == numchar || xlength == 0 ) )
                {
                  register_char_axis(numchar, charvals, axis_ids, chardim);
                  if ( ylength > 0 )
                    register_lat_axis(gridID, ylength, axis_ids);
                }
              else
                {
                  register_lon_axis(gridID, xlength, axis_ids);
                  register_char_axis(numchar, charvals, axis_ids, chardim);
                }
              if ( cdoVerbose )
                cdoPrint("Successfully defined a character axis '%s' instead of a grid axis.", chardim);
            }
        /*
              grid_ids[0] = 0;
              xcoord_vals = Malloc(xlength * sizeof(double));
              gridInqXvals(gridID, xcoord_vals);
              ycoord_vals = Malloc(ylength * sizeof(double));
              gridInqYvals(gridID, ycoord_vals);
              char yname[CDI_MAX_NAME];
              cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, yname);
              if ( strcmp(yname, "basin") != 0 )
                {
                  invert_ygriddes(kvl, vlistID, &gridID, ylength, ycoord_vals, ycell_bounds, &ynbounds);        
                  ycell_bounds = Malloc(2 * ylength * sizeof(double));        
                  ynbounds = gridInqYbounds(gridID, ycell_bounds);
                  check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
                  cmf = cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
                  Free(ycell_bounds);
                }
              else
                select_and_register_character_dimension(grid_file, axis_ids);

              char xname[CDI_MAX_NAME];
              cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, xname);
              if ( strcmp(xname, "basin") != 0 )
                {
                  xcell_bounds = Malloc(2 * xlength * sizeof(double));
                  xnbounds = gridInqXbounds(gridID, xcell_bounds);
                  check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
                  cmf = cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);
                  Free(xcell_bounds);
                }
              else
                select_and_register_character_dimension(grid_file, axis_ids);
              Free(xcoord_vals);
              Free(ycoord_vals);
            */
      else if ( type == GRID_CHARXY )
            {
              grid_ids[0] = 0;
              char *xname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
              char *yname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
              gridInqXname(gridID, xname);
              gridInqYname(gridID, yname);
              char *xdimname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
              char *ydimname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
              cdiGridInqKeyStr(gridID, 902, CDI_MAX_NAME, xdimname);
              cdiGridInqKeyStr(gridID, 912, CDI_MAX_NAME, ydimname);
              if ( strcmp(xdimname, "line") == 0 )
                strcpy(xdimname, "oline");
              int dimstrlen;   
              if ( dimstrlen = gridInqXIsc(gridID) )
                {
                  char **xchars = (char **)Malloc( (xlength+1) * sizeof(char *));
                  for ( int i = 0; i < xlength; i++ )
                    xchars[i] = (char *)Malloc( (dimstrlen+1) * sizeof(char));
                  gridInqXCvals(gridID, xchars);
                  for ( int j = 0; j < xlength; j++ )
                    xchars[j][dimstrlen] = 0;
                  xchars[xlength] = NULL;
                  register_char_axis(xlength, xchars, axis_ids, xdimname);
                  free_array(xchars);
                }
              else if ( xlength)
                register_lon_axis(gridID, xlength, axis_ids);

              if ( dimstrlen = gridInqYIsc(gridID) )
                {
                  char **ychars = (char **) Malloc( (ylength + 1) * sizeof(char));
                  for ( int i = 0; i < ylength; i++ )
                    ychars[i] = (char *)Malloc( (dimstrlen +1) * sizeof(char));
                  gridInqYCvals(gridID, ychars);
                  for ( int j = 0; j < ylength; j++ )
                    ychars[j][dimstrlen] = 0;
                  ychars[ylength] = NULL;
                  register_char_axis(ylength, ychars, axis_ids, ydimname);
                  free_array(ychars);
                }
              else if ( ylength )
                register_lat_axis(gridID, ylength, axis_ids);
              Free(xname); Free(yname); Free(xdimname); Free(ydimname);
            }
      else if ( type == GRID_PROJECTION )
            {
              cdoAbort("In grid registration:\n          For a 'rotated_lat_lon' projection, both grids, the unprojected lat/lon and the projected rlat/rlon are required.");            
            }
      else
            {
              grid_ids[0] = 0;
              cdoWarning("Registration of a grid is skipped. Either the grid type is unknown or a registration is not necessary.");
            }

      if ( projID != CDI_UNDEFID )
            {
              register_projection(grid_ids, projID, ycoord_vals, xcoord_vals, ycell_bounds, xcell_bounds, xlength, ylength);
              Free(xcoord_vals);
              Free(ycoord_vals);
              Free(xcell_bounds);
              Free(ycell_bounds);
            }
        
    }
  else
    grid_ids[0] = 0;
  Free(chardim);
  if ( cmf != 0 )
    cdoAbort("Function cmor_axis failed!");
}

static void register_variable(list_t *kvl, int vlistID, int varID, int *axis_ids,
                              struct mapping *var, int *grid_ids, char *name)
{
  int cmf = 0;
  if ( cdoVerbose )
    cdoPrint("Start to retrieve 'positive' and 'units'.");
  char *positive = get_txtatt(vlistID, varID, "positive");
  char *origname = get_txtatt(vlistID, varID, "original_name");
  char *history = get_txtatt(vlistID, varID, "history");
  char *varcom = get_txtatt(vlistID, varID, "variable_comment");
  char *units = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
  vlistInqVarUnits(vlistID, varID, units);
  char *attunits = kv_get_a_val(kvl, "u", NULL);
  char *attp = kv_get_a_val(kvl, "p", NULL);
  char *attorigname = kv_get_a_val(kvl, "original_name", NULL);
  char *attvarcom = kv_get_a_val(kvl, "vc", NULL);
  check_compare_set(&positive, attp, "positive", "");
  if ( strcmp(positive, " ") == 0 )
    strcpy(positive, "");
  check_compare_set(&units, attunits, "units", NULL);
  check_compare_set(&origname, attorigname, "original_name", "");
  if ( strcmp(origname, "") == 0 || strstr(origname, "var") )
    {
      Free(origname);
      origname = NULL;
    }
  check_compare_set(&varcom, attvarcom, "variable_comment", "");
  if ( strcmp(varcom, "") == 0 )
    {
      Free(varcom);
      varcom = NULL;
    }
  if ( cdoVerbose )
    cdoPrint("Successfully retrieved 'positive': '%s' and 'units' : '%s'.", positive, units);
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = vlistGridsizeMax(vlistID);
  int zsize = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->cdi_varID = varID;
  var->help_var = 0;
  if ( !var->data )
    {
      var->charvars = 0;
      if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
        {
          var->datatype = 'f';
          *(float *) missing_value = vlistInqVarMissval(vlistID, varID);
          var->data = Malloc(gridsize * zsize * sizeof(float));
        }
      else
        {
          var->datatype = 'd';
          *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
          var->data = Malloc(gridsize * zsize * sizeof(double));
        }
    }
  else
    *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
   
  if ( cdoVerbose )
    cdoPrint("Start to call cmor_variable.");
  if ( grid_ids[0] != 0 )
    {
      int *tmp_id = new_axis_id(axis_ids);
      *tmp_id = grid_ids[0];
      cmf = cmor_variable(&var->cmor_varID,
            name,units,(count_axis_ids(axis_ids)), axis_ids, var->datatype,
            (void *) missing_value, &tolerance, positive,
                        origname,
                        history,
                        kv_get_a_val(kvl, "vc", NULL));
    }
  else
    {
      cmf = cmor_variable(&var->cmor_varID,
           name, units, count_axis_ids(axis_ids),  axis_ids,   var->datatype,
          (void *) missing_value, &tolerance, positive,
                        origname,
                        history,
                        kv_get_a_val(kvl, "vc", NULL));
    }
  if ( cmf != 0 )
    cdoAbort("Function cmor_variable failed!");
  if ( cdoVerbose )
    cdoPrint("Successfully called cmor_variable.");
  if (positive) Free(positive); 
  if (origname) Free(origname); 
  if (history) Free(history); 
  if (units) Free(units);
}

static void register_all_dimensions(list_t *kvl, int streamID,
                             struct mapping vars[], int table_id, char *project_id, int miptab_freq, int *time_axis)
{
  int cmf = 0;
  int vlistID = pstreamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  char *time_units = kv_get_a_val(kvl, "required_time_units", NULL);

  if ( cdoVerbose )
    cdoPrint("7. Start to retrieve requested variables.");

  int numvals = 0;
  char **cmor_names = kv_get_vals(kvl, "cn", &numvals);

/* Cmdlinemapping: */
  char *mapname, *mapcode;
  if ( !kv_get_a_val(kvl, "mt", NULL) && numvals )
    {
      if ( mapname = kv_get_a_val(kvl, "n", NULL) )
        change_name_via_name(vlistID, mapname, cmor_names[0]);
      else if ( mapcode = kv_get_a_val(kvl, "c", NULL) )
        change_name_via_code(vlistID, mapcode, cmor_names[0]);
    }

  if ( cmor_names == NULL && vlistNvars(vlistID) > 1 )
    cdoPrint("In registration of all dimensions for the variables:\n          You have not requested any specific variable but there are several in infile! Notice that if attributes e.g. units are configured via cmdline, they will be used for every variable!");
  if ( cdoVerbose )
    cdoPrint("7. Successfully retrieved requested variables");
  int foundName = 0;
  int ps_required = 0;
  int ps_in_file = 0;
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      char name[CDI_MAX_NAME];
      vlistInqVarName(vlistID, varID, name);
      if ( !cmor_names || in_list(cmor_names, name, numvals) )
        {
          struct mapping *var = map_var(varID, vars);
          if ( !var )
            var = new_var_mapping(vars);
          int axis_ids[CMOR_MAX_AXES];
          axis_ids[0] = CMOR_UNDEFID;
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          if ( cdoVerbose )
            cdoPrint("8. Start to define variable with ID: '%d' and name: '%s'", varID, name);
          if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
            {
              cdoPrint("Since the zaxis of variable '%s' is of type HYBRID, surface pressure is required. An infile variable must have the name ps.", name);
              ps_required++;
            }
          foundName++;
          /* Time-Axis */
          if ( cdoVerbose )
            cdoPrint("8.1. Start to register a time axis");
          char cmor_time_name[CMOR_MAX_STRING]; cmor_time_name[0] = '\0';
          get_time_method(kvl, vlistID, varID, cmor_time_name, project_id, miptab_freq, time_axis);
          if ( strcmp(cmor_time_name, "none") != 0 )
            cmf = cmor_axis(new_axis_id(axis_ids),
                    cmor_time_name,
                    time_units,
                    0,NULL, 0, NULL, 0, NULL);
          if ( cdoVerbose && cmf == 0 )
            cdoPrint("8.1. Successfully handled time axis registration.");
          else if ( cmf != 0 )
            cdoAbort("Function cmor_axis failed!");
          /* Grid: */
          if ( cdoVerbose )
            cdoPrint("8.2. Start to register a grid");
          int grid_ids[CMOR_MAX_GRIDS];
          register_grid(kvl, vlistID, varID, axis_ids, grid_ids, project_id);
          cmf = cmor_set_table(table_id);
          if ( cmf != 0 )
            cdoAbort("Function cmor_set_table failed!");
          if ( cdoVerbose )
            cdoPrint("8.2. Successfully handled grid registration.");
          /* Z-Axis */
          if ( cdoVerbose )
            cdoPrint("8.3. Start to register a zaxis");
          register_z_axis(kvl, vlistID, varID, zaxisID, name, axis_ids, &var->zfactor_id, project_id, miptab_freq);
          if ( cdoVerbose )
            cdoPrint("8.3. Successfully handled zaxis registration.");
          /* Variable */
          register_variable(kvl, vlistID, varID, axis_ids, var, grid_ids, name);     
          if ( cdoVerbose )
            cdoPrint("8. Successfully defined variable with ID: '%d' and name: '%s'.", varID, name);
        }
    }
  if ( ps_required )
    {
      if ( cdoVerbose )
        cdoPrint("9. Start to find surface pressure.");
      for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
        if ( vlistInqVarCode(vlistID, varID) == 134 )
          {
            ps_in_file++;
            if ( cmor_names == NULL || in_list(cmor_names, "ps", numvals) )
              break;
            else
              {
                struct mapping *var = new_var_mapping(vars);
                size_t gridsize = vlistGridsizeMax(vlistID);
                var->cdi_varID = varID;
                var->help_var = 1;
                if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
                  {
                    var->datatype = 'f';
                    var->data = Malloc(gridsize * sizeof(float));
                  }
                else
                  {
                    var->datatype = 'd';
                    var->data = Malloc(gridsize * sizeof(double));
                  }
                break;
              }
          }
      if ( cdoVerbose )
        cdoPrint("9. Successfully registered surface pressure.");
    }
  if ( ps_required && !ps_in_file )
    cdoAbort("After registration of all dimensions for all variables:\n          No surface pressure found in infile but required for a hybrid sigma pressure z axis!");
  if ( !foundName && cmor_names )
    cdoAbort("After registration of all dimensions for all variables:\n          None of the given variables to process by attribute 'cmor_name' is found in infile.");
  if ( cdoVerbose )
    cdoPrint("Successfully registered all dimensions for %d variables successfully.", foundName);
}

static char *get_frequency(list_t *kvl, int streamID, int vlistID, int taxisID, int miptab_freq)
{
  char *frequency = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
  int ntsteps = vlistNtsteps(vlistID);
  int reccounter = 0;
  int recdummy = 0;

  switch ( miptab_freq )
    {
    case 11: strcpy(frequency, "yr"); break;
    case 2: strcpy(frequency, "yr"); break;
    case 12: strcpy(frequency, "mon"); break;
    case 3: strcpy(frequency, "mon"); break;
    case 13: strcpy(frequency, "day"); break;
    case 4: strcpy(frequency, "day"); break;
    case 14: strcpy(frequency, "6hr"); break;
    case 5: strcpy(frequency, "6hr"); break;
    case 6: strcpy(frequency, "6hr"); break;
    case 15: strcpy(frequency, "3hr"); break;
    default:
    {
      if ( cdoStreamName(0)->args[0] == '-' )
        {
            cdoAbort("No frequency could be determined from MIP-table and cdo cmor cannot check frequency of Ifile recs since you piped several cdo operators.");
/*          char *dummy;
          cdoWarning("Cdo cmor cannot check frequency of Ifile recs since you piped several cdo operators.\nIt is tried to use a configuration attribute frequency.");
          if ( !(dummy = kv_get_a_val(kvl, "frequency", NULL)) )
            cdoAbort("No attribute frequency is found.");
          else
            {
              strcpy(frequency, dummy);
              return frequency;
            }
*/
        } 
      
      int streamID2 = pstreamOpenRead(cdoStreamName(0));
      int vlistID2 = pstreamInqVlist(streamID2);
      int taxisID2 = vlistInqTaxis(vlistID2);
      if ( ntsteps < 0 )
        {
          while ( recdummy = pstreamInqTimestep(streamID2, reccounter++) );
          ntsteps = reccounter;
        }    
      ntsteps-=1;
      int fyear, lyear, fmonth, lmonth, dummyone, dummytwo;

      if ( ntsteps > 2 )
        {
          int recfirst = pstreamInqTimestep(streamID2, 0);
          cdiDecodeDate(taxisInqVdate(taxisID2), &fyear, &fmonth, &dummytwo);
          int reclast = pstreamInqTimestep(streamID2, ntsteps);    
          cdiDecodeDate(taxisInqVdate(taxisID2), &lyear, &lmonth, &dummytwo);

          double covered_years = lyear-fyear + 1.0;
          if ( DBL_IS_EQUAL(ntsteps / covered_years, 1) )
            strcpy(frequency, "yr");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 12) )
            strcpy(frequency, "mon");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366) )
            strcpy(frequency, "day");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365*4) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25*4) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366*4) )
            strcpy(frequency, "6hr");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365*8) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25*8) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366*8) )
            strcpy(frequency, "3hr");
          else 
            {
              int step_per_year = 0;
              reccounter = 0;
              if ( cdoVerbose )
                cdoPrint("Frequency is calculated by counting all timesteps in year %d\n          in order to calculate time bounds in case they are not given.", fyear, fmonth);
              while ( recdummy = pstreamInqTimestep(streamID2, reccounter++) )
                {
                  int reqyear;
                  cdiDecodeDate(taxisInqVdate(taxisID2), &reqyear, &lmonth, &dummytwo);
                  if ( reqyear == ( fyear + 1 ) )
                    break;
                  step_per_year++;
                } 
              int covered_months = lmonth-fmonth+1;
              if ( step_per_year > 366*8 )
                cdoAbort("In estimating frequency:\n          Frequency is sub-3hourly! Not yet enabled.");
              else
                {
                  if ( (double)step_per_year / (double)covered_months > 31*8 )
                    cdoAbort("Frequency is sub-3hourly! Not yet enabled.");
                  else if ( (double)step_per_year / (double)covered_months > 31*4 )
                    strcpy(frequency, "3hr");
                  else if ( (double)step_per_year / (double)covered_months > 31 )
                    strcpy(frequency, "6hr");
                  else if ( (double)step_per_year / (double)covered_months > 1 )
                    strcpy(frequency, "day");
                  else
                    strcpy(frequency, "mon");
                }
              if ( cdoVerbose )
                cdoPrint("Found %d time steps in year %d.\n          Therefore, the frequency is %s.", step_per_year, fyear, frequency);
            }
        }
      else
        {
          if ( !taxisHasBounds(taxisID2) && ntsteps > 0 )
            cdoAbort("In estimating frequency:\n          No time bounds are found in Ifile and for %d found timesteps no frequency can be computed - at least 3 timesteps are required.\n          Define time bounds before cdo cmor.", ntsteps);
          else
            cdoWarning("In frequency estimation:\n          For %d found timesteps no frequency can be computed - at least 3 timesteps are required.\n          Time bounds of the rec are used.", ntsteps);
        }
      pstreamClose(streamID2);
    }
    }
  return frequency;
}

static int get_tunitsec(int tunit)
{
  switch ( tunit )
    {
    case TUNIT_MINUTE: return 60; 
    case TUNIT_HOUR: return 3600; 
    case TUNIT_DAY: return 86400; 
    default: return 3600;
    }
}

static juldate_t get_cmor_time_val(int taxisID, juldate_t ref_date, int tunitsec, int calendar, char *frequency, int ts_id)
{
  int year, month, day;
  cdiDecodeDate(taxisInqVdate(taxisID), &year, &month, &day);
  juldate_t juldate = juldate_encode(calendar, taxisInqVdate(taxisID),
                                     taxisInqVtime(taxisID));

  if ( month == 0 || day == 0 || year == 0 )
    {
      int rdate, rtime;
      int ryear, rmonth, rday, addseconds = 0;
      juldate_decode(calendar, ref_date, &rdate, &rtime);
      cdiDecodeDate(rdate, &ryear, &rmonth, &rday);
      if ( ts_id < 2 )
        cdoWarning("In writing the data:\n          Time axis is incorrect. It is tried to calculate time values with frequency.\n          Note: These are only valid if\n           - cm=m \n           - a equally spaced monotonical time axis exist according to the frequency \n           - a correct calendar exist!");
      if ( strcmp(frequency, "yr") == 0 )
        {
          year = ryear+ts_id;
          month = 6; /* Is set to mid point by CMOR */
          day = 14; /* Is set to mid point by CMOR */
        }
      else if ( strcmp(frequency, "mon") == 0 )
        {
          year = ryear + floor(((double)(ts_id-1))/12);
          month = (ts_id % 12);
          if ( month == 0 )
            month = 12;
          day = 14; /* Is set to mid point by CMOR */
        }
      else if ( strcmp(frequency, "day") == 0 )
        {
          addseconds = ts_id * 24*60*60 + 60*60*12;
          juldate = juldate_add_seconds(addseconds, ref_date);
        }
      else if ( strcmp(frequency, "6hr") == 0 )
        {
          addseconds = ts_id * 6*60*60;
          juldate = juldate_add_seconds(addseconds, ref_date);
        }
      else if ( strcmp(frequency, "3hr") == 0 )
        {
          addseconds = ts_id * 3*60*60;
          juldate = juldate_add_seconds(addseconds, ref_date);
        }
      if ( addseconds == 0 )
        {
          int vdate = cdiEncodeDate(year, month, 1);
          int vtime = 0;
          juldate = juldate_encode(calendar, vdate, vtime);
        }
    }

  return juldate;
}

static double *get_time_bounds(list_t *kvl, int taxisID, char *frequency, juldate_t ref_date, juldate_t jtime_val, int calendar, int tunitsec, double *time_bnds, int time_axis)
{
  double time_val = juldate_to_seconds(juldate_sub(jtime_val, ref_date)) / tunitsec;
  int vdate0b, vdate1b, vtime0b, vtime1b, vdatecorr, vtimecorr;
  int year, month, day;
  int hour, min, sec;
  cdiDecodeDate(taxisInqVdate(taxisID), &year, &month, &day);
  if ( month == 0 || day == 0 )
    {
      juldate_decode(calendar, jtime_val, &vdatecorr, &vtimecorr);
      cdiDecodeDate(vdatecorr, &year, &month, &day);
    }
/***/
/* If file time axis has bounds use them, otherwise use cmor time axis deduced from miptable frequency and cell_methods or frequency itself*/
/***/

  if ( !taxisHasBounds(taxisID) || strcmp(kv_get_a_val(kvl, "tbnds_force", "n"), "y") == 0 )
    {
      vtime0b = 0;
      vtime1b = 0;
  
      if ( time_axis == 2 || time_axis == 3 )
        {
          int numdates;
          char **climyears = kv_get_vals(kvl, "climatology_interval", &numdates);
          if ( numdates != 2 )
            cdoAbort("In writing model output:\n          Could not calculate time bounds for climatology time axis because attribute 'climatology_interval' has not two values.");
          int expstartyear = atol(climyears[0]);
          int expendyear = atol(climyears[1]);
/***/
/* Climatologies */
/***/
          if ( time_axis == 2 )
            {
              vdate0b = cdiEncodeDate(expstartyear, month, 1);
              month++;
              if ( month > 12 ) { month = 1; year++; }
              vdate1b = cdiEncodeDate(expendyear, month, 1);
            }
/***/
/* Diurnal cycle */
/***/
          if ( time_axis == 3 )
            {
              vdate0b = cdiEncodeDate(expstartyear, month, day);
              cdiDecodeTime(taxisInqVtime(taxisID), & hour, &min, &sec);
              vtime0b = cdiEncodeTime(hour,0,0);

              hour++;
              if ( hour > 23 ) { hour = 0; day++; }
              vtime1b = cdiEncodeTime(hour,0,0);
              vdate1b = cdiEncodeDate(expendyear, month, day);
            }
        }
      else
        {
/***/
/* Frequency dependent: */
/***/
          if ( strcmp(frequency, "yr") == 0 )
            {
              vdate0b = cdiEncodeDate(year,   1, 1);
              vdate1b = cdiEncodeDate(year+1, 1, 1);
            }     
          else if ( strcmp(frequency, "mon") == 0 )
            {
              vdate0b = cdiEncodeDate(year, month, 1);
              month++;
              if ( month > 12 ) { month = 1; year++; }
              vdate1b = cdiEncodeDate(year, month, 1);
            }  
          else if ( strcmp(frequency, "day") == 0 )
            {
              time_bnds[0] = floor(time_val);
              time_bnds[1] = ceil(time_val);
              return time_bnds;
            }
/***/
/* Note that time_val must be correct in Infile for subdaily frequencies */
/***/  
          else if ( strcmp(frequency, "6hr") == 0 )
            {
              time_bnds[0] = time_val - 0.125;
              time_bnds[1] = time_val + 0.125;
              return time_bnds;
            }  
          else if ( strcmp(frequency, "3hr") == 0 )
            {
              time_bnds[0] = time_val - 0.0625;
              time_bnds[1] = time_val + 0.0625;
              return time_bnds;
            } 
        }
    }
  else
    {
      taxisInqVdateBounds(taxisID, &vdate0b, &vdate1b);
      taxisInqVtimeBounds(taxisID, &vtime0b, &vtime1b);
    }
  juldate_t juldate = juldate_encode(calendar, vdate0b, vtime0b);
  time_bnds[0] = juldate_to_seconds(juldate_sub(juldate, ref_date))
                / tunitsec;

  juldate = juldate_encode(calendar, vdate1b, vtime1b);
  time_bnds[1] = juldate_to_seconds(juldate_sub(juldate, ref_date))
              / tunitsec;
  return time_bnds;
}


static void read_record(int streamID, struct mapping vars[], int vlistID)
{
  int varID, levelID;
  pstreamInqRecord(streamID, &varID, &levelID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int type = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  double *buffer = (double *) Malloc(gridsize * sizeof(double));

  struct mapping *var = map_var(varID, vars);
  if ( var && var->charvars != 1 )
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int latdim = gridInqYsize(gridID);
      int levdim = zaxisInqSize(zaxisID);
      int chardim = gridsize/latdim;
      int nmiss;
      pstreamReadRecord(streamID, buffer, &nmiss);
      for ( size_t i = 0; i < gridsize; i++ )
        {
// Wrong:  (lat x basin, lev ) gridsize * levelID + i
// Wrong:  (basin x lat, lev) gridsize * levelID + i * chardim - ( int ) floor(i / latdim) * gridsize + ( int ) floor(i/latdim)
// Wrong:  (basin x lev, lat ) gridsize/latdim * levdim * ( i - ( int ) floor(i/latdim) * latdim ) + ( int ) floor(i/latdim) + gridsize/latdim * levelID;
// Wrong:  (lat x lev, basin ) latdim * levdim * ( int ) floor(i/latdim) + ( i - ( int ) floor(i/latdim) * latdim ) + levelID * latdim
// (lev x lat, basin )
          int newIndex;
          if ( levdim > 1 && type == GRID_CURVILINEAR )
            newIndex = i + gridsize*levelID;
          else if ( levdim > 1 )
            newIndex = i * levdim + levelID;
          else
            newIndex = i;
          if ( var->datatype == 'f' )
            {
              ((float *)var->data)[newIndex] = (float)buffer[i];
            }
          else
            {
              ((double *)var->data)[newIndex] = (double)buffer[i];
            }
        }
    }
  Free(buffer);
}

static void check_for_sfc_pressure(int *ps_index, struct mapping vars[], int vlistID, int timestep)
{
  int ps_required = 0;
  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++ )
    {
      if ( vlistInqVarCode(vlistID, vars[j].cdi_varID) == 134 )
        *ps_index = j;
      else if ( zaxisInqType(vlistInqVarZaxis(vlistID, vars[j].cdi_varID)) == ZAXIS_HYBRID )
        ps_required ++;
    }
  if ( *ps_index < 0 && ps_required )
    cdoAbort("In writing data with CMOR:\n          No surface pressure found for time step %d but required in Hybrid-sigma-pressure-coordinates. ", timestep);
}


static int check_append_and_size(list_t *kvl, int vlistID, char *testIn, int ifreq, int calendar)
{
  char *test = testIn;
  size_t filesize = fileSize((const char *)testIn);
  char old_start_date[CMOR_MAX_STRING];
  char old_end_date[CMOR_MAX_STRING];
  int i = 0, j = 0;
/* Get dates from chunk string */
  if ( cdoVerbose) cdoPrint("Start to retrieve dates from chunk string.");
  while ( *(test+i) != 0 )
    {
      if ( *(test+i) == '_' )
        {
          test+=(i+1);
          i = 0;
        }
      if ( *(test+i) == '-' )
        j = i;
      i++;
    }
  if ( !i || !j || *(test+j+1) == 0 || *(test+2*j) == 0 )
    {
      cdoWarning("In checking the last chunk:\n          Date from filename of the chunk cannot be read.\n          Switched to replace mode for this variable.");
      return 0;
    }

  strncpy(old_start_date, test, j);
  old_start_date[j] = 0;
  test += (j + 1);
  strncpy(old_end_date, test, j);
  old_end_date[j] = 0;

  if ( cdoVerbose) cdoPrint("Successfully retrieved start date: '%s' and end date: '%s' chunk string.\n", old_start_date, old_end_date);
/* Check frequency of chunk with frequency of file */

  if ( (j == 8 && ifreq !=3) || (ifreq == 3 && j != 8)
    || (j == 6 && ifreq !=2) || (ifreq == 2 && j != 6)
    || (j == 4 && ifreq !=1) || (ifreq == 1 && j != 4) )
    {
      cdoWarning("In checking last chunk:\n          Frequency of chunk file does not agree with frequency of the working file.\n          Switched to replace mode for this variable.");
      return 0;
    }


/* Encode in julseconds depending on frequency */
  if ( cdoVerbose) cdoPrint("Start to encode dates with frequencies to julseconds.");

  int old_start_year, old_start_month = 1, old_start_day = 1;
  int old_end_year, old_end_month = 1, old_end_day = 1;
  int new_end_year, new_end_month = 1, new_end_day = 1;

  switch ( j )
    {
    case ( 8 ):
      sscanf(old_start_date, "%04d%02d%02d", &old_start_year, &old_start_month, &old_start_day);
      sscanf(old_end_date, "%04d%02d%02d", &old_end_year, &old_end_month, &old_end_day);
      break;
    case ( 6 ):
      sscanf(old_start_date, "%04d%02d", &old_start_year, &old_start_month);
      sscanf(old_end_date, "%04d%02d", &old_end_year, &old_end_month);
      break;
    case ( 4 ):
      old_start_year = atol(old_start_date);
      old_end_year = atol(old_end_date);
      break;
    default:
      {
        cdoWarning("In checking last chunk:\n          Last chunk has subdaily frequency which is yet not enabled by cdo cmor.\n          Switched to replace mode for this variable.");
        return 0;
      }
    }

  int cdi_startdate = cdiEncodeDate(old_start_year, old_start_month, old_start_day);
  int cdi_enddate = cdiEncodeDate(old_end_year, old_end_month, old_end_day);
  int cdi_time = cdiEncodeTime(0, 0, 0);
  juldate_t julostart = juldate_encode(calendar, cdi_startdate, cdi_time);
  juldate_t juloend = juldate_encode(calendar, cdi_enddate, cdi_time);

  if ( cdoVerbose) cdoPrint("Successfully calculated juldates.", old_start_date);
/* Read in first vdate in case not piped */
  if ( cdoVerbose) cdoPrint("Start to calculate temporal gap between chunk and working file.");
  if ( cdoStreamName(0)->args[0] == '-' )
    {
      cdoWarning("Cdo cmor cannot enable append mode since you piped several cdo operators.\n          Switched to replace mode for this variable.");
      return 0;
    }
      
  int streamID2 = pstreamOpenRead(cdoStreamName(0));
  int vlistID2 = pstreamInqVlist(streamID2);
  int taxisID2 = vlistInqTaxis(vlistID2);
  juldate_t firstdate = juldate_encode(calendar, taxisInqVdate(taxisID2),
                                     taxisInqVtime(taxisID2));

/* Check temporal distance between last chunk date and first file date */
  double append_distance = juldate_to_seconds(juldate_sub(firstdate, juloend)) / 3600.0;

  if ( ( j == 8 && ( append_distance > 48.0 || append_distance < 0 ) )
     ||( j == 6 && ( append_distance/24.0 > 62.0 || append_distance < 0 ) )
     ||( j == 4 && ( append_distance/24.0/30.5 > 24.0 || append_distance < 0 ) ) )
    {
      cdoWarning("In checking the last chunk:\n          A temporal gap is diagnosed between end date of chunk file and first date of working file of: '%f' hours. Maximal valid gaps are:\n          48 hours for daily frequency\n          62 days for monthly frequency\n          24 month for yearly frequency.\n          Switched to replace mode for this variable.", append_distance);
      pstreamClose(streamID2);
      return 0;
    }

  if ( cdoVerbose) cdoPrint("Successfully checked temporal gap.");
/* Check file size */
  if ( cdoVerbose) cdoPrint("Start to check file size of chunk + working file.");
  double old_interval_sec = juldate_to_seconds(juldate_sub(juloend, julostart));
  double size_per_sec = (double) filesize / old_interval_sec;

  int maxsizegb = atol(kv_get_a_val(kvl, "ms", "2"));
  int maxsizeb = maxsizegb * 1024 * 1024 * 1024;

  int ntsteps = vlistNtsteps(vlistID2);
  if ( ntsteps < 0 )
    {
      ntsteps = 0;
      while ( pstreamInqTimestep(streamID2, ntsteps++)) ;
      if ( ntsteps == 0 )
        {
          cdoWarning("In checking whether append mode is possible:\n          No time steps found in infile.\n          Switched to replace mode for this variable.");
          pstreamClose(streamID2);
          return 0;
        }
    }
  
  double estimated_size;
  switch ( j )
    {
    case ( 8 ):
      estimated_size = ntsteps * 60 * 60 * 24 * size_per_sec + (double) filesize ;
      break;
    case ( 6 ):
      estimated_size = ntsteps * 60 * 60 * 24 * 30.5 * size_per_sec + (double) filesize;
      break;
    case ( 4 ):
      estimated_size = ntsteps * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize;
      break;
    default:
      {
        cdoWarning("In checking whether append mode is valid:\n          Selected chunk to append data has subdaily frequency which is yet not enabled by cdo cmor.\n          Switched to replace mode for this variable.");
        pstreamClose(streamID2);
        return 0;
      }
    }

  if ( (unsigned int)estimated_size > (unsigned int) maxsizeb )
    {
      cdoWarning("In checking whether append mode is valid:\n          Estimated file size of appended file is : '%f'gb and exceeds maximal allowed file size: '%d'gb.\n          Switched to replace mode for this variable.", estimated_size/1024.0/1024.0/1024.0, maxsizegb);
      pstreamClose(streamID2);
      return 0;
    }
  pstreamClose(streamID2);
  if ( cdoVerbose) cdoPrint("Successfully checked file size of chunk + working file.");
  return 1;
}

static char *use_chunk_des_files(list_t *kvl, int vlistID, int var_id, char *chunk_des_file, int ifreq, int calendar)
{
  char *chunk_file = (char *) Malloc(4096 * sizeof(char));
  if ( file_exist(chunk_des_file, 0, "chunk_description") )
    {
      FILE *fp = fopen(chunk_des_file, "r");
      size_t filesize = fileSize(chunk_des_file);
      char *buffer = (char*) Malloc(filesize);
      size_t nitems = fread(buffer, 1, filesize, fp);
      char *eof = readLineFromBuffer(buffer, &filesize, chunk_file, 4096);
      if ( eof != NULL )
        cdoWarning("In checking the last chunk:\n          Chunk description file contains more than one line.\n          All lines after line 1 are ignored.");
      fclose(fp);
      Free(buffer);
      if ( file_exist(chunk_file, 0, "chunk_description") && check_append_and_size(kvl, vlistID, chunk_file, ifreq, calendar) )
        return chunk_file;
      else
        cdoWarning("In checking the last chunk:\n          Chunk '%s' configured via chunk description file could either not be opened or is not suitable to be appended.\n          Switched to replace mode for this variable.", chunk_file);
    }
  else
    cdoWarning("Chunk description file '%s' could not be opened.\nSwitched to replace mode.", chunk_des_file);
  strcpy(chunk_file, " \0");
  return chunk_file;
}

static char **empty_array(struct mapping vars[], char ***chunk_files)
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    (*chunk_files)[i] = NULL;
  return *chunk_files;
}

static char **get_chunk_des_files(list_t *kvl, struct mapping vars[], char *miptab_freqptr, int nreq, int vlistID, char *charname)
{
  char **chunk_des_files = (char **) Malloc((nreq+1) * sizeof(char *));
  chunk_des_files[nreq] = NULL;

  char trunk[CMOR_MAX_STRING];
  const char *description_atts[] = {"model_id", "experiment_id", "member", NULL};
  strcpy(trunk, miptab_freqptr);
  for ( int i = 0; description_atts[i]; i++ )
    {
      strcat(trunk, "_");
      strcat(trunk, kv_get_a_val(kvl, description_atts[i], ""));
    }

  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++)
    {
      char *name = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
      if ( charname )
        strcpy(name, charname);
      else
        vlistInqVarName(vlistID, vars[j].cdi_varID, name);
      chunk_des_files[j] = (char *) Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(chunk_des_files[j], "CHUNK_FILE_%s_%s.txt\0", name, trunk);
      Free(name);
    }
  return chunk_des_files;  
}

static char **get_chunk_files(list_t *kvl, struct mapping vars[], int vlistID, int ifreq, int time_axis, int calendar, char *miptab_freqptr)
{
  int i = 0;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  char **chunk_files = (char **) Malloc((i+1) * sizeof(char *));
  chunk_files[i] = NULL;
  
  char *dummy = kv_get_a_val(kvl, "om", NULL);
  if ( !dummy || strcmp(dummy, "a") != 0 )
    return empty_array(vars, &chunk_files);
  else if ( time_axis == 4 )
    {
      cdoWarning("In validating append mode:\n          CMOR APPEND mode not possible for time independent variables.\n          Switched to replace mode for this variable");
      return empty_array(vars, &chunk_files);
    }

  if ( cdoVerbose )
    cdoPrint("Start to retrieve chunk files to append.\n");

  int num_aaf = 0;
  char **chunk_att_files = kv_get_vals(kvl, "lc", &num_aaf);
  char **chunk_des_files = NULL;
  if ( num_aaf != i && num_aaf > 0 )
    {
      cdoPrint("Number of chunk files '%d' disagree with number of requested variables '%d'.\n Switched to replace mode.\n", num_aaf, i); 
      return empty_array(vars, &chunk_files);
    }  
  else if ( num_aaf == 0 )
    {
      char *nd = kv_get_a_val(kvl, "d", "y");
/* For chunk description file : */
      if ( nd[0] == 'y' )
        chunk_des_files = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, NULL);
      else
        {
          cdoWarning("In getting chunk files:\n          Automatic chunk configuration via file not possible if DRS is not created.\n          Swichted to replace mode.");
          return empty_array(vars, &chunk_files);
        }
    }

  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++ )
    {
      if ( num_aaf != 0 )
        {
          if ( file_exist(chunk_att_files[j], 0, "chunk file") && check_append_and_size(kvl, vlistID, chunk_att_files[j], ifreq, calendar) )
            chunk_files[j] = strdup(chunk_att_files[j]);
          else
            {
              cdoWarning("Chunk '%s' could not be used.\n          Switched to replace mode for this variable.", chunk_att_files[j]);
              chunk_files[j] = strdup(" ");
            }   
        }
      else 
        {
          if ( cdoVerbose )
            cdoPrint("It is tried to open a chunk description file for varID: '%d': '%s'.", vars[j].cdi_varID, chunk_des_files[j]);
          chunk_files[j] = use_chunk_des_files(kvl, vlistID, vars[j].cdi_varID, chunk_des_files[j], ifreq, calendar);  
        }
      if ( cdoVerbose && strcmp(chunk_files[j], " ") != 0 )
        cdoPrint("Chunk file to append on var with CDI ID %d is: '%s'.", vars[j].cdi_varID, chunk_files[j]);
    }
  if ( chunk_des_files ) free_array(chunk_des_files);
  if ( cdoVerbose )
    cdoPrint("Successfully processed chunk file retrieval.");
  return chunk_files;
}

static void sigfunc(int sig)
{
  if ( sig == SIGTERM )
    cdoAbort("Program terminated by CMOR. A temporary ofile can outlive which needs to be deleted manually.");
}

static void write_variables(list_t *kvl, int *streamID, struct mapping vars[], int miptab_freq, int time_axis, int calendar, char *miptab_freqptr)
{
  int cmf = 0;
  int vlistID = pstreamInqVlist(*streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tsID = 0;
  int nrecs;
  size_t gridsize = vlistGridsizeMax(vlistID);

  if ( cdoVerbose )
    cdoPrint("10. Start to write variables via cmor_write.");
  if ( cdoVerbose )
    cdoPrint("10.1. Start to get frequency.");
  int sdate, stime, time_unit;
  get_taxis(kv_get_a_val(kvl, "rtu", NULL), &sdate, &stime, &time_unit);
  int tunitsec = get_tunitsec(time_unit);
  juldate_t ref_date = juldate_encode(calendar, sdate, stime);
  char *frequency = NULL;
  if ( time_axis != 4 )
    frequency = get_frequency(kvl, *streamID, vlistID, taxisID, miptab_freq);
  if ( cdoVerbose )
    cdoPrint("10.1. Successfully retrieved frequency.");


  int ifreq = 0;
  if ( frequency )
    {
      if ( strcmp(frequency,"yr") == 0 )
        ifreq = 1;
      if ( strcmp(frequency,"mon") == 0 )
        ifreq = 2;
      if ( strcmp(frequency,"day") == 0 )
        ifreq = 3;
    }

  if ( cdoVerbose )
    cdoPrint("10.2. Start to get chunk files.");
  char **chunk_files = get_chunk_files(kvl, vars, vlistID, ifreq, time_axis, calendar, miptab_freqptr);
  if ( cdoVerbose )
    cdoPrint("10.2. Successfully retrieved chunk files.");
  int i = 0;

  int zaxisID, zsize, pscheck = 1;
  char *charname = NULL;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    if ( vars[i].charvars )
      {
        if ( cdoVerbose )
          cdoPrint("10.3. Start to get auxiliary variables.");
        zaxisID = vlistInqVarZaxis(vlistID, vars[i].cdi_varID);
        zsize = zaxisInqSize(zaxisID);
        charname = (char *) Malloc(CDI_MAX_NAME * sizeof(char));
        vlistInqVarName(vlistID, vars[i].cdi_varID, charname);
        
        pstreamClose(*streamID);
        *streamID = pstreamOpenRead(cdoStreamName(0));
        pscheck = 0;
        if ( cdoVerbose )
          cdoPrint("10.3. Successfully retrieved auxiliary variables.");
        break;
      }
  if ( pscheck == 0 )
    cdoPrint("Since you defined a variable with character coordinate axis you cannot write another variable with zaxis of type ZAXIS_HYBRID.");

  if ( cdoVerbose )
    cdoPrint("10.4. Start to loop over time steps.");
  while ( (nrecs = pstreamInqTimestep(*streamID, tsID++)) )
    { 
      double time_bnds[2];
      double *time_bndsp;
      juldate_t jtime_val;
      double time_val;
      if ( time_axis != 4 )
        {
          jtime_val = get_cmor_time_val(taxisID, ref_date, tunitsec, calendar, frequency, tsID);
          time_val = juldate_to_seconds(juldate_sub(jtime_val, ref_date)) / tunitsec;
          time_bndsp = ( time_axis != 1 ) ? get_time_bounds(kvl, taxisID, frequency, ref_date, jtime_val, calendar, tunitsec, time_bnds, time_axis) : 0;
        }
      while ( nrecs-- )
        read_record(*streamID, vars, vlistID);

      int ps_index = -1;
      if ( pscheck )
        check_for_sfc_pressure(&ps_index, vars, vlistID, tsID);
      for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
        {
/*          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, vars[i].cdi_varID, name); */
          if ( !vars[i].help_var )
            {
              if ( time_axis != 4 )
                {
                  if ( vars[i].charvars )
                    {
                      void *dataslice = (void *) Malloc(gridsize * zsize * sizeof(double));
                      for ( int j = 0; j < gridsize * zsize; j++ )
                        ((double *)dataslice)[j] = ((double *)vars[i].data)[(tsID-1)*gridsize*zsize+j];
                      #if ( CMOR_VERSION_MAJOR == 2 )
                        cmf = cmor_write(vars[i].cmor_varID,
                       dataslice,
                       vars[i].datatype,
                       chunk_files[i],
                       1,
                       &time_val,
                       time_bndsp,
                       NULL);
                      Free(dataslice);
                      #elif ( CMOR_VERSION_MAJOR == 3 )
                        cmf = cmor_write(vars[i].cmor_varID,
                       dataslice,
                       vars[i].datatype,
                       1,
                       &time_val,
                       time_bndsp,
                       NULL);
                      Free(dataslice);
                      #endif
                    } 
                  else
                    {
                      #if ( CMOR_VERSION_MAJOR == 2 )
                        cmf = cmor_write(vars[i].cmor_varID,
                     vars[i].data,
                     vars[i].datatype,
                     chunk_files[i],
                     1,
                     &time_val,
                     time_bndsp,
                     NULL); 
                      #elif ( CMOR_VERSION_MAJOR == 3 )
                        cmf = cmor_write(vars[i].cmor_varID,
                     vars[i].data,
                     vars[i].datatype,
                     1,
                     &time_val,
                     time_bndsp,
                     NULL); 
                      #endif
                    }
                  if ( vars[i].zfactor_id > 0 )
                    {
                      #if ( CMOR_VERSION_MAJOR == 2 )
                        cmf = cmor_write(vars[i].zfactor_id,
                       vars[ps_index].data,
                       vars[ps_index].datatype,
                       chunk_files[i],
                       1,
                       &time_val,
                       time_bndsp,
                       &vars[i].cmor_varID);
                     #elif ( CMOR_VERSION_MAJOR == 3 )
                        cmf = cmor_write(vars[i].zfactor_id,
                       vars[ps_index].data,
                       vars[ps_index].datatype,
                       1,
                       &time_val,
                       time_bndsp,
                       &vars[i].cmor_varID);
                     #endif
                    }
                }
              else
                {
                  #if ( CMOR_VERSION_MAJOR == 2 )
                    cmf = cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   chunk_files[i], 0, 0, 0, NULL);
                  #elif ( CMOR_VERSION_MAJOR == 3 )
                    cmf = cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   0, 0, 0, NULL);
                  #endif 
                }
            }
        }
    }
  if ( cmf != 0 )
    cdoAbort("Function cmor_write failed!");
  if ( cdoVerbose )
    cdoPrint("10.4. Successfully looped over time steps.");
  if ( cdoVerbose )
    cdoPrint("10. Successfully written variables via cmor_write.");
  if ( cdoVerbose )
    cdoPrint("11. Start to close files, free allocated memory and, if necessary, write chunk files.");
  char **chunkdf = NULL;
  if ( strcmp(kv_get_a_val(kvl, "om", ""), "a") == 0 && strcmp(kv_get_a_val(kvl, "d", "y"), "y") == 0 )
    chunkdf = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, charname);

  char file_name[CMOR_MAX_STRING];
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    {
      if ( !vars[i].help_var )
        {
          cmf = cmor_close_variable(vars[i].cmor_varID, file_name, NULL);
          cdoPrint("     File stored in:  '%s' with cmor!", file_name);
          if ( chunkdf )
            {
              if ( cdoVerbose )
                cdoPrint("11.2. Start to write a chunk description file.");
              FILE *fp = fopen(chunkdf[i], "w+"); 
              if ( fp )
                fprintf(fp, "%s", file_name);
              else
                {
                  cdoPrint("Could not open a chunk description file '%s'.", chunkdf[i]);
                  continue;
                }
              fclose(fp);  
              if ( cdoVerbose )
                cdoPrint("11.2. Successfully written a chunk description file '%s'." , chunkdf[i]);            
            }         
        }
    }
  if ( cmf != 0 )
    cdoAbort("Function cmor_close_variable failed!");


  if (frequency) Free(frequency); if (chunk_files) free_array(chunk_files); if (chunkdf) free_array(chunkdf); if (charname) Free(charname);
  if ( cdoVerbose )
    cdoPrint("11. Successfully closed files and freed allocated memory.");
}

static list_t *check_for_charvars(list_t *maptab, const char *key)
{
/***/
/* If a mapping table variable selector (name or code) has more than one value, it must be a character coordinate*/
/* If it is given as a string and the string contains a ',', */
/* it must be divided into several values and is a variable with character coordinate */
/***/
  listNode_t *node = maptab->head;
  while ( node )
    {
      if ( node->data )
        {
          list_t *kvlist = *(list_t **)node->data;
          keyValues_t *kvn = NULL;
          if ( key )
            kvn = kvlist_search(kvlist, key);
          else
            {
              kvn = kvlist_search(kvlist, "name");
              if ( !kvn )
                kvn = kvlist_search(kvlist, "code");
            }
          if ( kvn && kvn->nvalues > 1 )
            return kvlist;

          if ( kvn && strstr(kvn->values[0], ",") && kvn->nvalues == 1 )
            {
              char *workchar2 = strdup(kvn->values[0]);
              Free(kvn->values[0]); Free(kvn->values);
              char *workchar = workchar2;
              char *thepoint = workchar2;
              int i = 0, j = 0;
              while ( *thepoint != '\0' )
                {
                  thepoint++;
                  if ( *thepoint == ',' )
                    j++;
                }
              j++;
              kvn->nvalues = j;
              kvn->values = (char **) malloc(kvn->nvalues*sizeof(char*)); 

              j = 0; thepoint = workchar;             
              while ( *thepoint != '\0' )
                {
                  if ( *thepoint == ',')
                    {
                      kvn->values[j] = (char *)Malloc( (i+1) * sizeof(char) );
                      strncpy(kvn->values[j], workchar, i);
                      kvn->values[j][i] = '\0';
                      j++; thepoint++; workchar+=i+1; i = 0;
                    }
                  else
                    {
                      thepoint++; i++;
                    }
                }
              if ( i > 0 )
                {
                  kvn->values[j] = (char *)Malloc( (i+1) * sizeof(char) );
                  strncpy(kvn->values[j], workchar, i);
                  kvn->values[j][i] = '\0';
                  workchar+=i; i = 0; j++; 
                }
              else
                {
                  Free(workchar2);
                  cdoWarning("In checking for variables with character coordinate:\n          Names in String for key '%s' could not be interpreted correctly due to a comma at end of line.", key);
                  return NULL;
                }
              Free(workchar2);
              return kvlist;
            }
        }
      node = node->next;
    } 
  return NULL;
}

static void read_maptab(list_t *kvl, int streamID, char *miptabfreq, struct mapping vars[])
{
/***/
/* Build mapping table from a combination of two attributes if mt does not begin with / and a directory path is given */
/***/
  if ( cdoVerbose )
    cdoPrint("5. Start to find, read and apply mapping table.");
  char *maptab = kv_get_a_val(kvl, "mt", NULL);
  char *maptabdir = kv_get_a_val(kvl, "mapping_table_dir", NULL);
  char *maptabbuild = NULL;
  keyValues_t *kvn = kvlist_search(kvl, "n");
  keyValues_t *kvc = kvlist_search(kvl, "c");
  keyValues_t *kvcn = kvlist_search(kvl, "cn");
  int byteorder;
  int filetype = cdiGetFiletype(cdoStreamName(0)->args, &byteorder);

  if ( maptab && maptabdir ) if ( maptab[0] != '/' )
    {
      maptabbuild = (char *) Malloc((strlen(maptab)+strlen(maptabdir)+2) * sizeof(char));
      sprintf(maptabbuild, "%s/%s\0", maptabdir, maptab);
    }
  if ( maptab )
    {
      if ( maptabbuild ) maptab = maptabbuild;
      int vlistID = pstreamInqVlist(streamID);

/***/
/* Parse the table as a fortran namelist wich contains lists (=lines) of keyvalues */
/***/
      if ( cdoVerbose )
        cdoPrint("5.1 Try to read mapping table: '%s'", maptab);
      kv_insert_a_val(kvl, "workfile4err", maptab, 1);
      list_t *pml = cdo_parse_cmor_file(maptab, kvl);
      if ( pml == NULL )
        {
          cdoWarning("5.1. In parsing the mapping table '%s':\n          Mapping table could not be parsed. Operator continues.", maptab);
          return;
        }
      const char *ventry[] = {"&parameter"};
      int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);
      list_t *charvarlist = NULL; 
/***/
/* If a variable selector name or code is given in cmdline, the corresponding variable is picked from Infile and mapped. */
/* Only the first value of name/code given in the cmdline is processed */
/* If no variable selector is given, process all variables and map via name and code */
/***/
/* However, if the mapping table contains a keyvalue pair for name or code with more than one value, */
/* the corresponding variable has a character coordinate and requires special treatment */
/* This is tested once before mapping. If the special variable equals the variable which is to map,
/* the special treatment begins with fct addcharvar */
/***/
/* Different CMOR variables are built with one model variable. */
/* Consequently, for one model variable more than one mapping table entry can exist */
/* As a second identification argument, the mapping table name (miptabfreq) is used */
/***/
/* If no variable selector is given in the mapping table, it is assumed that the infile variable is already named like cmor_name */
/***/
      if ( kvn )
        {
          if ( filetype == FILETYPE_GRB || filetype ==  FILETYPE_GRB2 )
            cdoPrint("5.1. In applying the mapping table:\n          Note that you use 'name' as selector keyword allthough the type of infile is GRB.");
          if ( charvarlist = check_for_charvars(pml, "name") )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          if ( kvn->nvalues > 1 )
            cdoWarning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline variable selection key 'name' is processed.", maptab);
          maptab_via_cmd(pml, kvn->values[0], vlistID, vlistNvars(vlistID),  "name", kvcn->values[0], miptabfreq);
          if ( cdoVerbose )
            cdoPrint("5. Successfully found, read and applied mapping table '%s'.", maptab);
        }
      else if ( kvc )
        {
          if ( charvarlist = check_for_charvars(pml, "code") )
            {
              keyValues_t *charkvc = kvlist_search(charvarlist, "code");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvc || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvc, vlistID, "code", vars);
            }
          if ( kvc->nvalues > 1 )
            cdoWarning("5.1. In applying the mapping table '%s':\n          Only the first value of commandline variable selection key 'code' is processed.", maptab);
          maptab_via_cmd(pml, kvc->values[0], vlistID, vlistNvars(vlistID), "code", kvcn->values[0], miptabfreq);
          if ( cdoVerbose )
            cdoPrint("5. Successfully found, read and applied mapping table '%s'.", maptab);
        }
      else if ( kvcn )
        { 
          if ( charvarlist = check_for_charvars(pml, NULL) )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( strcmp(charkvcn->values[0], kvcn->values[0]) == 0 )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          maptab_via_cn(pml, kvcn->values, vlistID, vlistNvars(vlistID), kvcn->nvalues, miptabfreq, filetype); 
          if ( cdoVerbose )
            cdoPrint("5. Successfully found, read and applied mapping table '%s'.", maptab);
        }
      else
        {
          if ( charvarlist = check_for_charvars(pml, NULL) )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
            {
/***/
/* Begin with Code in case infile is of type GRB */
/***/
              if ( filetype == FILETYPE_GRB || filetype ==  FILETYPE_GRB2 )
                if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "code", miptabfreq) )
                  {
                    if ( cdoVerbose )
                      cdoPrint("5.1. Successfully mapped varID '%d' via code.", varID);
                    continue;
                  }
              if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "name", miptabfreq) )
                {
                  if ( cdoVerbose )
                    cdoPrint("5.1. Successfully mapped varID '%d' via name.", varID);
                  continue;
                }
              if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "code", miptabfreq) )
                {
                  if ( cdoVerbose )
                    cdoPrint("5.1. Successfully mapped varID '%d' via code.", varID);
                  continue;
                }
/***/
/* In case corresponding mapping table entry does not contain a variable selector attribute */
/***/
              if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "cmor_name", miptabfreq) )
                {
                  if ( cdoVerbose )
                    cdoPrint("5.1. Successfully mapped varID '%d' via cmor_name.", varID);
                  continue;
                }
              cdoWarning("5.1. In applying the mapping table '%s':\n          Could not map variable with id '%d'.", maptab, varID);
            }
        }
/***/
/* In case a requested variable needs an auxilliary variable, the latter may be mapped later. */
/* If a mapping table exists is saved here */
/***/
      kv_insert_a_val(kvl, "mtproof", maptab, 1);
      list_destroy(pml);
      if ( maptabbuild ) Free(maptabbuild);
    }
  else if ( cdoVerbose )
    cdoPrint("5. No mapping table found.");
}

static void parse_cmdline(list_t *pml, char **params, int nparams, const char *ventry)
{
  list_t *kvl = NULL;
  kvl = list_new(sizeof(keyValues_t *), free_keyval, ventry);
  list_append(pml, &kvl);

  char *key = NULL, *eqpos = NULL;
  char **values = NULL;
  int i = 1, j = 0;
  for ( i = 1; i < nparams; i++ )
    {
      if ( eqpos = strchr(params[i], '=')  )
        {
          if ( key && values[0] )
            {
              const char *short_key = check_short_key(key);
              if ( short_key )
                {
                  if ( strcmp(short_key, key) != 0 )
                    {
                      Free(key);
                      key = strdup(short_key);
                    }
                  kvlist_append(kvl, (const char *)key, (const char **) values, j);
                }
              else
                cdoWarning("Unknown commandline keyword: '%s'\n", key);
              Free(key);
              free_array(values);
            }
          else if ( key )
            cdoAbort("Could not find values for commandline keyword: '%s'.", key);
          if ( strlen(eqpos) == 1 )
            cdoAbort("Could not find values for commandline parameter: '%s'.", params[i]);
          key = strdup(strtok(params[i], "="));
          values = (char **) Malloc(100 * sizeof(char *));
          j = 0;   
          int errh = copy_value(strtok(NULL, ""), values, &j);
          if ( errh > 0 )
            handleError(NULL, errh, NULL);
        }
      else
        {
          if ( !key )
            cdoAbort("Found no key for value '%s'.", params[i]);
          else
            {
              int errh = copy_value(params[i], values, &j);
              if ( errh > 0 )
                handleError(NULL, errh, NULL);
            }
        }
    }
  if ( key && values )
    {
      const char *short_key = check_short_key(key);
      if ( short_key )
        {
          if ( strcmp(short_key, key) != 0 )
            {
              Free(key);
              key = strdup(short_key);
            }
          kvlist_append(kvl, (const char *)key, (const char **) values, j);
        }
      Free(key);
      free_array(values);
    }
  else if ( values )
    cdoAbort("Found no commandline keyword for value '%s'.", params[i-1]);
}

static char *get_mip_table(char *params, list_t *kvl, char *project_id)
{
  char *miptab;
  if ( cdoVerbose )
    cdoPrint("2.2. Start to find a MIP table file.");
  if ( !params )
    cdoAbort("In finding the MIP table:\n          A mip table name or path is required as first argument. No first argument found.");
  if ( file_exist(params, 0, "MIP table") )
    {
      miptab = strdup(params);
      if ( cdoVerbose )
        cdoPrint("2.2. MIP table file '%s' exists.", miptab);
      return miptab;
    }
  else
    {
      cdoPrint("In finding the MIP table:\n          Your first argument is not an existing mip table file.\n          It is tried to build a path with additional configuration attributes 'mip_table_dir' and 'project_id'");
      char *miptabdir = kv_get_a_val(kvl, "mip_table_dir", NULL);
      if ( miptabdir && project_id )
        {
#if ( CMOR_VERSION_MAJOR == 2 )
          {
          miptab = (char *)Malloc((strlen(miptabdir)+strlen(project_id)+strlen(params)+3) * sizeof(char));
          sprintf(miptab, "%s/%s_%s\0", miptabdir, project_id, params);
          }
#elif ( CMOR_VERSION_MAJOR == 3 )
          {
          miptab = (char *)Malloc((strlen(miptabdir)+strlen(project_id)+strlen(params)+8) * sizeof(char));
          sprintf(miptab, "%s/%s_%s.json\0", miptabdir, project_id, params);
          }
#endif
          file_exist(miptab, 1, "MIP table");
          if ( cdoVerbose )
            cdoPrint("2.2. MIP table file '%s' exists.", miptab);
          return miptab;
        }
      else
        cdoAbort("In finding the MIP table:\n          Could not find attribute 'mip_table_dir'.");
    }   
}

static char *freq_from_path(char *mip_table)
{
  char *freq = mip_table;
  int fpos = 0, k = 0, j = 0;
  while ( *(mip_table + j) )
    {
      j++;
      if ( *(mip_table + j) == '/' )
        k = j + 1;
      if ( *(mip_table + j) == '_' && *(mip_table + j + 1) )
        fpos = j + 1;
    }
  freq += k;
  if ( fpos > k )
    freq += fpos-k;
  return freq;
}

static int get_miptab_freq(list_t *kvl, char *mip_table, char *project_id)
{
  int miptab_freq = 0;
  char *freq = freq_from_path(mip_table);
  if ( freq != NULL )
    {
      if ( strstr(freq, "yr") || strstr(freq, "Yr") )
        miptab_freq = 11;
      else if ( strstr(freq, "mon") || strstr(freq, "Mon") )
        miptab_freq = 12;
      else if ( strstr(freq, "day") || strstr(freq, "Day") )
        miptab_freq = 13;
      else if ( strstr(freq, "6hr") )
        miptab_freq = 14;
      else if ( strstr(freq, "3hr") )
        miptab_freq = 15;

      if ( strcmp(freq, "Oclim") == 0 )
        miptab_freq = 1;
      else if ( strcmp(freq, "Oyr") == 0 )
        miptab_freq = 2;
      else if ( strcmp(freq, "cfMon") == 0 )
        miptab_freq = 3;
      else if ( strcmp(freq, "day") == 0 )
        miptab_freq = 4;
      else if ( strcmp(freq, "6hrPlev") == 0 && strcmp(project_id, "CMIP5") == 0 )
        miptab_freq = 5;
      else if ( strcmp(freq, "6hrPlevPt") == 0 )
        miptab_freq = 5;
      else if ( strcmp(freq, "6hrLev") == 0 )
        miptab_freq = 6;
      else if ( strcmp(freq, "E1hrClimMon") == 0 )
        miptab_freq = 7;
    }
  return miptab_freq;
}

static void check_cmdline_mapping(list_t *kvl)
{
  char *name = kv_get_a_val(kvl, "n", NULL);
  char *code = kv_get_a_val(kvl, "c", NULL);
  char *cn = kv_get_a_val(kvl, "cn", NULL);
  if ( ( name && code ) )
    cdoAbort("Mapping via command line failed. Only one variable selector of 'name' and 'code' is allowed.");
  if ( ( name && !cn ) || ( code && !cn ) )
    cdoAbort("Mapping via command line failed. A corresponding 'cmor_name' is needed.");
}

static char *get_project_id(list_t *kvl)
{
  if ( cdoVerbose )
    cdoPrint("2.1. Start to check whether 'project_id' or 'mip_era' is denoted.");
  char *project_id, *dummy, *dummy2;
  dummy = kv_get_a_val(kvl, "project_id", NULL);
  dummy2 = kv_get_a_val(kvl, "mip_era", NULL);
#if defined(CMOR_VERSION_MAJOR)
#if ( CMOR_VERSION_MAJOR == 2 )
{
  if ( !dummy && !dummy2)
    cdoAbort("Attribute 'project_id' is required.");
  else if ( !dummy )
    cdoAbort("Cannot produce CMIP6 standard with CMOR2.\n          Value for attribute 'project_id' is required.");
  else
    project_id = strdup(dummy);
}
#elif ( CMOR_VERSION_MAJOR == 3 )
{
  if ( !dummy && !dummy2)
    cdoAbort("Attribute 'mip_era' or 'project_id' is required.");
  else if ( !dummy2 )
    {
      cdoWarning("You try to produce CMIP5 standard with CMOR3.\n          It is recommended to use CMOR2 for this job instead.");
      project_id = strdup(dummy);
    }
  else
    project_id = strdup(dummy2);
}
#endif
#else
  cdoAbort("Cannot check CMOR version: Missing makro CMOR_VERSION_MAJOR");
#endif

  if ( cdoVerbose )
    cdoPrint("2.1. Successfully found project_id / mip_era: '%s'.", project_id);
  return project_id;
}


static int cmor_load_and_set_table(list_t *kvl, char *param0, char *project_id, char **mip_table)
{
  int table_id = 0, cmf = 0;
#if ( CMOR_VERSION_MAJOR == 3 )
  Free(*mip_table);
  cdoPrint("Need to get the mip_table once more.");
  *mip_table = get_mip_table(param0, kvl, project_id);
#endif
  cmf = cmor_load_table(*mip_table, &table_id);
  if ( cmf != 0 )
    cdoAbort("Function cmor_load_table failed!");
  cmf = cmor_set_table(table_id);
  if ( cmf != 0 )
    cdoAbort("Function cmor_set_table failed!");
  return table_id;
}

static void removeDataset()
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  int procID = getpid();
  char *dataset_path = (char *) Malloc( (strlen(cwd) + 1 + strlen("dataset.json") + 3) * sizeof(char));;
  sprintf(dataset_path, "%s/dataset%d.json\0", cwd, procID);
  remove(dataset_path);
}
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  signal(SIGTERM, sigfunc);
  int nparams = operatorArgc();
  char **params = operatorArgv();
  char *miptableInput = strdup(params[0]);

  /* Definition of pml: */
  list_t *pml = list_new(sizeof(list_t *), free_kvlist, "pml");

  if ( nparams < 1 ) cdoAbort("Too few arguments!");

  /* Define kvl and read cmdline */
  parse_cmdline(pml, params, nparams, "cmdline");
  
  const char *pmlistHelper[] = {"cmdline"};
  /* Get kvl and use it from now on instead of pml */
  list_t *kvl = pmlist_get_kvlist_ventry(pml, 1, pmlistHelper);

  /* Check whether a command line mapping is active */
  check_cmdline_mapping(kvl);

  /* Config files are read with descending priority. */
  read_config_files(kvl);

  /* Get project_id, mip_table and mip_table frequency*/
  if ( cdoVerbose )
    cdoPrint("2. Start to find a MIP table and to deduce a frequency from MIP table file.");
  char *project_id = get_project_id(kvl);
  char *mip_table  = get_mip_table(miptableInput, kvl, project_id);
#if ( CMOR_VERSION_MAJOR == 3 )
  mip_table[strlen(mip_table)-5] = '\0';
#endif
  int miptab_freq  = get_miptab_freq(kvl, mip_table, project_id);

  char *miptab_freqptr = strdup(freq_from_path(mip_table));
  kv_insert_a_val(kvl, "miptab_freq", miptab_freqptr, 1);

  if ( cdoVerbose )
    cdoPrint("2. Successfully found a MIP table '%s' and deduced a MIP table frequency '%s'.", mip_table, miptab_freqptr);

  if ( cdoVerbose )
    cdoPrint("3. Start to open infile '%s'.", cdoStreamName(0)->args);
  int streamID = pstreamOpenRead(cdoStreamName(0));
  if ( cdoVerbose )
    cdoPrint("3. Successfully opened infile '%s'.", cdoStreamName(0)->args);
  /* Short keys from rtu, mt, gi must be included similar to global atts */
  add_globalhybrids(kvl);

  /* Check for attributes and member name */
  if ( cdoVerbose )
    cdoPrint("4. Start to check attributes.");
  check_attr(kvl, project_id);
  check_mem(kvl, project_id);
  if ( cdoVerbose )
    cdoPrint("4. Successfully checked global attributes.");

 /* dump_global_attributes(pml, streamID); */

  struct mapping *vars = construct_var_mapping(streamID);

 /* read mapping table */
  read_maptab(kvl, streamID, miptab_freqptr, vars);

  int time_axis = 0, calendar = 0;

  setup_dataset(kvl, streamID, &calendar);

  int table_id = cmor_load_and_set_table(kvl, miptableInput, project_id, &mip_table);

  register_all_dimensions(kvl, streamID, vars, table_id, project_id, miptab_freq, &time_axis);
  write_variables(kvl, &streamID, vars, miptab_freq, time_axis, calendar, miptab_freqptr);

#if ( CMOR_VERSION_MAJOR == 3 )
  removeDataset();  
#endif
  destruct_var_mapping(vars);
  Free(mip_table);
  Free(project_id); 
  Free(miptab_freqptr);
 /* Free(miptableInput); */
  list_destroy(pml); 

  pstreamClose(streamID);
#else
  cdoWarning("CMOR support not compiled in!");
#endif
  cdoFinish();
  return 0;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
