#ifndef __NAMELIST_H_
#define __NAMELIST_H_

enum namelisttype {
  NAMELIST_UNDEFINED = 0,
  NAMELIST_OBJECT    = 1,
  NAMELIST_KEY       = 2,
  NAMELIST_STRING    = 3,
  NAMELIST_WORD      = 4
};


enum namelisterr {
  NAMELIST_ERROR_INVAL = -1,   // Invalid character inside NAMELIST string/word
  NAMELIST_ERROR_PART  = -2,   // The string is not a full NAMELIST packet, more bytes expected 
  NAMELIST_ERROR_INKEY = -3,   // Invalid character inside NAMELIST key
  NAMELIST_ERROR_INTYP = -4,   // Invalid NAMELIST key type
  NAMELIST_ERROR_INOBJ = -5,   // Invalid NAMELIST object
  NAMELIST_ERROR_EMKEY = -6    // Empty key name
};

// NAMELIST token description.
typedef struct {
  int type;            // type (object, key, string word)
  int start;           // start position in NAMELIST buffer
  int end;             // end position in NAMELIST buffer
} namelisttok_t;


typedef struct {
  namelisttok_t *tokens;
  unsigned int num_tokens;
  unsigned int toknext;
  unsigned int pos;
  unsigned int lineno;
} namelist_parser;


namelist_parser *namelist_new(void);

void namelist_destroy(namelist_parser *parser);

int namelist_parse(namelist_parser *parser, const char *buf, size_t len);

void namelist_dump(namelist_parser *parser, const char *buf);

int namelist_verify(namelist_parser *parser, const char *buf);

#endif // __NAMELIST_H_
