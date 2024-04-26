#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "cdo_int.h"
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#if defined(HAVE_LIBXML2)
#include <libxml/parser.h>
#include <libxml/tree.h>
xmlNode *root_node;
xmlDoc *param_doc;
#endif


#define DBG_MSG 0 

void *magics_node, *results_node;


int init_XMLtemplate_parser( char *Filename )
{
#if defined(HAVE_LIBXML2)
  param_doc = xmlReadFile( Filename, NULL, 0 );
  if ( param_doc == NULL )
    {
      printf( "Error: Could not parse the file \"%s\"\n", Filename );
      return (1);
    }
  else
    {
      fprintf( stderr, "XML file %s being parsed \n", Filename );
      root_node = xmlDocGetRootElement( param_doc );
    }
#else
  
  cdoAbort("XML2 support not compiled in!");
  
#endif

  return 0;
}


int updatemagics_and_results_nodes(void)
{
#if defined(HAVE_LIBXML2)
  xmlNode *cur_node = NULL;
	
  if( root_node == NULL )
    {
      printf( "Invalid Root Node\n" );
      return 0;
    }

  for ( cur_node = root_node->children; cur_node; cur_node = cur_node->next )
    {   
      if ( cur_node->type == XML_ELEMENT_NODE )
        {   
#if DBG_MSG
          fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
          if( !strcmp( (const char*)cur_node->name, "magics" ) ) 
            {
              magics_node = (void*) cur_node;
#if DBG_MSG
              fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
	    }  

          if( !strcmp( (const char*)cur_node->name, "results" ) ) 
            {
              results_node = (void*) cur_node;
#if DBG_MSG
              fprintf( stdout, "Node Name: %s \n", cur_node->name );
#endif
	    }  
	}
    }
#else
  
  cdoAbort("XML2 support not compiled in!");
  
#endif

  return 0;
}


int quit_XMLtemplate_parser(void)
{
#if defined(HAVE_LIBXML2)
  xmlFreeDoc( param_doc );
  xmlCleanupParser( );
  if( param_doc == NULL )
    printf( "Cleaned XML parser\n" );
#if DBG_MSG
  fprintf( stdout, "Cleaned XML parser\n" );
#endif
#else
  
  cdoAbort("XML2 support not compiled in!");
  
#endif

  return 0;
}
