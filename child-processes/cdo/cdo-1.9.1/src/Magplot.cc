#if defined(HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


#if defined(HAVE_LIBMAGICS)

#include "magics_api.h"

#include "magics_template_parser.h"
#include "results_template_parser.h"
#include "StringUtilities.h"

#define DBG 0

/***** ADDED for handling plots with  defined lat lon min max *****/
/*** LAT_MIN,LAT_MAX, LON_MIN,LON_MAX ****/
/**** lat_min,lat_max,lon_min,lon_max ****/

/****  
subpage_lower_left_latitude 
subpage_lower_left_longitude
subpage_upper_right_latitude
subpage_upper_right_longitude
****/

int CONTOUR, SHADED, GRFILL;

const char  *contour_params[] = {"min","max","count","interval","list","colour","thickness","style","RGB","device", "step_freq","file_split","lat_min","lat_max","lon_min","lon_max","projection"};
int contour_param_count = sizeof(contour_params)/sizeof(char*);

const char  *shaded_params[] = {"min","max","count","interval","list","colour_min","colour_max","colour_table","RGB","colour_triad","device","step_freq","file_split","lat_min","lat_max","lon_min","lon_max","projection"};
int shaded_param_count = sizeof(shaded_params)/sizeof(char*);

const char  *grfill_params[] = {"min","max","count","interval","list","colour_min","colour_max","colour_table","resolution","RGB","colour_triad","device","step_freq","file_split","lat_min","lat_max","lon_min","lon_max","projection"};
int grfill_param_count = sizeof(grfill_params)/sizeof(char*);

const char  *STD_COLOUR_TABLE[] = {"red", "green", "blue", "yellow", "cyan", "magenta", "black", "avocado",
			     "beige", "brick", "brown", "burgundy",
			     "charcoal", "chestnut", "coral", "cream", 
			     "evergreen", "gold", "grey", 
			     "khaki", "kellygreen", "lavender",
			     "mustard", "navy", "ochre", "olive",
			     "peach", "pink", "rose", "rust", "sky",
			     "tan", "tangerine","turquoise",
			     "violet", "reddishpurple",
			     "purplered", "purplishred",
			     "orangishred", "redorange", "reddishorange",
			     "orange", "yellowishorange",
			     "orangeyellow", "orangishyellow", 
			     "greenishyellow", "yellowgreen",
			     "yellowishgreen", "bluishgreen",
			     "bluegreen", "greenishblue",
			     "purplishblue", "bluepurple",
			     "bluishpurple", "purple", "white"
			    };


char **USR_COLOUR_TABLE = NULL;

int  STD_COLOUR_COUNT = sizeof( STD_COLOUR_TABLE )/sizeof( char* );
int  USR_COLOUR_COUNT =0;


const char *STYLE_TABLE[] = { "SOLID","DASH","DOT","CHAIN_DASH","CHAIN_DOT"};
int STYLE_COUNT = sizeof( STYLE_TABLE )/ sizeof( char *);

const char *DEVICE_TABLE[] = { "PS","EPS","PDF","PNG","GIF","GIF_ANIMATION","JPEG","SVG","KML"};
int DEVICE_COUNT = sizeof( DEVICE_TABLE )/ sizeof( char *);


/*char *PROJECTION_TABLE[] = { "cylindrical", "polar_stereographic", "polar_north", "geos", "meteosat", "meteosat_57E", "goes_east", "lambert", "EPSG3857", "goode", "collignon", "mollweide", "robinson", "bonne", "google", "efas", "EPSG4326", "lambert_north_atlantic", "mercator", "cartesian", "taylor", "tephigram" };
*/
/** The following projections are having some issues to be clarified with Magics++ **/

const char *PROJECTION_TABLE[] = { "cylindrical", "polar_stereographic", "polar_north", "geos", "meteosat", "meteosat_57E", "lambert", "EPSG3857", "goode", "collignon", "mollweide", "robinson", "bonne", "google", "efas", "EPSG4326", "lambert_north_atlantic", "mercator" };
int PROJECTION_COUNT = sizeof( PROJECTION_TABLE )/ sizeof( char *);



int ANIM_FLAG = 0, STEP_FREQ = 0; /* '0' for static images like jpeg,ps, etc.. , '1' for animation formats */


int checkcolour( char *colour_in );
int ReadColourTable ( char *filepath );
int checkstyle( char *style_in );
int checkdevice( char *device_in );
int checkprojection( char *projection_in );

 /* Magics default values */
int COUNT = 10, isRGB = FALSE,   THICKNESS = 1, NUM_LEVELS = 0, FILE_SPLIT = FALSE;
double YMIN = 1.0e+200, YMAX = -1.0e+200, INTERVAL = 8.0, RESOLUTION = 10.0f, *LEV_LIST = NULL ;
double LAT_MIN = 1.0e+200, LAT_MAX = -1.e+200;
double LON_MIN = 1.0e+200, LON_MAX = -1.e+200;
const char *COLOUR = NULL, *COLOUR_MIN = NULL, *COLOUR_MAX = NULL, *STYLE = NULL, *DEVICE = NULL, *COLOUR_TRIAD = NULL, *PROJECTION = NULL;


static
void magplot( const char *plotfile, int operatorID, const char *varname, const char *units, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *array,  int nparam, char **params, char *datetime, bool lregular)

{
  long i;
  double dlon = 0, dlat = 0;
  char plotfilename[4096];
  char *titlename;
  int j, split_str_count;
  const char *sep_char = "=";
  char **split_str = NULL;
  char tempname[256];
  
  
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", params[i]);
      fflush( stderr );
  
      for( i = 0; i < nparam; ++i )
        {
           split_str_count = 0;
           sep_char = "=";
           split_str_count = StringSplitWithSeperator( params[i], sep_char, &split_str );
	
           if ( !strcmp( split_str[0],"min" ) )        fprintf(stderr,"Min Val %g\n",YMIN );
           if ( !strcmp( split_str[0],"max" ) )        fprintf(stderr,"Max Val %g\n",YMAX );
           // if ( !strcmp( split_str[0],"resolution" ) ) fprintf( stderr,"RESOLUTION %g\n",RESOLUTION );
           if ( !strcmp( split_str[0],"colour" ) )     fprintf(stderr,"COLOUR %s\n",COLOUR );
           if ( !strcmp( split_str[0],"colour_min" ) ) fprintf(stderr,"COLOUR %s\n",COLOUR_MIN );
           if ( !strcmp( split_str[0],"colour_max" ) ) fprintf(stderr,"COLOUR %s\n",COLOUR_MAX );
           if ( !strcmp( split_str[0],"interval" ) )   fprintf(stderr,"INTERVAL %f\n",INTERVAL );
           if( !strcmp( split_str[0],"count" ) )       fprintf(stderr,"COUNT %d\n",COUNT );
	
           if ( !strcmp( split_str[0],"list" ) ) 
	     {
               for( j = 0; j < split_str_count; j++ )  fprintf(stderr,"LIST %f\n",LEV_LIST[j] ); 
	     }
	
           if ( !strcmp( split_str[0],"thickness" ) )  fprintf(stderr,"THICKNESS %d\n",THICKNESS );
           if ( !strcmp( split_str[0],"style" ) )      fprintf(stderr,"STYLE %s\n",STYLE );
           if ( !strcmp( split_str[0],"device" ) )     fprintf(stderr,"DEVICE %s\n",DEVICE );
           if ( !strcmp( split_str[0],"step_freq" ) )  fprintf(stderr,"STEP_FREQ %d\n",STEP_FREQ );
           if ( !strcmp( split_str[0],"lat_min" ) )    fprintf(stderr,"Lat Min Val %g\n",LAT_MIN );
           if ( !strcmp( split_str[0],"lat_max" ) )    fprintf(stderr,"Lat Max Val %g\n",LAT_MAX );
           if ( !strcmp( split_str[0],"lon_min" ) )    fprintf(stderr,"Lon Min Val %g\n",LON_MIN );
           if ( !strcmp( split_str[0],"lon_max" ) )    fprintf(stderr,"Lon Max Val %g\n",LON_MAX );
           if ( !strcmp( split_str[0],"projection" ) ) fprintf(stderr,"PROJECTION %s\n",PROJECTION );

           Free(split_str);
        }
    }

  if ( nlon > 1 )
    {
      for ( i = 1; i < nlon; ++i ) dlon += (grid_center_lon[i] - grid_center_lon[i-1]);
      dlon /= (nlon-1);
    }
  if ( nlat > 1 )
    {
      for ( i = 1; i < nlat; ++i ) dlat += (grid_center_lat[nlon*i] - grid_center_lat[nlon*(i-1)]);
      dlat /= (nlat-1);
    }

  sprintf( plotfilename, "%s [%s] %s", varname, units, datetime );
  titlename = strdup( plotfilename );
  sprintf( plotfilename, "%s_%s", plotfile, varname );

  mag_setc ("output_name",      plotfilename);
  mag_new( "page");

  /* Set the input data arrays to magics++ */
   
  mag_set2r("input_field", array, nlon, nlat);

  if ( lregular )
    {
      mag_setc("input_field_organization", "REGULAR");
      // mag_setc("input_field_organization", "GAUSSIAN");

      mag_setr("input_field_initial_latitude", grid_center_lat[0]);
      mag_setr("input_field_latitude_step", dlat);
        
      mag_setr("input_field_initial_longitude", grid_center_lon[0]);
      mag_setr("input_field_longitude_step", dlon);
     }
  else
    {
      mag_setc("input_field_organization", "NONREGULAR");

      mag_set2r("input_field_latitudes", grid_center_lat, nlon, nlat);
      mag_set2r("input_field_longitudes", grid_center_lon, nlon, nlat);
    }
  
  /* magics_template_parser( magics_node ); */
  /* results_template_parser(results_node, varname ); */

  /* set up the coastline attributes */
  /* mag_setc ("map_coastline_colour", "khaki"); */
  /* mag_setc ("map_grid_colour",      "grey");  */ 

  /* Parameters common to all operators */
  if( DEVICE )                                
    {
      mag_setc ("output_format", DEVICE );
    }

  if( PROJECTION )
    {
      mag_setc ( "subpage_map_projection", PROJECTION );
    }

  mag_seti ("map_label_latitude_frequency",2);
  mag_seti ("map_label_longitude_frequency",2);
  /*mag_setr ("map_label_height",0.5);*/
  mag_setr ("map_label_height",0.4);

  /* define the contouring parameters */
  if ( operatorID == SHADED )
    {
      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );
      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );

      if( LAT_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_latitude", LAT_MIN );
        }

      if( LON_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_longitude", LON_MIN );
        }

      if( LAT_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_latitude", LAT_MAX );
        }

      if( LON_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_longitude", LON_MAX );
        }
      
      if( YMIN < 1.0e+200  )
        {
	   mag_setr( "contour_shade_min_level", YMIN );
	   mag_setr( "contour_min_level", YMIN );
        }
      
      if( YMAX > -1.0e+200 )
        {
	   mag_setr( "contour_shade_max_level", YMAX );
	   mag_setr( "contour_max_level", YMAX );
        }
      
      if( COLOUR_MAX )
	mag_setc( "contour_shade_max_level_colour", COLOUR_MAX );

      if( COLOUR_MIN )
	mag_setc( "contour_shade_min_level_colour", COLOUR_MIN );

      if( IS_NOT_EQUAL(INTERVAL, 8.0f) )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( USR_COLOUR_COUNT ) 
	{
	  mag_setc( "contour_shade_colour_method", "LIST" );
	  mag_set1c( "contour_shade_colour_list",( const char **)USR_COLOUR_TABLE, USR_COLOUR_COUNT ); 
	}
          
      if( COLOUR_TRIAD )                                
	{
	  mag_setc( "contour_shade_colour_direction", COLOUR_TRIAD );
	}
	
      
      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

      if( DBG )
        {
          mag_enqc ( "output_name", (char*)&tempname );
          fprintf( stderr, " SHADED Done %s!\n",tempname );
          fprintf( stderr, " SHADED Done!\n" );
        }
    }
  else if ( operatorID == CONTOUR )
    {
      mag_setc ("contour",                  "on");
      mag_setc ("contour_shade",            "off");
      mag_setc ("contour_label",            "on");
      mag_setc ("contour_highlight",        "off");
      
      if( LAT_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_latitude", LAT_MIN );
        }

      if( LON_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_longitude", LON_MIN );
        }

      if( LAT_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_latitude", LAT_MAX );
        }

      if( LON_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_longitude", LON_MAX );
        }
      
      if( YMIN < 1.0e+200  )
	mag_setr( "contour_min_level", YMIN );

      if( YMAX > -1.0e+200 )
	mag_setr( "contour_max_level", YMAX );

      
      if( COLOUR )
	mag_setc( "contour_line_colour", COLOUR );
      
      
      if( IS_NOT_EQUAL(INTERVAL, 8.0f) )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( THICKNESS != 1 )
	mag_seti( "contour_line_thickness", THICKNESS );
      
      if( STYLE )
      	  mag_setc( "contour_line_style", STYLE );
      
      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      if( DBG )
        fprintf( stderr, " CONTOUR Done!\n" );
    }
  else if ( operatorID == GRFILL )
    {
      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );

      // mag_setc ( "contour_shade_technique", "cell_shading" );
      mag_setc ( "contour_shade_technique", "grid_shading" );

      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );
      
      if( LAT_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_latitude", LAT_MIN );
        }

      if( LON_MIN < 1.0e+200  )
        {
	   mag_setr( "subpage_lower_left_longitude", LON_MIN );
        }

      if( LAT_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_latitude", LAT_MAX );
        }

      if( LON_MAX > -1.0e+200 )
        {
	   mag_setr( "subpage_upper_right_longitude", LON_MAX );
        }
      
      if( YMIN < 1.0e+200  )
        {
	   mag_setr( "contour_shade_min_level", YMIN );
	   mag_setr( "contour_min_level", YMIN );
        }

      if( YMAX > -1.0e+200 )
        {
	   mag_setr( "contour_shade_max_level", YMAX );
	   mag_setr( "contour_max_level", YMAX );
        }

      /*
      if( YMIN < 1.0e+200  )
	mag_setr( "contour_shade_min_level", YMIN );
      
      if( YMAX > -1.0e+200 )
	mag_setr( "contour_shade_max_level", YMAX );
      */
      
      if( COLOUR_MIN )
	mag_setc( "contour_shade_min_level_colour", COLOUR_MIN );
      
      if( COLOUR_MAX )
	mag_setc( "contour_shade_max_level_colour", COLOUR_MAX );
      
      if( IS_NOT_EQUAL(INTERVAL, 8.0f) )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( USR_COLOUR_COUNT )
	{
	  mag_setc( "contour_shade_colour_method", "LIST" );
	  mag_set1c( "contour_shade_colour_list",( const char ** ) USR_COLOUR_TABLE, USR_COLOUR_COUNT ); 
	}
      /*
      if( IS_NOT_EQUAL(RESOLUTION, 10.0f) )
	mag_setr( "contour_shade_cell_resolution", RESOLUTION );
      */
      if( COLOUR_TRIAD )                                
	  mag_setc( "contour_shade_colour_direction", COLOUR_TRIAD );

      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

      if( DBG )
        fprintf( stderr, " GrFILL Done!\n");
    }

  /* plot the title text and the coastlines */
  mag_cont ();
  mag_coast ();


  mag_set1c("text_lines", (const char **) &titlename, 1);
  mag_setc("text_colour", "black");

/*
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
*/

  mag_setc("text_justification", "left");
  mag_text();

  if ( LEV_LIST ) Free(LEV_LIST);
}


static
void init_MAGICS( )
{
  setenv( "MAGPLUS_QUIET","1",1 ); /* To suppress magics messages */

  mag_open();
/* Some standard parameters affectng the magics environment, moved from the xml file  ** begin ** */
  mag_setc ("page_id_line","off");
  mag_setc(  "output_name_first_page_number", "off" );
  if( FILE_SPLIT == TRUE )
    mag_setc(  "output_ps_split" , "on" );
}

static
void quit_MAGICS( )
{
  mag_close ();
  if( DBG )
    fprintf( stderr,"Exiting From MAGICS\n" );
}

static
void VerifyPlotParameters( int num_param, char **param_names, int opID )
{
  int i, j, k;
  int found = FALSE, syntax = TRUE, halt_flag = FALSE, /* file_found = TRUE, */ split_str_count;
  int param_count = 0;
  const char **params = NULL;
  char **split_str = NULL, **split_str1 = NULL;
  const char *sep_char = "=";
  char *temp_str;
  const char  orig_char = ';', rep_char = ',';
  FILE *fp;

/*  
  char  *contour_params[] = {"ymin","ymax","count","interval","list","colour","thickness","style"};
  char  *shaded_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colour_table","step_freq"};
  char  *grfill_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colour_table","resolution"};
*/


  for ( i = 0; i < num_param; ++i )
    {
      split_str_count = 0;
      found = FALSE;
      syntax = TRUE;
      split_str_count = StringSplitWithSeperator( param_names[i], sep_char, &split_str );
      
      if( DBG )
	fprintf( stderr, "Verifying params!\n");
      
      if( split_str_count > 1 ) 
	{
	  
	  if( opID == CONTOUR )
	    {
	      param_count = contour_param_count;
	      params = contour_params;
	    }
	  else if( opID == SHADED )
	    {
	      param_count = shaded_param_count;
	      params = shaded_params;
	    }
	  else if( opID == GRFILL )
	    {
	      param_count = grfill_param_count;
	      params = grfill_params;
	    }
	  
	  for ( j = 0; j < param_count; ++j )
	    {
	      if( !strcmp( split_str[0], params[j] ) )
		{
		  found = TRUE;
		  if( !strcmp( split_str[0],"colour" )     || !strcmp( split_str[0],"style" )       ||
		      !strcmp( split_str[0],"colour_min" ) || !strcmp( split_str[0],"colour_max" )  ||
		      !strcmp( split_str[0],"RGB" )        || !strcmp( split_str[0],"colour_triad" )||
		      !strcmp( split_str[0],"device")      || !strcmp( split_str[0],"file_split" )  ||
		      !strcmp( split_str[0],"projection") 
		    )
		    {
		      if( IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else
			{
			  if( !strcmp( split_str[0],"RGB" ) || !strcmp( split_str[0],"file_split") )
			    {
			      temp_str = strdup( split_str[1] );    
			      StrToUpperCase( temp_str );
			      if( strcmp( temp_str,"TRUE" ) && strcmp( temp_str,"FALSE" ) )
				syntax = FALSE;      
			      else
				{
                                  if( !strcmp( split_str[0],"RGB" ) )
                                    {
				      if( !strcmp( temp_str,"TRUE" ) )
				        isRGB = TRUE;
				      else
				        isRGB = FALSE;
				    }
                                  else if( !strcmp( split_str[0],"file_split" ) )
                                    {
				      if( !strcmp( temp_str,"TRUE" ) )
				        FILE_SPLIT = TRUE;
				      else
				        FILE_SPLIT = FALSE;
				    }
				}
			    }
			  else if( !strcmp( split_str[0],"style" ) )
			    {
			      if( checkstyle( split_str[1] ) )
				syntax = FALSE;
			    }
			  else if( !strcmp( split_str[0],"colour" ) || !strcmp( split_str[0],"colour_min" ) || !strcmp( split_str[0],"colour_max" ) )
			    {
			      
			      if( checkcolour( split_str[1] ) )
				syntax = FALSE;
                              else
                                {
      				   if( !strcmp( split_str[0],"colour" ) ) 
	                             {  
	                               temp_str = strdup( split_str[1] );  
	                               if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                    StrToUpperCase( temp_str );
	                                    StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
                                       COLOUR = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR );
	                             }
                                   if( !strcmp( split_str[0],"colour_min" ) ) 
	                             {  
	                               temp_str = strdup( split_str[1] );    
                                       if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                   StrToUpperCase( temp_str );
	                                   StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
	                               COLOUR_MIN = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR_MIN );
	                             }
                                   if( !strcmp( split_str[0],"colour_max" ) ) 
	                             {
	                               temp_str = strdup( split_str[1] );    
                                       if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                   StrToUpperCase( temp_str );
	                                   StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
	                               COLOUR_MAX = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR_MAX );
	                             }
			        }
			    }
			  else if( !strcmp( split_str[0],"device" ) )
			    {
			      if( checkdevice( split_str[1] ) )
  				    syntax = FALSE;
			    }
			  else if( !strcmp( split_str[0],"colour_triad" ) )
			    {
			      temp_str = strdup( split_str[1] );    
			      StrToUpperCase( temp_str );
			      if( strcmp( temp_str,"CW" ) && strcmp( temp_str,"ACW" ) )
				syntax = FALSE;      
			      else
				{
				   if( DBG )
				     fprintf( stderr, "TRIAD check  %s!\n",temp_str);
				   if( !strcmp( temp_str,"CW" ) )
				     COLOUR_TRIAD = "clockwise";
				   else
				     COLOUR_TRIAD = "anti_clockwise";
				}
			    }
			  else if( !strcmp( split_str[0],"projection" ) )
			    {
			      if( checkprojection( split_str[1] ) )
				syntax = FALSE;
			    }
			}
		    }
		      
		  if( !strcmp( split_str[0],"min" )       ||  !strcmp( split_str[0],"max" )            ||
		      !strcmp( split_str[0],"lat_min" )   ||  !strcmp( split_str[0],"lat_max" )        ||
		      !strcmp( split_str[0],"lon_min" )   ||  !strcmp( split_str[0],"lon_max" )        ||
		      !strcmp( split_str[0],"count" )     ||  !strcmp( split_str[0],"interval" )       ||
		      !strcmp( split_str[0],"thickness" ) ||  !strcmp( split_str[0],"resolution" )     ||
		      !strcmp( split_str[0],"step_freq" )
                    )
		    {
		      if( !IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else
			{
                           if( !strcmp( split_str[0],"min" ) )
	                     {
	                        YMIN = atof( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"max" ) )
	                     {
	  		        YMAX = atof( split_str[1] );
			     }
		           if( !strcmp( split_str[0],"count" ) )
			     {
	                        COUNT = atoi( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"interval" ) )
	                     {
	  		        INTERVAL = atof( split_str[1] );
	                     }
		           if( !strcmp( split_str[0],"thickness" ) )
			     {
	  		        THICKNESS = atoi( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"resolution" ) )
	                     {
	                        RESOLUTION = atoi( split_str[1] );
	                     }	
		           if( !strcmp( split_str[0],"step_freq" ) )
			     {
	  	                STEP_FREQ = atoi( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"lat_min" ) )
	                     {
	                        LAT_MIN = atof( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"lat_max" ) )
	                     {
	  		        LAT_MAX = atof( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"lon_min" ) )
	                     {
	                        LON_MIN = atof( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"lon_max" ) )
	                     {
	  		        LON_MAX = atof( split_str[1] );
			     }
	                }
		    }
		    
		  if( !strcmp( split_str[0],"colour_table" ) )
		    {
		      if( ( fp = fopen( split_str[1],"r") ) == NULL )
			{
			  fprintf( stderr,"Input Color Table File not found in specified path '%s'\n", split_str[1] );
			  halt_flag = TRUE;
			}
		      else
			{
			  ReadColourTable ( split_str[1] );
			}
		    }
		    
		  if( !strcmp( split_str[0],"list" ) )
		    {
		      sep_char = ";";
		      split_str_count = StringSplitWithSeperator( split_str[1], sep_char, &split_str1 );
		      if( !split_str_count )
			{
			  syntax = FALSE;
			}
		      else
		        {
			  for( k = 0; k < split_str_count; k++ )
			    {
			      if( !IsNumeric( split_str1[k] ) )
			        syntax = FALSE;
			    }
			  if( syntax == TRUE )
			    {
	                       NUM_LEVELS = split_str_count;
	                       LEV_LIST = (double*) Malloc( sizeof( double ) * split_str_count );
			       for( k = 0; k < split_str_count; k++ )
		                 {
		                    LEV_LIST[k] = atof( split_str1[k] );
		                 }
	                       Free(split_str1);
	                    }
			  }
		      }
		  sep_char = "=";
		} /*** if( !strcmp( split_str[0], params[j] ) )  ***/
	    } /*** Loop over param count ***/
	} /*** ( split_str_count > 1 ) ***/
      else
	{
	  syntax = FALSE;
	}
	
      if( found == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter  '%s'\n", param_names[i] );
	} 
      if( found == TRUE && syntax == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter specification  '%s'\n", param_names[i] );
	}
	
      if ( split_str ) 	  
        Free(split_str);
    } /*** Loop over params ****/
      
    if( halt_flag == TRUE )
      {
        exit(0);
      }
    
}


int checkcolour( char *colour_in )
{
    int i, n;
    int split_str_count;
    const char *sep_char =",";
    char **split_str = NULL;
    float  rgb_values[3];
    char temp[256];
    char *ref;
   
    ref = colour_in;
    
    if( isRGB )
      {
	if( strchr( colour_in,';') == NULL || strstr( colour_in,"RGB(") == NULL )
	  {
	    cdoWarning( "Found 'RGB=true',Specify Colour in 'RGB(r;g;b)' ( where r,g,b in [0.0,1.0] ) format!" );
	    Free(split_str);
	    return 1;
	  }
	  
	n = strlen( colour_in );
    
	if( DBG )
	  fprintf( stdout,"  count %d  original colour %s RGB %d\n", n, colour_in, isRGB  );
	
	for( i=0 ; i< n-1; i++ )
	  {
	    if( i > 3 )
	      { 
		temp[i-4] = *colour_in;
	      }
	    colour_in++; 
	  }
	  
	temp[i-4] = '\0';
	
	if( DBG )
	  fprintf( stdout,"  count %d  modified color %s \n", (int)strlen(temp), temp  );
	
	sep_char =";";
	split_str_count = StringSplitWithSeperator( temp, sep_char, &split_str );
    
	if(  split_str_count != 3 ) 
	  {
	    cdoWarning( " Colour specified in Improper format!" );
	    Free(split_str);
	    return 1;
	  }
    
	rgb_values[0] = atof( split_str[0] );
	rgb_values[1] = atof( split_str[1] );
	rgb_values[2] = atof( split_str[2] );
    
	if( rgb_values[0] + rgb_values[1] + rgb_values[2] > 3.0f  || 
	    rgb_values[0] + rgb_values[1] + rgb_values[2] < 0.0f 	   )
	  {
	    cdoWarning( " RGB Colour specified with Improper values!" );
	    Free(split_str);
	    return 1;
	  }
	  
        Free(split_str);  
      }
    else
      {
	if( strchr( colour_in,';') != NULL || strstr( colour_in,"RGB(") != NULL )
	  {
	    cdoWarning( "Found Colour with 'RGB(r;g;b)' format, set parameter RGB='true' !" );
	    Free(split_str);
	    return 1;
	  }
	  
	StrToLowerCase( colour_in );
	for ( i = 0 ; i < STD_COLOUR_COUNT; i++ )
	  {
	    if ( !strcmp( STD_COLOUR_TABLE[i], colour_in ) ) return 0;
	  }
        cdoWarning( "Specified Colour not in Standard colour list, resetting to blue(default colour)!" );
        return 1;
      }
      
    if( DBG )  
      cdoWarning( "Colour %s verified!",ref );

    return 0;
}


int ReadColourTable( char *filepath )
{
  int  i, num_colors = 0;
  char  orig_char = ';', rep_char = ',';
    
  FILE *fp = fopen( filepath,"r" );
    
  if ( !fp )
    {
      fprintf( stdout, "File Not available!" );
      return 1;
    }
    
  fscanf( fp, "%d", &num_colors );
    
  if( DBG )
    fprintf( stderr, "Num Colours %d\n", num_colors );
    
  if ( !num_colors )
    {
      cdoWarning("No colours found in File, proceeding with Standard Colour table!");
      fclose(fp);
      return 1;
    }
    
  USR_COLOUR_COUNT = 0;
  USR_COLOUR_TABLE = ( char **) Malloc( num_colors * sizeof( char* ));
  char **temp_table = ( char **) Malloc( num_colors * sizeof( char* ));
    
  for ( i = 0; i < num_colors; i++ )
    {
      temp_table[i] = (char *) Malloc(256 * sizeof( char ));
      fscanf( fp, "%s", temp_table[i] );
      if( DBG )
        fprintf( stdout, "%s\n", temp_table[i] );
    }
    
  fclose(fp);

  for ( i = 0; i < num_colors; i++ )
    {
      if( DBG )
        fprintf( stdout, "%s \n", temp_table[i] );
	  
      if ( !checkcolour( temp_table[i] ) )
        {
          if( isRGB )
            StrReplaceChar( temp_table[i], orig_char, rep_char ); /* replace ';' in RGB format to ',' */

          if( DBG )
            fprintf( stdout, "Before appending %s\n", temp_table[i] );
          
          USR_COLOUR_TABLE[ USR_COLOUR_COUNT ] = strdup( temp_table[i] );
	      
          /* strcpy( USR_COLOUR_TABLE[ USR_COLOUR_COUNT ], temp_table[i] ); */
          USR_COLOUR_COUNT++;
	      
          if( DBG )
            fprintf( stdout, "After appending %s\n", temp_table[i] );
        }
    }
    
  if( USR_COLOUR_COUNT < num_colors )
    {
      cdoWarning( " Discarding improper format colours and continuing!" );
    }

  for ( i = 0; i < num_colors; i++ ) Free(temp_table[i]);
  Free(temp_table);
    
  return 0;
}


int checkstyle( char *style_in )
{
    int i, found = FALSE;
    StrToUpperCase( style_in );
    for( i = 0 ; i < STYLE_COUNT; i++ )
      {
	if( DBG )
	  fprintf( stderr, "Input %s ref %s\n",style_in, STYLE_TABLE[i] );
	
	if( !strcmp( STYLE_TABLE[i], style_in ) )
	  {
	    found = TRUE;
	    STYLE = style_in;
	    return 0;
	  }
      }
      
    if( !found )
	 cdoWarning( " Style specified with Improper value!" );
    
    return 1; 
}


int checkdevice( char *device_in )
{
    int i, found = FALSE;
    StrToUpperCase( device_in );
    for( i = 0 ; i < DEVICE_COUNT; i++ )
      {
	if( DBG )
	  fprintf( stderr, "Input %s ref %s\n",device_in, DEVICE_TABLE[i] );
	
	if( !strcmp( DEVICE_TABLE[i], device_in ) )
	  {
	    found = TRUE;

	    DEVICE = device_in;
	    if( DBG )
	      fprintf( stderr,"DEVICE %s\n",DEVICE );

	    if( !strcmp( "GIF_ANIMATION" , device_in ) || !strcmp( "KML", device_in )  )
              {
	         ANIM_FLAG = 1;
		 STEP_FREQ = 1;
	      }
	    return 0;
	  }
      }
      
    if( !found )
	 cdoWarning( " Device specified with Improper value!" );
    
    return 1; 
}


int checkprojection( char *projection_in )
{
    int i, found = FALSE;

    /*StrToUpperCase( projection_in );*/

    for( i = 0 ; i < PROJECTION_COUNT; i++ )
      {
	if( DBG )
	  fprintf( stderr, "Input %s ref %s\n",projection_in, PROJECTION_TABLE[i] );
	
	if( !strcmp( PROJECTION_TABLE[i], projection_in ) )
	  {
	    found = TRUE;
	    PROJECTION = projection_in;
	    return 0;
	  }
      }
      
    if( !found )
      {	
	 cdoWarning( " Projection specified with Improper value!" );
	 cdoWarning( " Specify one of the following:" );
	 cdoWarning( " cylindrical polar_stereographic polar_north geos meteosat meteosat_57E geos_east lambert EPSG3857 goode collignon mollweide robinson bonne google efas EPSG4326 lambert_north_atlantic mercator cartesian taylor tephigram" );

      }	
    
    return 1; 
}
#endif


void *Magplot(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBMAGICS)
  int nrecs;
  int levelID;
  int nmiss;
  char varname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32], datetimestr[64];
  
  int nparam = operatorArgc();
  char **pnames = operatorArgv();
  
  CONTOUR = cdoOperatorAdd("contour", 0, 0, NULL);
  SHADED  = cdoOperatorAdd("shaded", 0, 0, NULL);
  GRFILL  = cdoOperatorAdd("grfill", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  if( nparam )
    {
      if( DBG )
	{
	  for( int i = 0; i < nparam; i++ )
	    fprintf( stderr,"Param %d is %s!\n",i+1, pnames[i] );
	}
      
      VerifyPlotParameters( nparam, pnames, operatorID );
    }

  int streamID = pstreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  int varID = 0;
  int gridID  = vlistInqVarGrid(vlistID, varID);
  // int zaxisID = vlistInqVarZaxis(vlistID, varID);

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_GME          ) cdoAbort("GME grid unspported!");
  if ( gridtype == GRID_UNSTRUCTURED ) cdoAbort("Unstructured grid unspported!");

  bool lregular = false;
  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) lregular = true;
  
  if ( gridtype != GRID_CURVILINEAR ) gridID = gridToCurvilinear(gridID, 1);

  int gridsize = gridInqSize(gridID);
  int nlon     = gridInqXsize(gridID);
  int nlat     = gridInqYsize(gridID);
  //int nlev     = zaxisInqSize(zaxisID);

  double *array           = (double*) Malloc(gridsize*sizeof(double));
  double *grid_center_lon = (double*) Malloc(gridsize*sizeof(double));
  double *grid_center_lat = (double*) Malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, grid_center_lon);
  gridInqYvals(gridID, grid_center_lat);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");
					
  int tsID = 0;

  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  /*
     init_XMLtemplate_parser( Filename );
     updatemagics_and_results_nodes( );
  */


  init_MAGICS( );

  while ( (nrecs = pstreamInqTimestep(streamID, tsID)) )
    {
      if( ANIM_FLAG )
        {
      	  if( nrecs > 1 )
	    {
	      cdoWarning("File has more than one variable! Animation creation not possible!!!");
	      break;
            }
      	  if( tsID % STEP_FREQ )
	    {
                tsID++;
		continue;
            }
	}
      else 	
        {
          if( STEP_FREQ )
	    {
          	if( tsID % STEP_FREQ )
	    	  {
                     tsID++;
	             cdoWarning("NOT PLOTTING STEP %d!!!",tsID);
	             continue;
	    	  }
            }
         else 
            {
		if( tsID )
		  {
	   		cdoWarning("File variables have values at more than one time step! Images created for first time step!!!");
           		cdoWarning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!!!");
           		cdoWarning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection operators!!!");
       	   		break;
		  }
	    }
        }
      
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      sprintf( datetimestr, "%s %s", vdatestr, vtimestr );
      if( DBG )
        fprintf( stderr,"Date %s Time %s\n",vdatestr, vtimestr );

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID, &varID, &levelID);
	  pstreamReadRecord(streamID, array, &nmiss);

          if ( nmiss ) cdoSetNAN(vlistInqVarMissval(vlistID, varID), gridsize, array);
	  
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarUnits(vlistID, varID, units);

	  if ( operatorID == SHADED || operatorID == CONTOUR || operatorID == GRFILL )
          {
                if( DBG )
                  {
                     if( operatorID == SHADED )
                       fprintf( stderr," Creating SHADED PLOT for %s\n",varname );
                     else if( operatorID == CONTOUR )
                       fprintf( stderr," Creating CONTOUR PLOT for %s\n",varname );
                     else if( operatorID == GRFILL )
                       fprintf( stderr," Creating GRFILL PLOT for %s\n",varname );
                  }

                if( DBG )
                  fprintf( stderr,"Plot %d\n",varID );
	  	magplot(cdoStreamName(1)->args, operatorID, varname, units, nlon, nlat, grid_center_lon, grid_center_lat, array, nparam, pnames, datetimestr, lregular);
          }
	  else
	  	fprintf(stderr,"operator not implemented\n");
	}

      if( DBG )
        fprintf( stderr,"TimeStep %d\n",tsID );

       
      tsID++;
      /*
      if( !STEP_FREQ  && tsID )
        {
	   cdoWarning("File variables have values at more than one time step! Images created for first time step!!!");
           cdoWarning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!!!");
           cdoWarning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection operators!!!");
       	   break;
	}
      else
        {
      	   tsID++;
           if( DBG )
             fprintf( stderr,"TimeStep %d\n",tsID );
	}
      */
    }

  if( ANIM_FLAG )
    {
      if( FILE_SPLIT == TRUE  ) 
        cdoWarning("File split parameter ignored!!!");
    }
  quit_MAGICS( );

  pstreamClose(streamID);

  if ( array  ) Free(array);
  if ( grid_center_lon ) Free(grid_center_lon);
  if ( grid_center_lat ) Free(grid_center_lat);

/*   quit_XMLtemplate_parser( ); */
#else
  
  cdoAbort("MAGICS support not compiled in!");

#endif

  cdoFinish();

  return 0;

}
