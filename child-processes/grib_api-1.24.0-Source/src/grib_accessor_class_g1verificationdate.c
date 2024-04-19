/*
 * Copyright 2005-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
 * virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
 */

#include "grib_api_internal.h"

/* 
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_long
   IMPLEMENTS = unpack_long
   IMPLEMENTS = init;dump
   MEMBERS=const char* date
   MEMBERS=const char* time
   MEMBERS=const char* step
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int unpack_long(grib_accessor*, long* val,size_t *len);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_g1verificationdate {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in g1verificationdate */
	const char* date;
	const char* time;
	const char* step;
} grib_accessor_g1verificationdate;

extern grib_accessor_class* grib_accessor_class_long;

static grib_accessor_class _grib_accessor_class_g1verificationdate = {
    &grib_accessor_class_long,                      /* super                     */
    "g1verificationdate",                      /* name                      */
    sizeof(grib_accessor_g1verificationdate),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    &dump,                       /* describes himself         */
    0,                /* get length of section     */
    0,              /* get length of string      */
    0,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    &unpack_long,                /* grib_unpack procedures long    */
    0,                /* grib_pack procedures double    */
    0,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    0,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    0,            /* pack_expression */
    0,              /* notify_change   */
    0,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_g1verificationdate = &_grib_accessor_class_g1verificationdate;


static void init_class(grib_accessor_class* c)
{
	c->next_offset	=	(*(c->super))->next_offset;
	c->string_length	=	(*(c->super))->string_length;
	c->value_count	=	(*(c->super))->value_count;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->pack_double	=	(*(c->super))->pack_double;
	c->unpack_double	=	(*(c->super))->unpack_double;
	c->pack_string	=	(*(c->super))->pack_string;
	c->unpack_string	=	(*(c->super))->unpack_string;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
	c->update_size	=	(*(c->super))->update_size;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->compare	=	(*(c->super))->compare;
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long l, grib_arguments* c)
{
	grib_accessor_g1verificationdate* self = (grib_accessor_g1verificationdate*)a;
	int n = 0;

	self->date = grib_arguments_get_name(a->parent->h,c,n++);
	self->time = grib_arguments_get_name(a->parent->h,c,n++);
	self->step = grib_arguments_get_name(a->parent->h,c,n++);

	a->flags |= GRIB_ACCESSOR_FLAG_READ_ONLY; 
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
	grib_dump_long(dumper,a,NULL);
}


static int unpack_long(grib_accessor* a, long* val, size_t *len)
{   
	grib_accessor_g1verificationdate* self = (grib_accessor_g1verificationdate*)a;
  int ret=0;
	long date = 0;
	long time = 0;
	long cdate = 0;
	long step = 0;
	long vtime = 0;
	long vdate = 0;
	long vd = 0;

	if ((ret=grib_get_long_internal(a->parent->h, self->date,&date))!=GRIB_SUCCESS) return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->time,&time))!=GRIB_SUCCESS) return ret;
  if ((ret=grib_get_long_internal(a->parent->h, self->step,&step))!=GRIB_SUCCESS) return ret;

  time /= 100;

	cdate = (long)grib_date_to_julian (date);
	vtime = cdate * 24 + time + step;
	vd    = vtime / 24;
	vdate = grib_julian_to_date (vd);

	  /* printf("\n********\n date %d, time %d, step %d, vdate: %d, cdate %d, vd %d\n********\n", date, time, step, vdate, cdate, vd);    */

	if(*len < 1)
		return GRIB_ARRAY_TOO_SMALL;

	*val   = vdate;

	  /* fprintf(stdout,"\n********\n %d cdate %d vd %d\n********\n", vdate, cdate, step); */
	return GRIB_SUCCESS;
}

