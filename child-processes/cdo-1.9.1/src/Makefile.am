## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libcdo.la
libcdo_la_SOURCES =            \
			   argument.h		\
			   argument.cc		\
			   array.h         \
               array.cc         \
               cdo_int.h       \
               compare.h       \
               cdo_pthread.cc   \
               cdo_vlist.cc     \
               cdo_getopt.cc    \
               cdo_getopt.h    \
               cdo_task.cc      \
               cdo_task.h      \
               cdo_history.cc   \
               cdo_read.cc      \
               cdi_uuid.h      \
               cmortable_parser.cc   \
               after_sptrans.cc \
               after_fctrans.cc \
               after_dvtrans.cc \
               after_vertint.cc \
               after_vertint.h \
               after_namelist.cc\
               afterburnerlib.cc\
               afterburner.h   \
               vct_l191.h      \
               calendar.h      \
               constants.h     \
               constants.cc     \
               color.cc         \
               color.h         \
               commandline.cc   \
               const.h         \
               convert_units.cc \
               convert_units.h \
               counter.h       \
               datetime.cc      \
               datetime.h      \
               dmemory.h       \
               ecacore.cc       \
               ecacore.h       \
               ecautil.cc       \
               ecautil.h       \
               error.h         \
               etopo.h         \
               temp.h          \
               mask.h          \
               exception.cc     \
               expr.cc          \
               expr.h          \
               expr_fun.cc      \
               expr_fun.h      \
               expr_lex.cc      \
               expr_yacc.cc     \
               expr_yacc.h     \
               features.cc      \
               field.cc         \
               field.h         \
               field2.cc        \
               fieldc.cc        \
               fieldmem.cc      \
               fieldmer.cc      \
               fieldzon.cc      \
               functs.h        \
               gradsdeslib.cc   \
               gradsdeslib.h   \
               grid.cc          \
               grid.h          \
               grid_proj.cc     \
               grid_proj.h     \
               grid_area.cc     \
               grid_define.cc   \
               grid_gme.cc      \
               grid_rot.cc      \
               grid_from_name.cc   \
               grid_read.cc        \
               grid_read_pingo.cc  \
               grid_print.cc    \
               gridreference.cc \
               griddes.cc       \
               griddes.h       \
               griddes_h5.cc    \
               griddes_nc.cc    \
               hetaeta.cc       \
               hetaeta.h       \
               institution.cc   \
               interpol.cc      \
               interpol.h      \
               job.cc           \
               juldate.cc       \
               grid_search.cc   \
               grid_search.h   \
               listarray.cc     \
               listarray.h     \
               list.cc          \
               list.h          \
               listbuf.cc       \
               listbuf.h       \
	       merge_sort2.cc   \
	       merge_sort2.h   \
               modules.cc       \
               modules.h       \
               namelist.cc      \
               namelist.h      \
               namelist_parser.cc     \
               normal.cc        \
               nth_element.cc   \
               nth_element.h   \
               operator_help.h \
               par_io.cc        \
               par_io.h        \
               parse_literal.cc      \
               percentiles_hist.cc   \
               percentiles_hist.h   \
               percentiles.cc   \
               percentiles.h   \
               pipe.cc          \
               pipe.h          \
               pmlist.cc        \
               pmlist.h        \
               sellist.cc       \
               sellist.h       \
               pragma_omp_atomic_update.h \
               printinfo.h     \
               process.cc       \
               process.h       \
               pstream.cc       \
               pstream.h       \
               pstream_write.h \
               pstream_int.h   \
               pthread_debug.cc \
               pthread_debug.h \
               readline.cc      \
               realtime.cc      \
               remap.h         \
               remaplib.cc      \
               remapsort.cc     \
               remap_scrip_io.cc \
               remap_search_reg2d.cc \
               remap_search_latbins.cc \
               remap_store_link.cc \
               remap_store_link.h \
               remap_store_link_cnsrv.cc \
               remap_store_link_cnsrv.h \
               remap_conserv.cc \
               remap_conserv_scrip.cc \
               remap_distwgt.cc \
               remap_bicubic_scrip.cc \
               remap_bilinear_scrip.cc \
               stdnametable.cc  \
               stdnametable.h  \
               specspace.cc     \
               specspace.h     \
               statistic.cc     \
               statistic.h     \
               table.cc         \
               text.cc          \
               text.h          \
               timebase.h      \
               timer.cc         \
               userlog.cc       \
               uthash.h        \
               util.cc          \
               util.h          \
			   varlist.h		\
               zaxis_print.cc   \
               zaxis.cc
libcdo_la_SOURCES +=           \
               json/jsmn.h     \
               json/jsmn.c
libcdo_la_SOURCES +=           \
               kdtreelib/kdtree.h             \
               kdtreelib/kdtree_cartesian.cc  \
               kdtreelib/kdtree_common.cc     \
               kdtreelib/kdtree_spherical.cc  \
               kdtreelib/qsort.cc             \
               kdtreelib/pmergesort.cc        \
               kdtreelib/pqueue.cc            \
               kdtreelib/pqueue.h
libcdo_la_SOURCES +=           \
               clipping/clipping.c           \
               clipping/clipping.h           \
               clipping/area.c               \
               clipping/area.h               \
               clipping/ensure_array_size.c  \
               clipping/ensure_array_size.h  \
               clipping/geometry.h           \
               clipping/grid.h               \
               clipping/points.h             \
               clipping/dep_list.h           \
               clipping/grid_cell.c          \
               clipping/grid_cell.h          \
               clipping/intersection.c       \
               clipping/utils.c              \
               clipping/utils.h 
#
bin_PROGRAMS = cdo
#
cdo_SOURCES  = cdo.cc
cdo_SOURCES += Adisit.cc        \
               Afterburner.cc   \
               Arith.cc         \
               Arithc.cc        \
               Arithdays.cc     \
               Arithlat.cc      \
               CDItest.cc       \
               CDIread.cc       \
               CDIwrite.cc      \
               Cat.cc           \
               Change.cc        \
               Change_e5slm.cc  \
               Cloudlayer.cc    \
               CMOR.cc          \
               CMOR_lite.cc     \
               CMOR_table.cc    \
               Collgrid.cc      \
               Command.cc       \
               Comp.cc          \
               Compc.cc         \
               Complextorect.cc \
               Cond.cc          \
               Cond2.cc         \
               Condc.cc         \
               Consecstat.cc    \
               Copy.cc          \
               Deltat.cc        \
               Deltime.cc       \
               Derivepar.cc     \
               Detrend.cc       \
               Diff.cc          \
               Distgrid.cc      \
               Duplicate.cc     \
               EOFs.cc          \
	       Eof3d.cc         \
               EcaIndices.cc    \
               Echam5ini.cc     \
               Enlarge.cc       \
               Enlargegrid.cc   \
               Ensstat.cc       \
               Ensstat3.cc      \
               Ensval.cc        \
               Eofcoeff.cc      \
	       Eofcoeff3d.cc    \
               EstFreq.cc       \
               Exprf.cc         \
               FC.cc            \
               Filedes.cc       \
               Fillmiss.cc      \
               Filter.cc        \
               Fldrms.cc        \
               Fldstat.cc       \
               Fldstat2.cc      \
               Fourier.cc       \
               Gengrid.cc       \
               Gradsdes.cc      \
               Gridboxstat.cc   \
               Gridcell.cc      \
               Gridsearch.cc    \
               Harmonic.cc      \
               Hi.cc            \
               Histogram.cc     \
               Importamsr.cc    \
               Importbinary.cc  \
               Importcmsaf.cc   \
               Importobs.cc     \
               Info.cc          \
               Input.cc         \
               Intgrid.cc       \
               Intgridtraj.cc   \
               Intlevel.cc      \
               Intlevel3d.cc    \
               Intntime.cc      \
               Inttime.cc       \
               Intyear.cc       \
               Invert.cc        \
               Invertlev.cc     \
               Isosurface.cc    \
               Log.cc           \
	       MapReduce.cc     \
               Maskbox.cc       \
               Mastrfu.cc       \
               Math.cc          \
               Merge.cc         \
               Mergegrid.cc     \
               Mergetime.cc     \
               Merstat.cc       \
               Monarith.cc      \
               Mrotuv.cc        \
               Mrotuvb.cc       \
               Ninfo.cc         \
               Nmldump.cc       \
               Output.cc        \
               Outputgmt.cc     \
               Pack.cc          \
               Pardup.cc        \
               Pinfo.cc         \
               Pressure.cc      \
               Regres.cc        \
               Remap.cc         \
               Remapeta.cc      \
               Replace.cc       \
               Replacevalues.cc \
               Rhopot.cc        \
               Rotuv.cc         \
               Runpctl.cc       \
               Runstat.cc       \
               Samplegrid.cc    \
               Samplegridicon.cc \
               Seascount.cc     \
               Seaspctl.cc      \
               Seasstat.cc      \
               Selbox.cc        \
               Selgridcell.cc   \
               Select.cc        \
               Selmulti.cc      \
               Seloperator.cc   \
               Selrec.cc        \
               Seltime.cc       \
               Selvar.cc        \
               Set.cc           \
               Setattribute.cc  \
               Setbox.cc        \
               Setgatt.cc       \
               Setgrid.cc       \
               Sethalo.cc       \
               Setmiss.cc       \
               Setpartab.cc     \
               Setrcaname.cc    \
               Settime.cc       \
               Setzaxis.cc      \
               Shiftxy.cc       \
               Showinfo.cc      \
               Showattribute.h  \
               Showattribute.cc \
               Sinfo.cc         \
               Smooth.cc        \
               Sort.cc          \
               Sorttimestamp.cc \
               Specinfo.cc      \
               Spectral.cc      \
               Spectrum.cc      \
               Split.cc         \
               Splitrec.cc      \
               Splitsel.cc      \
               Splittime.cc     \
               Splityear.cc     \
               Subtrend.cc      \
               Tee.cc           \
               Templates.cc     \
               Test.cc          \
               Tests.cc         \
               Timcount.cc      \
               Timcumsum.cc     \
               Timpctl.cc       \
               Timselpctl.cc    \
               Timselstat.cc    \
               Timsort.cc       \
               Timstat.cc       \
               Timstat2.cc      \
               Timstat3.cc      \
               Tinfo.cc         \
               Tocomplex.cc     \
               Transpose.cc     \
               Trend.cc         \
               Trms.cc          \
               Tstepcount.cc    \
               Vargen.cc        \
               Varrms.cc        \
               Vertintml.cc     \
               Vertintap.cc     \
               Vertstat.cc      \
               Vertcum.cc       \
               Vertwind.cc      \
               Verifygrid.cc    \
               Wct.cc           \
               Wind.cc          \
               WindTrans.cc     \
               Writegrid.cc     \
               Writerandom.cc   \
               XTimstat.cc      \
               YAR.cc           \
               Yearmonstat.cc   \
               Ydayarith.cc     \
               Ydaypctl.cc      \
               Ydaystat.cc      \
               Ydrunpctl.cc     \
               Ydrunstat.cc     \
               Yhourarith.cc    \
               Yhourstat.cc     \
               Ymonarith.cc     \
               Ymonpctl.cc      \
               Ymonstat.cc      \
               Yseaspctl.cc     \
               Yseasstat.cc     \
               Zonstat.cc       \
               cdo.h

cdo_SOURCES += nearpt3c.h
#if ENABLE_NEARPT3
#cdo_SOURCES +=           \
#               nearpt3x.h      \
#               nearpt3c.h      \
#               nearpt3c.cc     \
#               cellsearchorder.h
#endif

#if ENABLE_MAGICS
cdo_SOURCES += Magplot.cc       \
               Magvector.cc     \
               Maggraph.cc      \
               template_parser.h   \
               template_parser.cc   \
               results_template_parser.h   \
               results_template_parser.cc   \
               magics_template_parser.h   \
               magics_template_parser.cc   \
               StringUtilities.h   \
               StringUtilities.cc   \
               CdoMagicsMapper.h   \
               CdoMagicsMapper.cc
#endif

cdo_CPPFLAGS = -I$(top_srcdir)/libcdi/src
cdo_LDADD    = libcdo.la $(top_builddir)/libcdi/src/libcdi.la
cdo_LDFLAGS  =

if ENABLE_ALL_STATIC
cdo_LDFLAGS  += -all-static
endif

libcdo_la_CPPFLAGS = $(cdo_CPPFLAGS)
noinst_PROGRAMS  = cdotest
cdotest_SOURCES  = cdo_int.h	\
	           cdotest.cc
cdotest_LDADD    = $(cdo_LDADD)
cdotest_CPPFLAGS = $(cdo_CPPFLAGS)
cdotest_LDFLAGS  = $(cdo_LDFLAGS)
if ENABLE_ALL_STATIC
cdotest_LDFLAGS  += -all-static
endif

#cdo-userlog.o: userlog.cc config.h
#	$(COMPILE) -DLOGPATH=${exec_prefix}/log -c -o cdo-userlog.o `test -f 'userlog.cc' || echo '$(srcdir)/'`userlog.cc
#cdo_static-userlog.o: userlog.cc config.h
#	$(COMPILE) -DLOGPATH=${exec_prefix}/log -c -o cdo_static-userlog.o `test -f 'userlog.cc' || echo '$(srcdir)/'`userlog.cc
#
CLEANFILES  = `ls *~ 2> /dev/null`
#
clean-local: clean-local-dirs
.PHONY: clean-local-dirs
clean-local-dirs:
	-rm -rf *.dSYM
