## Process this file with automake to produce Makefile.in
#
SUBDIRS = src interfaces app examples examples/pio tests
#
EXTRA_DIST=config/default \
	doc/cdi_cman.pdf \
	doc/cdi_fman.pdf \
	src/cfortran.doc \
	util/sunf95preproc-wrapper \
	util/sxpreproc-wrapper \
	util/xlfpreproc-wrapper \
	config/interface.rb
#
ACLOCAL_AMFLAGS = -I m4
#
CLEANFILES  = `ls *~`
#
doc/cdi_cman.pdf: 
	(cd $(top_srcdir)/doc/tex ; ./makepdf_c ; mv cdi_cman.pdf .. ; ./cleanup)
#
doc/cdi_fman.pdf: 
	(cd $(top_srcdir)/doc/tex ; ./makepdf_f ; mv cdi_fman.pdf .. ; ./cleanup)
