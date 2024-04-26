## Process this file with automake to produce Makefile.in
#
EXTRA_DIST = cdi_read_f2003.f90 cdi_write_f2003.f90
#
noinst_PROGRAMS = cdi_write cdi_write_relativ cdi_write_ens cdi_write_hybrid cdi_read cdi_copy
if CREATE_ISOC
  noinst_PROGRAMS += cdi_read_f2003 cdi_write_f2003
endif
if ENABLE_NETCDF
  noinst_PROGRAMS += cdi_write_const
endif
#
AM_CPPFLAGS              = -I$(top_srcdir)/src
AM_FCFLAGS               = $(FPP_INCOPT)$(top_srcdir)/src $(FC_MOD_FLAG)../src
LDADD = $(top_builddir)/src/libcdi.la
#
cdi_write_SOURCES        = cdi_write.c
#
cdi_write_relativ_SOURCES = cdi_write_relativ.c
#
cdi_write_ens_SOURCES    = cdi_write_ens.c
#
cdi_write_const_SOURCES  = cdi_write_const.c
#
cdi_write_hybrid_SOURCES = cdi_write_hybrid.c
#
cdi_read_SOURCES         = cdi_read.c
#
cdi_copy_SOURCES         = cdi_copy.c
#
cdi_read_f2003_SOURCES   = cdi_read_f2003.f90
cdi_read_f2003_LDADD     = $(top_builddir)/src/mo_cdi.o $(top_builddir)/src/libcdi.la
cdi_read_f2003.$(OBJEXT): $(top_builddir)/src/mo_cdi.$(FCMODEXT)
#
cdi_write_f2003_SOURCES  = cdi_write_f2003.f90
cdi_write_f2003_LDADD    = $(top_builddir)/src/mo_cdi.$(OBJEXT) $(top_builddir)/src/libcdi.la
cdi_write_f2003.$(OBJEXT): $(top_builddir)/src/mo_cdi.$(FCMODEXT)

#
if ENABLE_ALL_STATIC
cdi_write_LDFLAGS        = -all-static
cdi_read_LDFLAGS         = -all-static
cdi_copy_LDFLAGS         = -all-static
cdi_read_f2003_LDFLAGS   = -all-static
cdi_write_f2003_LDFLAGS  = -all-static
endif
#
clean-local: clean-local-dirs
.PHONY: clean-local-dirs
clean-local-dirs:
	-rm -rf *.dSYM
