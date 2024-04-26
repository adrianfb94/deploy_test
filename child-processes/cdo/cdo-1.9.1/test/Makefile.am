CLEANFILES =

export

# COLOR ALWAYS also if no color capable terminal detected
#     this is needed i.e. to run the tests via a oneline ssh cmd:
#     ssh host '... make check'
AM_COLOR_TESTS=always

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/config/tap-driver.sh

# tests which should pass
TESTS = threads.test tsformat.test wildcard.test File.test Comp.test Compc.test \
        Read_grib.test Read_netcdf.test Copy_netcdf.test Cat.test Gridarea.test  \
        Detrend.test Genweights.test Remap.test Remap2.test Remapeta.test Select.test Spectral.test \
        Ensstat.test Enspctl.test Gridboxstat.test \
        Vertstat.test Fldstat.test Fldpctl.test Merstat.test Zonstat.test \
        Timstat.test Timselstat.test Seasstat.test \
        Runstat.test Multiyearstat.test Ydrunstat.test \
        Mergetime.test EOF.test \
        Vertint.test Afterburner.test Arithc.test Arith.test Expr.test \
        Gradsdes.test Collgrid.test MapReduce.test Ninfo.test

# tests which should fail
XFAIL_TESTS = 

#        $(top_srcdir)/test/test_Remap.sh \
#        $(top_srcdir)/test/test_info.py
#        $(top_srcdir)/test/test_diff.py
#        $(top_srcdir)/test/test_Arith.py
#        $(top_srcdir)/test/test_Arithc.py
#        $(top_srcdir)/test/test_Selvar.py
#        $(top_srcdir)/test/test_intgridbil.py

#EXTRA_DIST = $(TESTS) $(top_srcdir)/test/testStreams.py

CDO        = $(top_builddir)/src/cdo

DATAPATH   = $(top_srcdir)/test/data

PYTHONPATH = $(top_srcdir)/contrib/python:$(top_srcdir)/test

#.PHONY: $(TESTS)

CLEANFILES += `ls *.pyc`

#AUTOMAKE_OPTIONS = color-tests
