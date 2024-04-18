# Process this file with autoconf to produce a configure script.
#
#  automake 1.14
#  autoconf 2.68
#  libtool  2.4.2

AC_INIT([cdo], [1.9.1], [http://mpimet.mpg.de/cdo])

AC_DEFINE_UNQUOTED(CDO, ["$PACKAGE_VERSION"], [CDO version])

CONFIG_ABORT=yes
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall])
AC_CONFIG_HEADERS([src/config.h])
AM_MAINTAINER_MODE([disable])

# Set up libtool.
AM_PROG_AR
#
AC_MSG_NOTICE([setting up libtool])
LT_INIT

# bash
AC_PATH_TOOL([BASH],[bash],[:])

# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_CHECK_TOOL([CXX],[g++],[:])
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_CXX
AC_C_RESTRICT
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([CPP],[cpp],[:])
AC_CHECK_TOOL([LD],[ld],[:])
AC_CHECK_TOOL([NM],[nm],[:])
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([AS],[as],[:])
AC_CHECK_TOOL([DLLTOOL],[dlltool],[:])
AC_CHECK_TOOL([OBJDUMP],[objdump],[:])
AC_CHECK_TOOL([STRIP],[strip],[:])
AC_CHECK_TOOL([RANLIB],[ranlib],[:])
AC_PROG_INSTALL
#
FC=no
export FC
#  ----------------------------------------------------------------------
# Check large file support on 32 bit system
AC_SYS_LARGEFILE
#  ----------------------------------------------------------------------
# Checks for library functions.
AC_FUNC_MMAP
#  ----------------------------------------------------------------------
# Checks for library malloc.
AC_CHECK_LIB(malloc, malloc)
#  ----------------------------------------------------------------------
# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])
#  ----------------------------------------------------------------------
# Checks for header files
AC_CHECK_HEADERS(sys/resource.h)
AC_CHECK_HEADERS(sys/times.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_HEADERS(glob.h)
AC_CHECK_HEADERS(fnmatch.h)
AC_CHECK_HEADERS(wordexp.h)
AC_CHECK_HEADERS(execinfo.h)
#  ----------------------------------------------------------------------
# Checks for the availability of functions
AC_CHECK_FUNCS(mallinfo)
AC_CHECK_FUNCS(backtrace)
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(getrlimit)
AC_CHECK_FUNCS(gethostname)

# Check compiler version
case "$CXX" in
  pgc*)   CXX_VERSION=`$CXX -V | head -2 | tail -n 1`;;
  *gcc*)  CXX_VERSION=`$CXX --version | head -n 1`;;
  g++*)   CXX_VERSION=`$CXX --version | head -n 1`;;
  clang*) CXX_VERSION=`$CXX --version | head -n 1`;;
  icpc*)  CXX_VERSION=`$CXX --version | head -n 1`;;
  sxc*)   CXX_VERSION=`$CXX -V 2>&1   | tail -n 1`;;
  xlc*)   CXX_VERSION=`$CXX -qversion 2>&1   | head -n 1`;;
  *)      CXX_VERSION=`$CXX -V 2>&1   | head -n 1 | grep -v error`;;
esac

case "$CC" in
  pgc*)   C_VERSION=`$CC -V | head -2 | tail -n 1`;;
  *gcc*)  C_VERSION=`$CC --version | head -n 1`;;
  g++*)   C_VERSION=`$CC --version | head -n 1`;;
  clang*) C_VERSION=`$CC --version | head -n 1`;;
  icc*)   C_VERSION=`$CC --version | head -n 1`;;
  sxc*)   C_VERSION=`$CC -V 2>&1   | tail -n 1`;;
  xlc*)   C_VERSION=`$CC -qversion 2>&1   | head -n 1`;;
  *)      C_VERSION=`$CC -V 2>&1   | head -n 1 | grep -v error`;;
esac

if test -z "$CXX_VERSION" ; then CXX_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(CXX_VERSION, ["$CXX_VERSION"], [CXX Compiler version])

if test -z "$C_VERSION" ; then C_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(C_VERSION, ["$C_VERSION"], [C Compiler version])

# Checks for username, hostname and system type
USERNAME=$LOGNAME
if test -z "$USERNAME" ; then USERNAME=$USER; fi;
if test -z "$USERNAME" ; then USERNAME="unknown"; fi;
AC_DEFINE_UNQUOTED([USER_NAME],["$USERNAME"], [User name])
AC_SUBST([USER_NAME],["$USERNAME"])

AS_IF([test -z "$HOST"],
      [HOST=unknown
       AS_IF([test -x /bin/hostname],[HOST=$(hostname)],
             [AS_IF([test -x /bin/uname],
                    [HOST=$(uname -n)])])])
AC_DEFINE_UNQUOTED([HOST_NAME],["$HOST"],[Host name])
AC_SUBST([HOST_NAME],["$HOST"])

AC_DEFINE_UNQUOTED([SYSTEM_TYPE],["$ac_cv_build"], [System type])
AC_SUBST([SYSTEM_TYPE],["$ac_cv_build"])
#  ----------------------------------------------------------------------
#  Check for math library
AC_CHECK_LIB([m],[floor])
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 math functions
AC_CHECK_DECLS([isnan],,,[AC_INCLUDES_DEFAULT
@%:@include <math.h>])
AC_CHECK_FUNCS(sqrtl)
AC_CHECK_FUNCS(feenableexcept)
#
AC_CHECK_MEMBERS([fenv_t.__control, fenv_t.__mxcsr],,,[[#include <fenv.h>]])
#  ----------------------------------------------------------------------
#  Enable NEARPT3 support
# AC_MSG_CHECKING([for nearpt3 support])
# AC_ARG_ENABLE([nearpt3], AS_HELP_STRING([--enable-nearpt3],[nearpt3 support [default=no]]))
# AS_IF([test "x$enable_nearpt3" = "xyes"], [
#   AC_DEFINE(ENABLE_NEARPT3, [1], [Define to 1 for nearpt3 support])
# ])
# AC_MSG_RESULT([$enable_nearpt3])
# AC_SUBST([ENABLE_NEARPT3],[$enable_nearpt3])
# AM_CONDITIONAL([ENABLE_NEARPT3],[test x$enable_nearpt3 = 'xyes'])
#  ----------------------------------------------------------------------
#  Enable DATA support
AC_MSG_CHECKING([for DATA support])
AC_ARG_ENABLE([data], AS_HELP_STRING([--enable-data],[DATA support [default=yes]]))
AS_IF([test "x$enable_data" != 'xno'], [
  AC_DEFINE(ENABLE_DATA, [1], [Define to 1 for DATA support])
])
AC_MSG_RESULT([$enable_data])
AC_SUBST([ENABLE_DATA],[$enable_data])
#  ----------------------------------------------------------------------
CFLAGS="$CFLAGS ${OPENMP_CFLAGS}"
CXXFLAGS="$CXXFLAGS ${OPENMP_CFLAGS}"
#  ----------------------------------------------------------------------

# Add configure options
ACX_OPTIONS

#  ----------------------------------------------------------------------

# configure code from valgrind

# does the x86/amd64 assembler understand SSE 4.2 instructions?
# Note, this doesn't generate a C-level symbol.  It generates a
# automake-level symbol (BUILD_SSE42_TESTS), used in test Makefile.am's
AC_MSG_CHECKING([if x86/amd64 assembler speaks SSE4.2])

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
  do { long long int x; 
   __asm__ __volatile__(
      "crc32q %%r15,%%r15" : : : "r15" );
   __asm__ __volatile__(
      "pblendvb (%%rcx), %%xmm11" : : : "memory", "xmm11"); 
   __asm__ __volatile__(
      "aesdec %%xmm2, %%xmm1" : : : "xmm2", "xmm1"); }
  while (0)
]])], [
ac_have_as_sse42=yes
AC_MSG_RESULT([yes])
], [
ac_have_as_sse42=no
AC_MSG_RESULT([no])
])

AM_CONDITIONAL(BUILD_SSE42_TESTS, test x$ac_have_as_sse42 = xyes)


# does the x86/amd64 assembler understand AVX instructions?
# Note, this doesn't generate a C-level symbol.  It generates a
# automake-level symbol (BUILD_AVX_TESTS), used in test Makefile.am's
AC_MSG_CHECKING([if x86/amd64 assembler speaks AVX])

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
  do { long long int x; 
   __asm__ __volatile__(
      "vmovupd (%%rsp), %%ymm7" : : : "xmm7" );
   __asm__ __volatile__(
      "vaddpd %%ymm6,%%ymm7,%%ymm8" : : : "xmm6","xmm7","xmm8"); }
  while (0)
]])], [
ac_have_as_avx=yes
AC_MSG_RESULT([yes])
], [
ac_have_as_avx=no
AC_MSG_RESULT([no])
])

AM_CONDITIONAL(BUILD_AVX_TESTS, test x$ac_have_as_avx = xyes)


# does the x86/amd64 assembler understand AVX2 instructions?
# Note, this doesn't generate a C-level symbol.  It generates a
# automake-level symbol (BUILD_AVX2_TESTS), used in test Makefile.am's
AC_MSG_CHECKING([if x86/amd64 assembler speaks AVX2])

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
  do { long long int x; 
   __asm__ __volatile__(
      "vpsravd (%%rsp), %%ymm8, %%ymm7" : : : "xmm7", "xmm8" );
   __asm__ __volatile__(
      "vpaddb %%ymm6,%%ymm7,%%ymm8" : : : "xmm6","xmm7","xmm8"); }
  while (0)
]])], [
ac_have_as_avx2=yes
AC_MSG_RESULT([yes])
], [
ac_have_as_avx2=no
AC_MSG_RESULT([no])
])

AM_CONDITIONAL(BUILD_AVX2_TESTS, test x$ac_have_as_avx2 = xyes)

#  ----------------------------------------------------------------------
AC_SUBST([CXXFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([FCFLAGS])
AC_SUBST([CLDFLAGS])
AC_SUBST([CLIBS])
#AC_SUBST([INCLUDES])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AR])

# Checks for compiler
C_COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(C_COMPILER, ["$C_COMPILER"], [C Compiler])
CXX_COMPILER="$CXX $CXXFLAGS"
AC_DEFINE_UNQUOTED(CXX_COMPILER, ["$CXX_COMPILER"], [CXX Compiler])

AC_REQUIRE_AUX_FILE([tap-driver.sh])
AC_PROG_AWK

AC_CONFIG_FILES([test/File.test test/Read_grib.test test/Read_netcdf.test test/Copy_netcdf.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Comp.test test/Compc.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Cat.test test/Gridarea.test test/Genweights.test test/Remap.test test/Remap2.test test/Remapeta.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/EOF.test test/Select.test test/Spectral.test test/Vertint.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Timstat.test test/Timselstat.test test/Seasstat.test test/Runstat.test test/Multiyearstat.test test/Ydrunstat.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Gridboxstat.test test/Vertstat.test test/Fldstat.test test/Fldpctl.test test/Ensstat.test test/Enspctl.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Merstat.test test/Zonstat.test test/Mergetime.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Afterburner.test test/Detrend.test test/Arithc.test test/Arith.test test/Expr.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/Gradsdes.test test/Collgrid.test test/threads.test test/tsformat.test test/wildcard.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([test/MapReduce.test test/Ninfo.test],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([Makefile src/Makefile contrib/Makefile test/Makefile test/data/Makefile cdo.spec cdo.settings])
AC_OUTPUT

#  ----------------------------------------------------------------------
#  Show configuration
AC_MSG_NOTICE([CDO is configured with the following options:])
cat cdo.settings

AS_IF([test "$ac_cv_prog_cc_c99" = "no"],
      [AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])])
AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDO package and 'make install' to install it afterwards.
])
