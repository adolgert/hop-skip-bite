#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([hopskip], [0.1], [ajd27@cornell.edu])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([src/rcpp_hello_world.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([src/Makevars])

# Checks for programs.
R_CXX_CMD=`R CMD CONFIG CXX1X`
AC_PROG_CXX([$R_CXX_CMD])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AX_PTHREAD
AX_BOOST_BASE([1.55])
AX_BOOST_CHRONO
AX_BOOST_FILESYSTEM
AX_BOOST_LOG
AX_BOOST_LOG_SETUP
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_INT64_T
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_CHECK_FUNCS([pow])



AC_OUTPUT
