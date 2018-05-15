# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([network], [1.0])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign parallel-tests subdir-objects])
AC_CONFIG_SRCDIR([src/network.cc])
AC_CONFIG_HEADERS([config.h])
m4_pattern_allow([AC_OPENMP])


# Enable use of Google Sparsehash and get include location if requested.
AC_ARG_ENABLE([sparsehash],
    AS_HELP_STRING([--enable-sparsehash], [Enable the use of Google Sparsehash (must be installed).]))

AS_IF([test "x$enable_sparsehash" = "xyes"], [
AC_DEFINE([HAVE_SPARSEHASH], [1], [Enable compilation with Google Sparsehash])
AC_ARG_WITH([sparsehash-include-path],
  [AS_HELP_STRING([--with-sparsehash-include-path],
    [location of Google Sparsehash headers])],
  [SPARSEHASH_CFLAGS="-I$withval"],
  [SPARSEHASH_CFLAGS=""])
AC_SUBST([SPARSEHASH_CFLAGS])

])

# Checks for programs.
AC_PROG_CXX
AM_PROG_CC_C_O
AX_CXX_COMPILE_STDCXX_11

# Checks for libraries.
AC_CHECK_LIB([gomp], [omp_set_num_threads],, [AC_MSG_WARN([Unable to locate OpenMP library, you should probably specify '--disable-openmp'.])])
AC_CHECK_LIB([z], [gzread],, [AC_MSG_ERROR([Zlib not found, reading gzipped files will not be possible.])]) 

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])

# Check for OpenMP parallel execution support
AC_OPENMP

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt])

# For test harness
AC_PROG_AWK

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
