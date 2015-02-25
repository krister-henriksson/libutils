#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libutils], [1], [krister.henriksson@helsinki.fi])
AC_CONFIG_SRCDIR([src/funcfit-basics.cpp])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC


AM_PROG_AR
LT_INIT



# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])
# FIXME: Replace `main' with a function in `-lrt':
AC_CHECK_LIB([rt], [clock_gettime])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
