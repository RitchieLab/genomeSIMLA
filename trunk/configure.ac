#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(genomeSIMLA, 1.2.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

RL_CONFIGURE([1.40.0],[HAVE_PNG HAVE_BOOST HAVE_GSL])

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/genomesim/genomesim.cpp])
AC_CONFIG_HEADER([config.h])

AC_CONFIG_FILES([Makefile
                 src/genomesim/Makefile
                 src/ga/Makefile
                 src/genetics/Makefile
                 src/pngwriter/Makefile
                 src/random/Makefile
                 src/simla/Makefile
                 src/simpen/Makefile
                 src/simulation/Makefile
                 src/utility/Makefile
                 ])
AC_OUTPUT
