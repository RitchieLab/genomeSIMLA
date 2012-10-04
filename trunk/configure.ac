#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(genomeSIMLA, 1.1.1, software@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/genomesim/genomesim.cpp])
AC_CONFIG_HEADER([config.h])

RL_CONFIGURE([1.40.0],[HAVE_PNG HAVE_BOOST HAVE_GSL])

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
