#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(genomeSIMLA, 1.1.0, bugs@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/genomesim/genomesim.cpp])
AC_CONFIG_HEADER([config.h])

m4_include([make/configure.ac.base])

AC_CONFIG_FILES([Makefile
		 src/ga/Makefile
                 src/genetics/Makefile
                 src/genomesim/Makefile
                 src/pngwriter/Makefile
                 src/random/Makefile
                 src/simla/Makefile
                 src/simpen/Makefile
                 src/simulation/Makefile
                 src/utility/Makefile
                 ])
AC_OUTPUT
