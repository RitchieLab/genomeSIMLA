#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
APPHOME=apps
ROOTPATH=../..

PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

#DEFINE THE sources
SRCS=config.cpp \
	genomesim.cpp 

#OVERRIDE the path to sources if necessary
SRCPATHS=.. .

#indicate what the executable will be named
PROJECTNAME=genomesim

#list the projects which build into libraries
LOCAL_LIBS=simulation simpen utility ga

#EXTERNAL_LINKS=/home/torstenson/projects/MdrPdt/bin/lib/libga.a


MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base

