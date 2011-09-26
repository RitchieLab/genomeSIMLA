#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
APPHOME=apps
ROOTPATH=../../..

PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

#DEFINE THE sources
SRCS=config.cpp \
	converter.cpp 

#OVERRIDE the path to sources if necessary
SRCPATHS=../.. . /usr/include/freetype2 $(3RDPARTY)/include/freetype2

#indicate what the executable will be named
PROJECTNAME=HapMapConverter

#list the projects which build into libraries
LOCAL_LIBS=utility

#EXTERNAL_LINKS=/home/torstenson/projects/MdrPdt/bin/lib/libga.a
#EXTERNAL_LINKS=`freetype-config --cflags` `freetype-config --libs` -lz -lpng -lpngwriter -lpthread


MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base

