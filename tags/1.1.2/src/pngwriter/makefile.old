PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=pngwriter

HDRS=pngwriter.h
#DEFINE THE sources
SRCS=pngwriter.cpp

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..

#Uncomment this line for executables
LIBRARY_MAKE=1

#Uncomment this line for test applications
IS_TEST=0

FREETYPE=1

MAKEPATH=$(ROOTPATH)
#LIB_LINK:=$(shell freetype-config --cflags)
FREETYPE=1
include $(MAKEPATH)/extlibs.make

include $(MAKEPATH)/make.base

OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



