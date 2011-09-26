PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=genomeSIMLA

HDRS=config.h \
	genomesim.h \
	timestamp.h \
	timestamp.l

#DEFINE THE sources
SRCS=config.cpp  \
	genomesim.cpp 

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
LOCAL_LIBS=simulation simla simpen genetics utility ga pngwriter random
#appinterface.l

#Uncomment this line for executables
#LIBRARY_MAKE=0

#Uncomment this line for test applications
IS_TEST=0

#appinterface.l: 
	#inc_version.sh timestamp.h

MAKEPATH=$(ROOTPATH)
FREETYPE=1
EXTERNAL_LINKS=-lz -lpthread -lpng
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base



OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



