PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=paris

HDRS=analyzer.h\
	appconfiguration.h\
	bin.h\
	chromosome.h\
	feature.h\
	gene.h\
	genomicregion.h\
	knowledgebase.h\
	main.h\
	pathway.h\
	parisapp.h\
	timestamp.l

#ADLDEPS=
#timestamp.l

#DEFINE THE sources
SRCS=analyzer.cpp\
	appconfiguration.cpp\
	bin.cpp\
	chromosome.cpp\
	feature.cpp\
	gene.cpp\
	genomicregion.cpp\
	knowledgebase.cpp\
	main.cpp\
	pathway.cpp\
	parisapp.cpp

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
LOCAL_LIBS=utility random
#appinterface.l

#Uncomment this line for executables
#LIBRARY_MAKE=0

#Uncomment this line for test applications
IS_TEST=0
SOCI=1
MAKEPATH=$(ROOTPATH)
EXTERNAL_LINKS=
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



