PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=biofilter
DEADHRD=snptranslator.h
HDRS=appconfiguration.h \
	bioapplication.h \
	biofilter.h \
	genegenemodel.h \
	genegenemodelarchive.h \
	ldcorrection.h \
	kbentity.h \
	kbgroup.h \
	kbmetagroup.h \
	kbregion.h \
	modelreport.h \
	regionspline.h \
	snpmanager.h \
	ldsplineimporter.h \
	timestamp.l

#ADLDEPS=
#timestamp.l

#DEFINE THE sources
SRCS=appconfiguration.cpp \
	bioapplication.cpp \
	biofilter.cpp \
	genegenemodel.cpp \
	genegenemodelarchive.cpp \
	ldcorrection.cpp \
	kbentity.cpp \
	kbgroup.cpp \
	kbmetagroup.cpp \
	kbregion.cpp \
	modelreport.cpp \
	regionspline.cpp \
	snpmanager.cpp \
	ldsplineimporter.cpp 

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
LOCAL_LIBS=utility ldspline
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
	



