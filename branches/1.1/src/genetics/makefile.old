PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=genetics

HDRS=allelicdecoder.h \
	binarygenoparser.h \
	eliminategenotypeerrors.h \
	familyevaluation.h \
	familynode.h \
	familymember.h \
	familyrepoevaluationbasic.h \
	familyrepoevaluation.h \
	familyrepository.h \
	fixtrios.h \
	genotypedata.h \
	gtevaluation.h \
	gtfileparser.h \
	gtfileparserbuffered.h \
	gtlineparser.h \
	modelstatistics.h \
	ptestdistribution.h \
	resultsrepository.h \
	snpaligned.h  \
	snpevalsuite.h \
	snppool.h \
	snpcontainer.h \
	snprecipient.h \
	snprepository.h \
	snprepostxtsorted.h \
	snpverificationmethod.h
 
#DEFINE THE sources
SRCS=allelicdecoder.cpp \
	binarygenoparser.cpp \
	eliminategenotypeerrors.cpp \
	familyevaluation.cpp \
	familynode.cpp \
	familymember.cpp \
	familyrepoevaluationbasic.cpp \
	familyrepoevaluation.cpp \
	familyrepository.cpp \
	fixtrios.cpp \
	genotypedata.cpp \
	gtevaluation.cpp \
	gtfileparser.cpp \
	gtfileparserbuffered.cpp \
	gtlineparser.cpp \
	modelstatistics.cpp \
	ptestdistribution.cpp \
	resultsrepository.cpp \
	snpaligned.cpp  \
	snpevalsuite.cpp \
	snppool.cpp \
	snpcontainer.cpp \
	snprecipient.cpp \
	snprepository.cpp \
	snprepostxtsorted.cpp \
	snpverificationmethod.cpp 


#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..

#indicate what the executable will be named
#EXENAME=ese

#list the projects which build into libraries
#LOCAL_LIBS=genetics utility

#Uncomment this line for executables
LIBRARY_MAKE=1

#Uncomment this line for test applications
IS_TEST=0



MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



