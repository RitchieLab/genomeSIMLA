PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=genetics

#DEFINE THE sources
SRCS=binarygenoparser.cpp \
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


#list the projects which build into libraries
#LOCAL_LIBS=genetics utility

#Uncomment this line for executables
LIBRARY_MAKE=1

#Uncomment this line for test applications
#IS_TEST=0


MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base


#I need to figure out how to use wildcards without catching multiple files model and 
#$(OBJS): $(SRCS)
#	$(SUPRESS)mkdir -p $(OBJECT_PATH)
#	$(SUPRESS)$(CPP) -c $(CPPFLAGS) $(filter $(*F).cpp, $^) -o $(OBJECT_PATH)/$(*F).o




