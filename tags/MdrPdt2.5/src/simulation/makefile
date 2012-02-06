PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=simulation

#DEFINE THE sources
SRCS=basicsample.cpp \
	chromosome.cpp \
	chrompool.cpp \
	diseasemodel.cpp \
	distmappingfn.cpp \
	exponentialgrowth.cpp \
	gamodel.cpp \
	growthrate.cpp \
	locus.cpp \
	logisticgrowth.cpp \
	lineargrowth.cpp \
	individual.cpp \
	modelmanager.cpp \
	pedigreesample.cpp \
	penetrancemodel.cpp \
	poolmanager.cpp 
	
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




