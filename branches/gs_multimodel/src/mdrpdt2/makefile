PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=mdrpdt

HDRS=appconfig.h\
	appmdrpdt.h\
	evaluationmethod.h\
	familyrepository.h\
	foldproduction.h\
	foldstatistic.h\
	genotypecounter.h\
	genotyperepository.h\
	individual.h\
	matchedoddsratio.h\
	mdrpdt.h\
	missingdataevaluation.h\
	orderreport.h\
	pdtmodel.h\
	pedigree.h\
	pedigreecleaner.h\
	pedigreeexclusion.h\
	pedigreerepositorytraverser.h\
	pedigreestatistics.h\
	ptestdistribution.h\
	sibship.h\
	tfinalreport.h\
	timestamp.h\
	tstatistic.h\
	timestamp.l

#DEFINE THE sources
SRCS=appconfig.cpp\
	appmdrpdt.cpp\
	familyrepository.cpp\
	foldproduction.cpp\
	foldstatistic.cpp\
	genotyperepository.cpp\
	individual.cpp\
	matchedoddsratio.cpp\
	mdrpdt.cpp\
	missingdataevaluation.cpp\
	orderreport.cpp\
	pdtmodel.cpp\
	pedigree.cpp\
	pedigreecleaner.cpp\
	pedigreeexclusion.cpp\
	pedigreerepositorytraverser.cpp\
	pedigreestatistics.cpp\
	ptestdistribution.cpp\
	sibship.cpp\
	tfinalreport.cpp\
	tstatistic.cpp
	
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

#appinterface.l: 
	#inc_version.sh timestamp.h

MAKEPATH=$(ROOTPATH)

EXTERNAL_LINKS=-lpthread
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base



OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



