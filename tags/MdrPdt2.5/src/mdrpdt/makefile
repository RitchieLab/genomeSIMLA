#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
APPHOME=apps
ROOTPATH=../..

PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

#DEFINE THE sources
SRCS=eseconfiguration.cpp \
	evalbalancedaccuracy.cpp \
	evalbalancedaccuracypdt.cpp \
	evalmaxdifference.cpp \
	gtlineparsermdrpdt.cpp \
	locuslogascii.cpp \
	matchedoddsratio.cpp \
	modelfinder.cpp \
	pdtfold.cpp \
	reportlog.cpp \
	reportlogascii.cpp \
	reportlogbinary.cpp \
	sibshipfoldproduction.cpp \
	snpsearchapplication.cpp

#OVERRIDE the path to sources if necessary
SRCPATHS=.. .

#indicate what the executable will be named
PROJECTNAME=mdr-pdt

#list the projects which build into libraries
LOCAL_LIBS=genetics utility


MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base

