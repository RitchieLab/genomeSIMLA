#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
APPHOME=apps
ROOTPATH=..

ifdef USE_DOPT
	PROJECT_COMPILER_FLAGS=-D USE_DOPTIMIZATION
endif

#DEFINE THE sources
SRCS=appconfiguration.cpp \
	evalbalancedaccuracypdt.cpp \
	gtlineparsermdrpdt.cpp \
	locuslogascii.cpp \
	matchedoddsratio.cpp \
	modelfinder.cpp \
	reportlog.cpp \
	reportlogascii.cpp \
	reportlogbinary.cpp 

#OVERRIDE the path to sources if necessary
#SRCPATHS=.

#indicate what the executable will be named
PROJECTNAME=mdr-pdt

#list the projects which build into libraries
LOCAL_LIBS=genetics utility


MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base

