APPHOME=apps/testing
PROJECTNAME=unitese

#Define the sources
MAIN_SRCS=testpdttstat.cpp \
	testpdttstatreal.cpp \
	snprepostest.cpp 
	
MAIN_PATH=./mdr

GEN_TEST_PATH=./genetics
GEN_TEST_SRCS=snpalignedtest.cpp snppooltest.cpp fileparsertests.cpp
#peddatatest.cpp 

UTIL_TEST_PATH=./utility
UTIL_TEST_SRCS=binarrayparsertest.cpp lineparsing.cpp testcasecontrolstatus.cpp

ESE_PATH=$(ROOTPATH)/src
ESE_SRC=appconfiguration.cpp \
	evalbalancedaccuracypdt.cpp \
	gtlineparsermdrpdt.cpp \
	locuslogascii.cpp \
	matchedoddsratio.cpp \
	reportlog.cpp \
	reportlogascii.cpp \
	reportlogbinary.cpp 

SRCS=main.cpp $(patsubst %,$(MAIN_PATH)/%,$(MAIN_SRCS)) $(patsubst %,$(UTIL_TEST_PATH)/%,$(UTIL_TEST_SRCS)) $(patsubst %, $(GEN_TEST_PATH)/%, $(GEN_TEST_SRCS)) $(patsubst %, $(ESE_PATH)/%, $(ESE_SRC)) $(patsubst %, $(MDRPDT_TEST_PATH)/%, $(MDRPDT_TEST_SRCS))
#SRCS=main.cpp $(MAIN_SRCS) $(UTIL_TEST_SRCS) $(GEN_TEST_SRCS) $(ESE_SRC)
#VPATH=$(GEN_TEST_PATH):$(MAIN_PATH):$(UTIL_TEST_PATH):$(ESE_PATH)
SRCPATHS=./mdr ./genetics ./utility $(ROOTPATH)/src

PROJECT_COMPILER_FLAGS=
PROJECT_LINKER_FLAGS=

TARGET_DIR=..

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=..

#OVERRIDE the path to sources if necessary
#SRCPATH=.

ADDITIONAL_INC=$(ROOTPATH)/src

#list the projects which build into libraries
LOCAL_LIBS=genetics utility

TARGET_OBJS=$(patsubst %,$(TARGET_DIR)/$(OBJECT_PATH)/%, $(TEST_OBJS))

#Uncomment this line for executables
LIBRARY_MAKE=0

#Uncomment this line for test applications
TEST=1



MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base

#INCLUDES=-I.. -I$(BINPATH)/inc -I$(GEN_TEST_PATH) -I$(UTIL_TEST_PATH) $(patsubst %, -I$(ROOTPATH)/dist/inc/%, $(LOCAL_LIBS))


