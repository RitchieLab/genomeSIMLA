BINPATH=./bin
ROOTPATH=.

include $(ROOTPATH)/platform.make
include $(ROOTPATH)/color.make

ifeq ($(INC_VERSION), 1)
	INCREMENT=INC_BUILD=1
endif 

ifeq ($(DEBUG), 1)
	LIBSUFFIX=d
	RELEASEMODE=RELEASE=0
	RELEASE_DESC="Debug Mode"
	INCREMENT=
else
	LIBSUFFIX=
	RELEASEMODE=RELEASE=1
	RELEASE_DESC="Release Mode"
endif

ifeq ($(USE_DOPT), 1)
	OPTION_DOPT=USE_DOPT=1
endif

ifeq ($(ENABLE_CV), 0)
	OPTION_CV=
	CV_SUFF=
else
	OPTION_CV=ENABLE_CV=1
	CV_SUFF="_cv"
endif
GNUCC=g++

ifdef COMPAT
	CPP=$(GNUCC)$(COMPAT)
else
	CPP=$(GNUCC)
endif




ifeq ($(USE_MPI), 1)
	FLAGS_MPI=-D USE_MPI
	MPI_LINK=$(MPI_LIBS)
	EXEC_PREFIX=p
	MPIPATH=parallel
	MPI_DESC="w/ MPI"
else
	MPI_LINK=
	EXEC_PREFIX=
	MPIPATH=serial
	MPI_DESC="w/o MPI"
endif

include $(ROOTPATH)/platform.make
#LIBPATH=$(ROOTPATH)/bin/lib

LIBPATH=$(ROOTPATH)/bin/lib/$(MPIPATH)/$(PLATFORM)-$(BITCOUNT)bit/$(CPP)

GCC_VERSION=$(shell $(CXX) --version | grep gcc | cut -f3 -d\ | cut -f1-2 -d. )


ifdef VAMPIRE
    MPI_LIBS= -lmpich -lpthread
    MPI_INC=/usr/local/mpich-ether/p4/absoft_gcc/include/
    CL=mpiCC
	VMP_SUFF="_vmp"
else
	MPI_LIBS=-lmpichcxx -lmpich -lpthread
    MPI_INC=
    CL=$(CPP)
	VMP_SUFF=
endif



PROJECT_OPTIONS=-s $(OPTION_DOPT) $(OPTION_CV) $(RELEASE_MODE) $(INCREMENT)

LOCAL_LIBS=ldspline pngwriter utility genetics simulation simla simpen ga random


DEPENDANTLIBS=$(patsubst %, $(LIBPATH)/lib%$(EXEC_PREFIX)$(LIBSUFFIX).a, $(LOCAL_LIBS))

all: $(BINPATH)/genomeSIMLA$(LIBSUFFIX)

.PHONY: clean test serial parallel genomesim


mdrpdt: $(BINPATH)/mdr-pdt$(LIBSUFFIX)$(EXEC_PREFIX)

bc:
	@echo "Building Binary Compression"
	@$(MAKE) -C $(ROOTPATH)/src/tools -fmake.binarycompress $(RELEASEMODE) $(PROJECT_OPTIONS)
	


serial: $(BINPATH)/mdr-pdt
parallel: $(BINPATH)/mdr-pdtp
vampire:
	@echo "Building Parallel for vampire cluster";
	@$(MAKE) -C $(ROOTPATH) $(BINPATH)/mdr-pdtp USE_MPI=1 $(PROJECT_OPTIONS) VAMPIRE=1
	
power: 
	@echo "Building Basic power application";
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt -fmake.power  $(PROJECT_OPTIONS)

powerp: $(BINPATH)/mdr-pdtp
	@echo "Building Basic parallel power application";
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt -fmake.power USE_MPI=1 $(PROJECT_OPTIONS)

power-vampire: $(BINPATH)/mdr-pdtp
	@echo "Building Basic power application for vampire";
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt -fmake.power USE_MPI=1 $(PROJECT_OPTIONS) VAMPIRE=1

gs: $(BINPATH)/gs$(LIBSUFFIX)
gs-gui: $(BINPATH)/wxGenomeSIMLA$(LIBSUFFIX)
simla: $(LIBPATH)/libsimla$(LIBSUFFIX).a
ga: $(LIBPATH)/libga$(LIBSUFFIX).a
simpen: $(LIBPATH)/libsimpen$(LIBSUFFIX).a
simulation: $(LIBPATH)/libsimulation$(LIBSUFFIX).a
utility: $(LIBPATH)/libutility$(LIBSUFFIX).a
random: $(LIBPATH)/librandom$(LIBSUFFIX).a
genomesim: $(LIBPATH)/libgenomesim$(LIBSUFFIX).a
pngwriter: $(LIBPATH)/libpngwriter$(LIBSUFFIX).a
ldspline: $(BINPATH)/ldspline$(LIBSUFFIX)
biofilter: $(BINPATH)/biofilter$(LIBSUFFIX)
paris: $(BINPATH)/paris$(LIBSUFFIX)
paris-test: $(BINPATH)/paris-test$(LIBSUFFIX)


test: $(BINPATH)/test/mdrpdt.test
	@echo "Building Basic power application";
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(PROJECT_OPTIONS) 

release:
	make COMPAT=32 

debug: 
	make DEBUG=1 

clean:
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) clean

cleantmps:
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) cleantmps
$(BINPATH)/gs$(LIBSUFFIX): $(DEPENDANTLIBS)
	$(MAKE) -C $(ROOTPATH)/src/gs $(RELEASEMODE) $(PROJECT_OPTIONS) 
$(BINPATH)/paris-test$(LIBSUFFIX): $(DEPENDANTLIBS)
	$(MAKE) -C $(ROOTPATH)/src/paris -f make.test $(RELEASEMODE) $(PROJECT_OPTIONS) INC_BUILD=0
$(BINPATH)/paris$(LIBSUFFIX): $(DEPENDANTLIBS)
	$(MAKE) -C $(ROOTPATH)/src/paris $(RELEASEMODE) $(PROJECT_OPTIONS)
$(BINPATH)/biofilter$(LIBSUFFIX): $(DEPENDANTLIBS)
	$(MAKE) -C $(ROOTPATH)/src/biofilter $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/wxGenomeSIMLA$(LIBSUFFIX): $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/wxgenomesimla $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/genomeSIMLA$(LIBSUFFIX):  $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/genomeSIMLAp:  	$(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/mdr-pdt-power: 	$(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt2 -f$(ROOTPATH)/src/make.power $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/mdr-pdt$(LIBSUFFIX): $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt2 $(RELEASEMODE) $(PROJECT_OPTIONS) USE_MPI=0 


$(BINPATH)/mdr-pdt$(LIBSUFFIX)p: $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt2 $(RELEASEMODE) $(PROJECT_OPTIONS) USE_MPI=1

$(BINPATH)/ldspline$(LIBSUFFIX): $(DEPENDENTLIBS) $(LIBPATH)/libldspline$(EXEC_PREFIX)$(LIBSUFFIX).a
	@$(MAKE) -C $(ROOTPATH)/src/ldspline -f main.make $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libldspline$(EXEC_PREFIX)$(LIBSUFFIX).a:
	@$(MAKE) -C $(ROOTPATH)/src/ldspline $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libgenetics$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/librandom$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/random $(RELEASEMODE) $(PROJECT_OPTIONS)
	
$(LIBPATH)/libpngwriter$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/pngwriter $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimpen$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simpen $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimla$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simla $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libga$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/ga $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimulation$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libutility$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) $(PROJECT_OPTIONS)


$(BINPATH)/test/mdrpdt.test: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/test $(RELEASEMODE) $(PROJECT_OPTIONS)
