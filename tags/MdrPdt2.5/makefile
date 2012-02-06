BINPATH=./bin/apps
ROOTPATH=.

LIBPATH=$(ROOTPATH)/bin/lib

ifeq ($(DEBUG), 1)
	LIBSUFFIX=d
	RELEASEMODE="RELEASE=0"
else
	LIBSUFFIX=
	RELEASEMODE="RELEASE=1"
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


ifdef COMPAT
	CPP=$(GNUCC)$(COMPAT)
else
	CPP=$(GNUCC)
endif



#OK, gcc version 3 seems to uses drand for stl random functions (random_shuffle), while
#version 4 uses drand48....nice! Let's set the right macro so we can set the right seed
GCC_VERSION=$(shell $(CXX) --version | grep gcc | cut -f3 -d\ | cut -f1 -d. )
ifeq ($(GCC_VERSION), 3)
	#This is required for certain distributions of the stdlib it seems- but not all
	FLAGS_RAND=-D MDRPDT_STL_USES_RAND48
endif


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



PROJECT_OPTIONS=$(OPTION_DOPT) $(OPTION_CV) $(DEBUG=1)

LOCAL_LIBS=utility genetics simulation

DEPENDANTLIBS=$(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))

all: $(BINPATH)/mdr-pdts

.PHONY: clean test serial parallel genomesim



simpen:
	@echo "Building SimPen\n";
	@$(MAKE) -C $(ROOTPATH)/src/simpen $(PROJECT_OPTIONS)

serial: 
	@echo "Building Serial version of MDR-PDT";
	@$(MAKE) -C $(ROOTPATH) $(BINPATH)/mdr-pdts USE_MPI=0 $(PROJECT_OPTIONS)

parallel: 
	@echo "Building Parallel (standard MPI. Not for VAMPIRE)";
	@$(MAKE) -C $(ROOTPATH) $(BINPATH)/mdr-pdtp USE_MPI=1 $(PROJECT_OPTIONS)

vampire:
	@echo "Building Parallel for vampire cluster";
	@$(MAKE) -C $(ROOTPATH) $(BINPATH)/mdr-pdtp USE_MPI=1 $(PROJECT_OPTIONS) VAMPIRE=1
	
power:
	@echo "Building Basic power application";
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt -fmake.power  $(PROJECT_OPTIONS)

genomesim:
	@echo "Building genomesim";
#	@$(MAKE) -C $(ROOTPATH)/src/simpen $(RELEASE) $(PROJECT_OPTIONS)
#	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASE) $(PROJECT_OPTIONS)

	@echo Executing: $(MAKE) $(BINPATH)/genomesims  $(PROJECT_OPTIONS) 
	@$(MAKE) -C $(ROOTPATH) $(BINPATH)/genomesims $(PROJECT_OPTIONS) 

test: $(BINPATH)/test/mdrpdt.test
	@echo "Building Basic power application";
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(PROJECT_OPTIONS) 

release:
	make COMPAT=32 

debug: 
	make DEBUG=1 USE_DRAND48=1

clean:
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) cleantmps

cleantmps:
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) cleantmps

$(BINPATH)/genomesims:  $(DEPENDANTLIBS)
	echo "LIBPATH: $(LIBPATH)"	
	echo $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/genomesimp:  	$(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/src/genomesim $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/mdr-pdt-power: 	$(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt -f$(ROOTPATH)/src/make.power $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/mdr-pdts: $(DEPENDANTLIBS)
#$(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	echo "LIBPATH: $(LIBPATH)"	
	echo $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt $(RELEASEMODE) $(PROJECT_OPTIONS) USE_MPI=0 

$(BINPATH)/mdr-pdtp: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/src/mdrpdt $(RELEASEMODE) $(PROJECT_OPTIONS) USE_MPI=1

$(LIBPATH)/libgenetics$(LIBSUFFIX).a:
	echo "Building libgenetics"
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimulation$(LIBSUFFIX).a:
	echo "Building libsimulation"
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libutility$(LIBSUFFIX).a:
	echo "Building libutility"
	echo $(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) $(PROJECT_OPTIONS)
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/test/mdrpdt.test: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/test $(RELEASEMODE) $(PROJECT_OPTIONS)
