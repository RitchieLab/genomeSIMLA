BINPATH=./bin/apps/
ROOTPATH=.

ifeq ($(DEBUG), 1)
	LIBSUFFIX=d
	RELEASE="RELEASE=0"
else
	LIBSUFFIX=
	RELEASE="RELEASE=1"
endif

ifeq ($(USE_DOPT), 1)
	PROJECT_OPTIONS=USE_DOPT=1
endif

LOCAL_LIBS=utility genetics



all: $(BINPATH)/ese

.PHONY: clean test 

test: $(BINPATH)/test/mdrpdt.test

debug: 
	make DEBUG=1

clean:
	@$(MAKE) -C $(ROOTPATH)/src $(RELEASE) clean
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASE) clean
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASE) clean

cleantmps:
	@$(MAKE) -C $(ROOTPATH)/src $(RELEASE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASE) cleantmps
	

$(BINPATH)/ese: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/src $(RELEASE) $(PROJECT_OPTIONS)

$(LIBPATH)/libgenetics$(LIBSUFFIX).a:
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASE) $(PROJECT_OPTIONS)

$(LIBPATH)/libutility$(LIBSUFFIX).a:
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASE) $(PROJECT_OPTIONS)

$(BINPATH)/test/mdrpdt.test: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/test $(RELEASE) $(PROJECT_OPTIONS)