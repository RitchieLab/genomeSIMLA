PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=simulation
HDRS=allelesource.h \
	basicsample.h \
	blocklistnode.h \
	blocklistnodefourgammetes.h \
	chromosome.h \
	chrompool.h \
	cpair.h \
	cpool.h \
	cpoolxy.h \
	diseaselocus.h \
	diseasemodel.h \
	distmappingfn.h \
	exponentialgrowth.h \
	gamodel.h \
	growthrate.h \
	haplotypeblock.h \
	individual.h \
	kinshipcalculator.h \
	ldcalculator.h \
	ldfooterregion.h \
	ldlabelregion.h \
	ldlineplot.h \
	ldlocationbar.h \
	ldmafregion.h \
	ldplotmanager.h \
	ldplotter.h \
	ldpngcomponent.h \
	ldprimepngwriter.h \
	ldwriter.h \
	lineargrowth.h \
	locreport.h \
	locselection.h \
	locus.h \
	locusassociation.h \
	locusmanager.h \
	locusxy.h \
	logisticgrowth.h \
	meantable.h \
	minchrom.h \
	modelmanager.h \
	par_region.h \
	pedigreereferencesample.h \
	pedigreesample.h \
	penetranceeval.h \
	penetrancemodel.h \
	penfilemodel.h \
	poolmanager.h\
	richardslogistic.h \
	simlamodel.h \
	simulation.h \
	templatedpedigree.h

#DEFINE THE sources
SRCS=basicsample.cpp \
	blocklistnode.cpp \
	blocklistnodefourgammetes.cpp \
	chromosome.cpp \
	chrompool.cpp \
	diseasemodel.cpp \
	distmappingfn.cpp \
	exponentialgrowth.cpp \
	gamodel.cpp \
	growthrate.cpp \
	kinshipcalculator.cpp \
	ldfooterregion.cpp \
	ldprimepngwriter.cpp \
	ldlabelregion.cpp \
	ldlineplot.cpp \
	ldlocationbar.cpp \
	ldmafregion.cpp \
	ldpngcomponent.cpp \
	ldplotter.cpp \
	ldwriter.cpp \
	locus.cpp \
	locusassociation.cpp \
	locreport.cpp \
	locselection.cpp \
	logisticgrowth.cpp \
	lineargrowth.cpp \
	individual.cpp \
	meantable.cpp \
	minchrom.cpp \
	pedigreereferencesample.cpp \
	pedigreesample.cpp \
	penetranceeval.cpp \
	penetrancemodel.cpp \
	penfilemodel.cpp \
	poolmanager.cpp \
	richardslogistic.cpp \
	simlamodel.cpp \
	templatedpedigree.cpp

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

#FT_CPPFLAGS := $(shell freetype-config --cflags)

#FT_LINK := $(shell freetypetype-config --libs)
#EXTERNAL_LINKS=-l png $(FT_LINK) -lz -lpthread 
#ALT_INCLUDES=$(FT_CPPFLAGS)

#EXTERNAL_LINKS=-l png -L /usr/X11R6/lib/libfreetype.dylib -lz -lpthread
FREETYPE=1

MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



