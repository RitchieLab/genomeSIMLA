PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=ga

#DEFINE THE sources
CSRCS=garandom.C \
	gaerror.C \
	GAParameter.C \
	GAStatistics.C \
	GABaseGA.C \
	GASStateGA.C \
	GASimpleGA.C \
	GAIncGA.C \
	GADemeGA.C \
	GADCrowdingGA.C \
	GASelector.C \
	GAScaling.C \
	GAPopulation.C \
	GAGenome.C \
	GABinStr.C \
	gabincvt.C \
	GAAllele.C \
	GAStringGenome.C \
	GARealGenome.C \
	GA1DBinStrGenome.C \
	GA2DBinStrGenome.C \
	GA3DBinStrGenome.C \
	GABin2DecGenome.C \
	GA1DArrayGenome.C \
	GA2DArrayGenome.C \
	GA3DArrayGenome.C \
	GATreeBASE.C \
	GATree.C \
	GATreeGenome.C \
	GAListBASE.C \
	GAList.C \
	GAListGenome.C 

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



MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



