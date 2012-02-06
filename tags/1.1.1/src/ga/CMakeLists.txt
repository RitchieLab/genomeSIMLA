###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 

set (GaFiles 
	garandom.C 
	gaerror.C 
	GAParameter.C 
	GAStatistics.C 
	GABaseGA.C 
	GASStateGA.C 
	GASimpleGA.C 
	GAIncGA.C 
	GADemeGA.C 
	GADCrowdingGA.C
	GASelector.C 
	GAScaling.C
	GAPopulation.C
	GAGenome.C
	GABinStr.C
	gabincvt.C
	GAAllele.C
	GAStringGenome.C
	GARealGenome.C
	GA1DBinStrGenome.C
	GA2DBinStrGenome.C
	GA3DBinStrGenome.C
	GABin2DecGenome.C
	GA1DArrayGenome.C
	GA2DArrayGenome.C
	GA3DArrayGenome.C
	GATreeBASE.C
	GATree.C
	GATreeGenome.C
	GAListBASE.C
	GAList.C
	GAListGenome.C
)


ADD_LIBRARY (ga ${GaFiles})

SET_SOURCE_FILES_PROPERTIES(${GaFiles}
	PROPERTIES
	LANGUAGE CXX)

