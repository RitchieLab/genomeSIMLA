###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
IF (USE_XY)
	set (XY-Files 
		allelesource.cpp allelesource.h
		cpair.cpp cpair.h
		cpool.cpp cpool.h
		cpoolxy.cpp cpoolxy.h
		ldcalculator.cpp ldcalculator.h
		ldhtmlreport.cpp ldhtmlreport.h
		ldplotmanager.cpp ldplotmanager.h
		locusxy.cpp locusxy.h
		locusmanager.cpp locusmanager.h
		par_region.cpp par_region.h
	)
ELSE (USE_XY)
	SET (XY-Files )
ENDIF (USE_XY)

set (SimulationFiles 
	basicsample.cpp basicsample.h
	blocklistnode.cpp blocklistnode.h
	blocklistnodefourgammetes.cpp blocklistnodefourgammetes.h
	chromosome.cpp chromosome.h
	chrompool.cpp chrompool.h
	diseasemodel.cpp diseasemodel.h
	distmappingfn.cpp distmappingfn.h
	exponentialgrowth.cpp exponentialgrowth.h
	gamodel.cpp gamodel.h
	growthrate.cpp growthrate.h
	kinshipcalculator.cpp kinshipcalculator.h
	ldfooterregion.cpp ldfooterregion.h
	ldprimepngwriter.cpp ldprimepngwriter.h
	ldlabelregion.cpp ldlabelregion.h
	ldlineplot.cpp ldlineplot.h
	ldlocationbar.cpp ldlocationbar.h
	ldmafregion.cpp ldmafregion.h
	ldpngcomponent.cpp ldpngcomponent.h
	ldplotter.cpp ldplotter.h
	ldwriter.cpp ldwriter.h
	locus.cpp locus.h
	locusassociation.cpp locusassociation.h
	locreport.cpp locreport.h
	locselection.cpp locselection.h
	logisticgrowth.cpp logisticgrowth.h
	lineargrowth.cpp lineargrowth.h
	individual.cpp individual.h
	meantable.cpp meantable.h
	minchrom.cpp minchrom.h
	pedigreereferencesample.cpp pedigreereferencesample.h
	pedigreesample.cpp pedigreesample.h
	penetranceeval.cpp penetranceeval.h
	penetrancemodel.cpp penetrancemodel.h
	penfilemodel.cpp penfilemodel.h
	poolmanager.cpp poolmanager.h
	richardslogistic.cpp richardslogistic.h
	simlamodel.cpp simlamodel.h
)
#par_region.cpp par_region.h

set (simulationTestFiles
	main.cpp
)
IF (USE_XY)
	SET_SOURCE_FILES_PROPERTIES ( ${SimulationFiles} ${XY-Files}
			PROPERTIES COMPILE_FLAGS "-DUSE_XY=1" )
ENDIF (USE_XY)

ADD_LIBRARY (simulation ${SimulationFiles} ${XY-Files} )
SET_TARGET_PROPERTIES( simulation PROPERTIES DEBUG_POSTFIX d)

IF (MAKE_TESTS)
ADD_EXECUTABLE (unit ${SimulationFiles} ${XY-Files} ${simulationTestFiles})
SET_TARGET_PROPERTIES( unit PROPERTIES DEBUG_POSTFIX d)
TARGET_LINK_LIBRARIES (unit simla simpen ga genetics utility random cppunit dl pngwriter ${LIBPNG} z pthread ${FREETYPE_LIBRARIES} )
ENDIF (MAKE_TESTS)
