noinst_LTLIBRARIES=libsimulation.la

libsimulation_la_LIBADD=\
	../utility/libutility.la \
	$(GSL_LIBS)

libsimulation_la_SOURCES=basicsample.cpp \
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
	templatedpedigree.cpp \
	allelesource.h \
	basicsample.h \
	blocklistnode.h \
	blocklistnodefourgammetes.h \
	chromosome.h \
	chrompool.h \
	cpair.h \
	cpool.h \
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
	templatedpedigree.h \
	ContinuousMultivarSample.h \
	ContinuousMultivarSample.cpp 
	
AM_CPPFLAGS=-I$(top_srcdir)/src $(FREETYPE_CPPFLAGS) $(GSL_CFLAGS)

AM_LDFLAGS=$(GSL_LDFLAGS)