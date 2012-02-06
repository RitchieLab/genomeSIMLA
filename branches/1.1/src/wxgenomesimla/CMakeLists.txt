###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
MESSAGE ("-------------------- GenomeSIMLA GUI--------------------")
SET (TARGET_APPLICATION wxGenomeSIMLA)
PROJECT(wxGenomeSIMLA)
SET (wxGenomeSIMLAFiles
	appconfig.cpp 
	appcontroller.cpp 
	appinterface.cpp 
	diseasemodeldetails.cpp 
	dlggrowthreview.cpp 
	genomesimla.cpp 
	growthrateconfig.cpp 
	locusmanager.cpp 
	mainframe.cpp 
	modellist.cpp 
	populationcontroller.cpp 
	treelistctrl.cpp
	wxbeginsimulation.cpp 
	wxchromcfgdialog.cpp 
	wxdialogtaskalert.cpp
	wxdlgconfigureccdatasets.cpp 
	wxdlgconfigurediseasemodel.cpp 
	wxdlgconfigurepeddatasets.cpp 
	wxdlgeditlocusselector.cpp 
	wxdlgbasichtmlreport.cpp 
	wxdlgprocessmonitor.cpp 
	wxdlgsavecurrent.cpp 
	wxdlgselectlocusrange.cpp 
	wxdlgsimla.cpp 
	wxgridloci.cpp 
	wximggrowthchart.cpp 
	wxldanalysis.cpp 
	wxlocuslistcontrol.cpp 
	wxpageinfo.cpp 
	wxpagebbchrom.cpp 
	wxpagegeneralsettings.cpp 
	wxpageloadchromosome.cpp 
	wxpagelocicreation.cpp 
	wxpagelocusreporting.cpp 
	wxpagepenetrancemodel.cpp 
	wxpagereporting.cpp 
	wxpageselectgeneration.cpp 
	wxpagesimpen.cpp 
	wxpanelchromosomerep.cpp 
	wxpanelreports.cpp 
	wxuser.cpp 
	wxwizardsimulatedata.cpp 
	wxwizardcreatechromosome.cpp 
	wxwizpagebbchrom.cpp 
	wxwizpagedefinedatasets.cpp 
	wxwizpagelaunchconfirmation.cpp 
	wxwizpagelaunchdroppoints.cpp 
	wxwizpagelaunchproject.cpp 
	wxwizpagelaunchselectgeneration.cpp 
	wxwizpagelaunchsimulation.cpp 
	wxwizpagereview.cpp 
	wxwizpagereviewsimulation.cpp 
	wxwizpageselectgeneration.cpp 
	wxwizpageselectloci.cpp 
	wxwizrunsimulation.cpp  
)
#	treelistctrl.cpp 

#MAC OS X Bundle
IF (APPLE) 
	MESSAGE ("Setting up details for creating MAC Bundle")
	SET (MACOSX_BUNDLE_BUNDLE_NAME wxGenomeSIMLA)
	SET (MACOSX_BUNDLE_BUNDLE_GUI_INDENTIFIER "edu.vanderbilt.wxgenomesimla")
	SET (PLATFORM_EXEC_TYPE MACOSX_BUNDLE)
ENDIF(APPLE)

SET (wxWidgets_USE_UNICODE, FALSE)
SET (wxWidgets_USE_LIBS aui richtext adv html core xml base xrc qa net)
SET (wxWidgets_FIND_COMPONENTS true)
FIND_PACKAGE(wxWidgets REQUIRED)
INCLUDE (${wxWidgets_USE_FILE})

IF (ZINCREMENT_BUILD) 
	MESSAGE ("ADDING increment build command")
	ADD_CUSTOM_COMMAND ( OUTPUT appinterface.l
		COMMAND build.sh ARGS ${PROJECT_SOURCE_DIR}/appinterface.cpp
		COMMENT "Incrementing build"
	)
	SET (wxGenomeSIMLAFiles ${wxGenomeSIMLAFiles} appinterface.l)
ENDIF (ZINCREMENT_BUILD)

ADD_EXECUTABLE (wxGenomeSIMLA ${PLATFORM_EXEC_TYPE} ${wxGenomeSIMLAFiles})
SET_TARGET_PROPERTIES (wxGenomeSIMLA PROPERTIES DEBUG_POSTFIX d)
SET_SOURCE_FILES_PROPERTIES (appinterface.l PROPERTIES SYMBOLIC true)
	SET (LIBPNG png)
	SET (LIBZ z)
	SET (LIBPTHREAD pthread)

TARGET_LINK_LIBRARIES (wxGenomeSIMLA simulation simla simpen ga genetics utility random ${wxWidgets_LIBRARIES} pngwriter ${LIBPNG} ${LIBZ} ${LIBPTHREAD} ${FREETYPE_LIBRARIES})
MESSAGE ("** Source DIR: ${PROJECT_SOURCE_DIR}")





