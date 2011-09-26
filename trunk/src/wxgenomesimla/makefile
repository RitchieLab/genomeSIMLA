PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=wxGenomeSIMLA
HDRS=appconfig.h \
	appcontroller.h \
	appinterface.h \
	diseasemodeldetails.h \
	dlggrowthreview.h \
	genomesimla.h \
	growthrateconfig.h \
	locusmanager.h \
	mainframe.h \
	modellist.h \
	populationcontroller.h \
	treelistctrl.h \
	wxbeginsimulation.h \
	wxchromcfgdialog.h \
	wxdialogtaskalert.h\
	wxdlgconfigureccdatasets.h \
	wxdlgconfigurediseasemodel.h \
	wxdlgconfigurepeddatasets.h \
	wxdlgeditlocusselector.h \
	wxdlgbasichtmlreport.h \
	wxdlgprocessmonitor.h \
	wxdlgsavecurrent.h \
	wxdlgselectlocusrange.h \
	wxdlgsimla.h \
	wxgridloci.h \
	wximggrowthchart.h \
	wxldanalysis.h \
	wxlocuslistcontrol.h \
	wxpageinfo.h \
	wxpagebbchrom.h \
	wxpagegeneralsettings.h \
	wxpageloadchromosome.h \
	wxpagelocicreation.h \
	wxpagelocusreporting.h \
	wxpagepenetrancemodel.h \
	wxpagereporting.h \
	wxpageselectgeneration.h \
	wxpagesimpen.h \
	wxpanelchromosomerep.h \
	wxpanelreports.h \
	wxuser.h \
	wxwizardsimulatedata.h \
	wxwizardcreatechromosome.h \
	wxwizpagebbchrom.h \
	wxwizpagedefinedatasets.h \
	wxwizpagelaunchconfirmation.h \
	wxwizpagelaunchdroppoints.h \
	wxwizpagelaunchproject.h \
	wxwizpagelaunchselectgeneration.h \
	wxwizpagelaunchsimulation.h \
	wxwizpagereview.h \
	wxwizpagereviewsimulation.h \
	wxwizpageselectgeneration.h \
	wxwizpageselectloci.h \
	wxwizrunsimulation.h
#DEFINE THE sources
SRCS=appconfig.cpp \
	appcontroller.cpp \
	appinterface.cpp \
	diseasemodeldetails.cpp \
	dlggrowthreview.cpp \
	genomesimla.cpp \
	growthrateconfig.cpp \
	locusmanager.cpp \
	mainframe.cpp \
	modellist.cpp \
	populationcontroller.cpp \
	treelistctrl.cpp \
	wxbeginsimulation.cpp \
	wxchromcfgdialog.cpp \
	wxdialogtaskalert.cpp\
	wxdlgconfigureccdatasets.cpp \
	wxdlgconfigurediseasemodel.cpp \
	wxdlgconfigurepeddatasets.cpp \
	wxdlgeditlocusselector.cpp \
	wxdlgbasichtmlreport.cpp \
	wxdlgprocessmonitor.cpp \
	wxdlgsavecurrent.cpp \
	wxdlgselectlocusrange.cpp \
	wxdlgsimla.cpp \
	wxgridloci.cpp \
	wximggrowthchart.cpp \
	wxldanalysis.cpp \
	wxlocuslistcontrol.cpp \
	wxpageinfo.cpp \
	wxpagebbchrom.cpp \
	wxpagegeneralsettings.cpp \
	wxpageloadchromosome.cpp \
	wxpagelocicreation.cpp \
	wxpagelocusreporting.cpp \
	wxpagepenetrancemodel.cpp \
	wxpagereporting.cpp \
	wxpageselectgeneration.cpp \
	wxpagesimpen.cpp \
	wxpanelchromosomerep.cpp \
	wxpanelreports.cpp \
	wxuser.cpp \
	wxwizardsimulatedata.cpp \
	wxwizardcreatechromosome.cpp \
	wxwizpagebbchrom.cpp \
	wxwizpagedefinedatasets.cpp \
	wxwizpagelaunchconfirmation.cpp \
	wxwizpagelaunchdroppoints.cpp \
	wxwizpagelaunchproject.cpp \
	wxwizpagelaunchselectgeneration.cpp \
	wxwizpagelaunchsimulation.cpp \
	wxwizpagereview.cpp \
	wxwizpagereviewsimulation.cpp \
	wxwizpageselectgeneration.cpp \
	wxwizpageselectloci.cpp \
	wxwizrunsimulation.cpp

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..

#indicate what the executable will be named
EXENAME=wxGenomeSIMLA

#list the projects which build into libraries
LOCAL_LIBS=simulation simla simpen genetics utility ga random pngwriter 
#appinterface.l

#Uncomment this line for executables
#LIBRARY_MAKE=0

#Uncomment this line for test applications
IS_TEST=0
#appinterface.l: 
	#inc_version.sh timestamp.h
FREETYPE=1
WX_WIDGETS=1
MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base

#WX_CPPFLAGS := $(shell wx-config --cppflags)
#FT_CPPFLAGS := $(shell freetype-config --cflags)
#WX_LINK := $(shell wx-config --libs)

#FT_LINK := $(shell freetype-config --libs)
#EXTERNAL_LINKS=$(WX_LINK) -lpng $(FT_LINK) -lz -lpthread 
#ALT_INCLUDES=$(WX_CPPFLAGS) $(FT_CPPFLAGS)

OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



