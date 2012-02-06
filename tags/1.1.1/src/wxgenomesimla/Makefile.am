bin_PROGRAMS=wxGenomeSIMLA

wxGenomeSIMLA_SOURCES=appconfig.cpp \
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
	
AM_CPPFLAGS=-I.. $(FREETYPE_CPPFLAGS)

wxGenomeSIMLA_LDADD=../simulation/libsimulation.la ../simla/libsimla.la \
	../simpen/libsimpen.la ../genetics/libgenetics.la ../utility/libutility.la \
	../ga/libga.la ../pngwriter/libpngwriter.la ../random/librandom.la
