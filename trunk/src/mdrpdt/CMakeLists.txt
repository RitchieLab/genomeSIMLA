###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
include_directories (${Parent} ${Parent}/galib)

SET (AppFiles
	eseconfiguration.cpp 
	evalbalancedaccuracy.cpp 
	evalbalancedaccuracypdt.cpp 
	evalmaxdifference.cpp 
	gtlineparsermdrpdt.cpp 
	locuslogascii.cpp 
	matchedoddsratio.cpp 
	modelfinder.cpp 
	pdtfold.cpp 
	reportlog.cpp 
	reportlogascii.cpp 
	reportlogbinary.cpp 
	sibshipfoldproduction.cpp 
	snpsearchapplication.cpp
)

SET (PowerFiles 
	eseconfiguration.cpp 
	evalbalancedaccuracy.cpp 
	evalbalancedaccuracypdt.cpp
	evalmaxdifference.cpp
	gtlineparsermdrpdt.cpp
	locuslogascii.cpp
	matchedoddsratio.cpp
	pdtfold.cpp
	powerstudy.cpp
	reportlog.cpp
	reportlogascii.cpp
	reportlogbinary.cpp
	sibshipfoldproduction.cpp
	snpsearchapplication.cpp
)


ADD_EXECUTABLE (mdr-pdt ${AppFiles})
TARGET_LINK_LIBRARIES (mdr-pdt genetics utility random)
SET_TARGET_PROPERTIES( mdr-pdt PROPERTIES DEBUG_POSTFIX d)

ADD_EXECUTABLE (mdr-pdt-power ${PowerFiles})
TARGET_LINK_LIBRARIES (mdr-pdt-power genetics utility random)
SET_TARGET_PROPERTIES (mdr-pdt-power PROPERTIES DEBUG_POSTFIX d)
