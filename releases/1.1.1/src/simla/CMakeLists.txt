###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
set (SimlaFiles 
	dcdflib.cpp 
	penetrance.cpp	
	random.cpp 
)

ADD_LIBRARY (simla ${SimlaFiles})
SET_TARGET_PROPERTIES( simla PROPERTIES DEBUG_POSTFIX d)
