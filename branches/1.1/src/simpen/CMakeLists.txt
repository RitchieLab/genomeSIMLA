###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
set (SimpenFiles 
	prefexcept.cpp 
	simpen.cpp 
	simpenmath.cpp 
	simpenmdr.cpp 
	simpreferences.cpp 
)

ADD_LIBRARY (simpen ${SimpenFiles})
SET_TARGET_PROPERTIES( simpen PROPERTIES DEBUG_POSTFIX d)
