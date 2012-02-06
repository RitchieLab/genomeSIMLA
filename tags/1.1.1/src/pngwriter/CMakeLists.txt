###########################################################################
# PNG Writer Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
MESSAGE ("---------------------PNG WRITER---------------")
set (PngWriterFiles 
	pngwriter.cpp
)

add_library (pngwriter ${PngWriterFiles})
SET_TARGET_PROPERTIES( pngwriter PROPERTIES DEBUG_POSTFIX d)
