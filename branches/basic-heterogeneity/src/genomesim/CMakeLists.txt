###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
MESSAGE ("-------------------- GenomeSIMLA (Command Line)--------------------")

INCLUDE (${PROJECT_SOURCE_DIR}/FreeType.cmake)


SET (GenomeSIMLAFiles
	appinterface.l
	config.cpp 
	genomesim.cpp 
)

SET (WrapperMainFiles 
	genomesimlawrapper.cpp
)

IF (USE_XY)
	SET_SOURCE_FILES_PROPERTIES (${GenomeSIMLAFiles} 
			PROPERTIES COMPILE_FLAGS "-DUSE_XY=1")
ENDIF (USE_XY)


ADD_EXECUTABLE (genomeSIMLA ${GenomeSIMLAFiles})
SET_TARGET_PROPERTIES( genomeSIMLA PROPERTIES DEBUG_POSTFIX d)# PROPERTIES DEBUG_POSTFIX d)
TARGET_LINK_LIBRARIES (genomeSIMLA simulation pngwriter simla simpen ga genetics utility random pngwriter ${LIBPNG} z pthread ${FREETYPE_LIBRARIES} )
IF (INCREMENT_BUILD) 
	MESSAGE ("ADDING increment build command")
	ADD_CUSTOM_COMMAND ( OUTPUT appinterface.l
		COMMAND inc_version.sh ARGS ${PROJECT_SOURCE_DIR}/src/genomesim/timestamp.h
		COMMENT "Incrementing build"
	)
	SET (genomeSIMLA ${GenomeSIMLAFiles} appinterface.l)
ENDIF (INCREMENT_BUILD)

#ADD_EXECUTABLE (wrapperMains ${WrapperMainFiles})
#TARGET_LINK_LIBRARIES (wrapperMains pthread)