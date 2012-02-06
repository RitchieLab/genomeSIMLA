##################################################################
# This is just a simple inclusion that will add Freetype to a
# CMake application
#

# FindPackage(Freetype) isn't working on windows
IF (WIN32) 
	MESSAGE ("Using canned paths for freetype. If you have it installed elsewhere, you will want to edit FreeType.cmake and rerun cmake before continuing")
	SET (FREETYPE_INCLUDE_DIRS /c/unx/3rdparty/include /c/unx/3rdparty/include/freetype2)
	SET (FREETYPE_LIBRARIES libfreetype.a)
	#SET (FREETYPE_LIBRARIES freetype)
ELSE (WIN32)
	MESSAGE ("Attempting to add FreeType support")
	FIND_PACKAGE (Freetype REQUIRED)
ENDIF (WIN32)

MESSAGE ("-   ${FREETYPE_INCLUDE_DIRS}")
MESSAGE ("-   ${FREETYPE_LIBRARIES}")	
