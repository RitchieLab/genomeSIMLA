##################################################################
# This is just a simple inclusion that will add Freetype to a
# CMake application
#
#   SET(wxWidgets_USE_LIBS base core gl net)
#   FIND_PACKAGE(wxWidgets)
#   IF(wxWidgets_FOUND)
#     INCLUDE(${wxWidgets_USE_FILE})
#     # and for each of your dependant executable/library targets:
#     TARGET_LINK_LIBRARIES(<YourTarget> ${wxWidgets_LIBRARIES})
#   ENDIF(wxWidgets_FOUND)

MESSAGE ("- Searching for WxWidgets")
SET (wxWidgets_USE_LIBS aui richtext adv html core xml base xrc qa net)
SET (wxWidgets_FIND_COMPONENTS true)
FIND_PACKAGE(wxWidgets)
IF (wxWidgets_FOUND)
	INCLUDE (${wxWidgets_USE_FILE})
	MESSAGE ("wx Libs:  ${wxWidgets_LIBRARIES}" )
	TARGET_LINK_LIBRARIES(${ProjectName} ${wxWidgets_LIBRARIES})
ENDIF (wxWidgets_FOUND)

