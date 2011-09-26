# - Find the native SOCI includes and library
#
# This module defines
#  SOCI_INCLUDE_DIR, where to find soci.h, etc.
#  SOCI_LIBRARIES, the libraries to link against to use SOCI.
#  SOCI_FOUND, If false, do not try to use SOCI.
# also defined, but not for general use are
#  SOCI_LIBRARY, where to find the SOCI library.

FIND_PATH(SOCI_INCLUDE_DIR soci.h /usr/include/soci /usr/local/include/soci)
FIND_LIBRARY(SOCI_LIBRARY_SQLITE3 NAMES soci_sqlite3-gcc-2_2 PATH /usr/lib /usr/local/lib)
FIND_LIBRARY(SOCI_LIBRARY_CORE NAMES soci_core-gcc-2_2 PATH /usr/lib /usr/local/lib)
IF (SOCI_INCLUDE_DIR AND SOCI_LIBRARY_SQLITE3 AND SOCI_LIBRARY_CORE)
	SET (SOCI_FOUND TRUE)
	SET (SOCI_LIBRARY ${SOCI_LIBRARY_SQLITE3} ${SOCI_LIBRARY_CORE})
ENDIF (SOCI_INCLUDE_DIR AND SOCI_LIBRARY_SQLITE3 AND SOCI_LIBRARY_CORE)

MESSAGE ("SOCI_INCLUDE_DIR ${SOCI_INCLUDE_DIR}")
MESSAGE ("SOCI_LIBRARY ${SOCI_LIBRARY}")
SET(SOCI_NAMES ${SOCI_NAMES} soci libsoci)
FIND_LIBRARY(SOCI_LIBRARY NAMES ${SOCI_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set SOCI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SOCI DEFAULT_MSG SOCI_LIBRARY SOCI_INCLUDE_DIR)

MARK_AS_ADVANCED(SOCI_INCLUDE_DIR SOCI_LIBRARY )
