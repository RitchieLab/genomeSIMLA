
#
# CMake makefile.
#
# This file is part of the build system used to build 
# genomeSIMLA (cmd and gui) and mdr-pdt on windows, mac
# and windows. 
#
# If you don't have cmake 2.7 installed, you will need to
# Download and install it freely at:
# www.cmake.org
IF (GCC34)
	SET(CXX g++34)
ENDIF (GCC34)

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
CMAKE_POLICY (SET CMP0003 NEW)

SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/src)
MESSAGE ("CMAKE MODULE PATH: ${CMAKE_MODULE_PATH}")

SET (CMAKE_FIND_LIBRARY_PREFIXES lib)
SET (CMAKE_FIND_LIBRARY_SUFFIXES .a)
FIND_LIBRARY(LIBPNG png PATHS /usr/local/lib /usr/lib/ /usr/X11R6/lib NO_DEFAULT_PATH
NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH)
IF (LIBPNG-NOTFOUND)
	MESSAGE("Unable to find libpng. Unable to continue")
ELSE (LIBPNG-NOTFOUND)
	MESSAGE("LIBPNG: ${LIBPNG}")
ENDIF (LIBPNG-NOTFOUND)

MESSAGE ("Build Types Available: ${CMAKE_CONFIGURATION_TYPES}")
SET (ProjectName GenomeSIMLA)
#SET (MAKE_GENOMESIM true)

IF (MAKE_GENOMESIM)	
	SET (ProjectName wxGenomeSIMLA)
	INCLUDE(FreeType.cmake)
ELSEIF (MAKE_MDRPDT)
	SET (ProjectName mdr-pdt)
ELSEIF (MAKE_MDRPDT2)
	SET (ProjectName mdr-pdt2)
ELSEIF (MAKE_BIOFILTER)
	MESSAGE ("Making Biofilter")
	SET (ProjectName biofilter)
	SET (USE_SOCI 1)
ELSE (MAKE_GENOMESIM)
	ADD_DEFINITIONS(-DCPPUNIT)
	SET (ProjectName unit)
	INCLUDE(FreeType.cmake)
ENDIF (MAKE_GENOMESIM)
	
PROJECT(${ProjectName})

#SET (LIBPNG /usr/local/lib/libpng.a)

IF (APPLE) 
	SET (CMAKE_OSX_ARCHITECTURES "i386")		#x86_64 ppc ppc64
	SET (MACOSX_DEPLOYMENT_TARGET "10.4")
	SET (CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.4u.sdk")
#	SET (PLAT_CFLAGS "-I/Xcode2.5/usr/include")
#	SET (SDK "/Developer/SDKs/MacOSX10.4u.sdk")
#	SET (PLAT_CFLAGS "-isysroot ${SDK} -B /usr/local/include -I/usr/local/include -I/usr/incude")
#	SET (PLAT_LFLAGS "${PLAT_CFLAGS} -Wl,-syslibroot,${SDK}")
	SET (PlatformBuildFlags "-mmacosx-version-min=10.4 -I/usr/include ")
ENDIF(APPLE)
#-I/sw/include -L/sw/lib

SET (PLATFORM "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
SET (SourceDIR ${PROJECT_SOURCE_DIR}/src)
MESSAGE ("Cache File Directory: ${CMAKE_CACHEFILE_DIR}")
SET (BuildDIR ${CMAKE_CACHEFILE_DIR}/../bin/${PLATFORM})
SET (LibSearchPaths /usr/local /usr/local/lib)

MESSAGE ("Source Dir: 		${SourceDIR}")

SET (EXECUTABLE_OUTPUT_PATH ../..)

IF (USE_XY)
	SET (XY_FLAG "-DUSE_XY=1")
ELSE (USE_XY)
	SET (XY_FLAG "")
ENDIF (USE_XY)



IF (NOFREEZE)
	SET (INCREMENT_BUILD true)
	MESSAGE ("Setting Increment Build number to true")
ELSE (NOFREEZE)
	SET (INCREMENT_BUILD false)
ENDIF (NOFREEZE)


MESSAGE ("${XY_FLAG}")

IF (MING) 
	MESSAGE ("Setting up details for creating Windows application")
	SET (EXTERNAL_INC C:/unx/3rdparty/include)
	SET (EXTERNAL_LINK C:/unx/3rdparty/lib)
	SET (PLATFORM_EXEC_TYPE WIN32)
	SET (LIBPNG /c/unx/3rdParty/lib/libpng12.a)
ENDIF (MING)

#SET (SOCI_HOME /home/torstees/sqlite3)
#SET (SOCI_HOME /usr)
IF (USE_SOCI)
	FIND_PACKAGE(SOCI REQUIRED)
	INCLUDE (${SOCI_INCLUDE_DIR})

	SET (EXTERNAL_INC ${EXTERNAL_INC} ${SOCI_INCLUDE_DIR} ${SOCI_INCLUDE_DIR}/sqlite3 ${SOCI_INCLUDE_DIR}/mysql)
ENDIF (USE_SOCI)
#SET (CMAKE_CXX_FLAGS "-fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations ${PlatformBuildFlags} ${PLAT_CFLAGS} ${PLAT_LFLAGS}")
SET (CMAKE_CXX_FLAGS "${PlatformBuildFlags} ${PLAT_CFLAGS} ${PLAT_LFLAGS} ${XY_FLAGS}")

INCLUDE_DIRECTORIES (${SourceDIR} $ENV{rdParty}/include ${BOOST_INC} ${EXTERNAL_INC} ${PTHREAD_LIB} ${FREETYPE_INCLUDE_DIRS})
LINK_DIRECTORIES (${SOCI_LIBRARY/} $ENV{rdParty}/lib /usr/local/lib ${EXTERNAL_LINK} /usr/local/lib ${PTHREAD_LIB})

#MESSAGE (${INCLUDE_DIRECTORIES})

ADD_SUBDIRECTORY (${SourceDIR}/utility)
#ADD_SUBDIRECTORY (${SourceDIR}/genetics)
ADD_SUBDIRECTORY (${SourceDIR}/random)

IF (MAKE_BIOFILTER)
	MESSAGE("Compiling Biofilter")
	ADD_SUBDIRECTORY (${SourceDIR}/biofilter)
ENDIF (MAKE_BIOFILTER)

IF (MAKE_TESTS)
	MESSAGE("Compiling Unit Tests")
	ADD_SUBDIRECTORY (${SourceDIR}/simpen)
	ADD_SUBDIRECTORY (${SourceDIR}/ga)
	ADD_SUBDIRECTORY (${SourceDIR}/simulation)
	ADD_SUBDIRECTORY (${SourceDIR}/simla)
ENDIF(MAKE_TESTS)

IF (MAKE_GENOMESIM)	
	MESSAGE ("Compiling GenomeSIMLA")
	ADD_SUBDIRECTORY (${SourceDIR}/pngwriter)
	ADD_SUBDIRECTORY (${SourceDIR}/genomesim)
	ADD_SUBDIRECTORY (${SourceDIR}/simpen)
	ADD_SUBDIRECTORY (${SourceDIR}/ga)
	ADD_SUBDIRECTORY (${SourceDIR}/genetics)
	ADD_SUBDIRECTORY (${SourceDIR}/simulation)
	ADD_SUBDIRECTORY (${SourceDIR}/simla)
IF (NOT SKIP_GUI)
	ADD_SUBDIRECTORY (${SourceDIR}/wxgenomesimla)
ENDIF (NOT SKIP_GUI)

ELSEIF (MAKE_MDRPDT)
	ADD_SUBDIRECTORY (${SourceDIR}/mdrpdt)
	LINK_LIBRARIES (z pthread)
ELSEIF (MAKE_MDRPDT2)
	ADD_SUBDIRECTORY (${SourceDIR}/mdrpdt2)
	LINK_LIBRARIES (z pthread)
ENDIF (MAKE_GENOMESIM)
