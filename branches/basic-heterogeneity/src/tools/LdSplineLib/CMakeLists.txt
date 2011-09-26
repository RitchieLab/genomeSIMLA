PROJECT(ldlibrary)

MESSAGE ("Creating makefile required for the LDSpling library")
SET (CMAKE_VERBOSE_MAKEFILES)

SET (CMAKE_CXX_FLAGS "-fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations")

# PATH FOR MYSQL INCLUDE DIR
INCLUDE_DIRECTORIES(/usr/include/mysql /usr/include/)

#ADD_DEFINITIONS("-DHAVE_DLOPEN")
#ldlibrary.def

SET ( Sources
	ldspline.cpp
	ldstringency.cpp
)

ADD_LIBRARY(ldspline SHARED ${Sources} )

