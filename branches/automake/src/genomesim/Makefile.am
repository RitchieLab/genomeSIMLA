bin_PROGRAMS=genomeSIMLA

genomeSIMLA_SOURCES=config.cpp  \
	genomesim.cpp
	
AM_CPPFLAGS=-I.. $(FREETYPE_CPPFLAGS)
	
genomeSIMLA_LDADD=../simulation/libsimulation.la ../simla/libsimla.la \
	../simpen/libsimpen.la ../genetics/libgenetics.la ../utility/libutility.la \
	../ga/libga.la ../pngwriter/libpngwriter.la ../random/librandom.la 
