bin_PROGRAMS=genomeSIMLA

genomeSIMLA_SOURCES=config.cpp  \
	genomesim.cpp \
	config.h \
	genomesim.h \
	timestamp.h \
	defaults.h
	
AM_CPPFLAGS=-I$(top_srcdir)/src $(FREETYPE_CPPFLAGS)
	
genomeSIMLA_LDADD=../simulation/libsimulation.la ../simla/libsimla.la \
	../simpen/libsimpen.la ../genetics/libgenetics.la ../utility/libutility.la \
	../ga/libga.la ../pngwriter/libpngwriter.la ../random/librandom.la 
