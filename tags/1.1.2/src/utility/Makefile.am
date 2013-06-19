noinst_LTLIBRARIES=libutility.la

libutility_la_SOURCES=application.cpp \
	basiclog.cpp \
	binarrayparser.cpp \
	casecontrolstatus.cpp \
	configurationparser.cpp \
	executionlog.cpp \
	filebuffer.cpp \
	genolookup.cpp \
    genobplookup.cpp \
	generatereport.cpp \
	genotypeparser.cpp \
	lineparser.cpp \
	random.cpp \
	stringlookup.cpp \
	strings.cpp \
	typelookup.cpp \
	types.cpp\
	application.h \
	array2d.h \
	basiclog.h \
	binarrayparser.h \
	casecontrolstatus.h \
	configurationparser.h \
	exception.h \
	executionlog.h \
	filebuffer.h \
	functors.h \
	generatereport.h \
	genobplookup.h \
	genolookup.h \
	genotypeparser.h \
	lineparser.h \
	random.h\
	rbtree.h \
	stat.h \
	stringlookup.h \
	strings.h \
	typelookup.h \
	types.h \
	utility.h \
	tsrandom.h
	
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS)

libutility_la_LDFLAGS=$(BOOST_LDFLAGS)
libutility_la_LIBADD=$(BOOST_THREAD_LIB)
