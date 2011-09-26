###########################################################################
# Utility Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 
SET (ARCHIVE_OUTPUT_DIRECTORY "${BuildDIR}/utility")
MESSAGE ("Archive Output Directory: 	${ARCHIVE_OUTPUT_DIRECTORY}")
SET (UtilityFiles 
	application.cpp 		application.h
							array2d.h
	basiclog.cpp 			basiclog.h
	binarrayparser.cpp 		binarrayparser.h
	casecontrolstatus.cpp 	casecontrolstatus.h
	configurationparser.cpp configurationparser.h
	executionlog.cpp 		executionlog.h
	genolookup.cpp 			genolookup.h
	genobplookup.cpp		genobplookup.h
	generatereport.cpp 		generatereport.h
	genotypeparser.cpp 		genotypeparser.h
	lineparser.cpp 			lineparser.h
	random.cpp 				random.h
	strings.cpp 			strings.h
	stringlookup.cpp 		stringlookup.h
	typelookup.cpp			typelookup.h
)	

ADD_LIBRARY (utility ${UtilityFiles})
SET_TARGET_PROPERTIES( utility PROPERTIES DEBUG_POSTFIX d)

