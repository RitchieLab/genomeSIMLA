###########################################################################
# Genetics Library used for Ritchie Lab software 
#
# This library contains various generic classes and tools used by several 
# applications produced at the Ritchie Lab. 
#
#
# 

set (GeneticsFiles 
	allelicdecoder.cpp 
	binarygenoparser.cpp 
	eliminategenotypeerrors.cpp 
	familyevaluation.cpp 
	familynode.cpp 
	familymember.cpp 
	familyrepoevaluationbasic.cpp 
	familyrepoevaluation.cpp 
	familyrepository.cpp 
	fixtrios.cpp 
	genotypedata.cpp 
	gtevaluation.cpp 
	gtfileparser.cpp 
	gtfileparserbuffered.cpp 
	gtlineparser.cpp 
	modelstatistics.cpp 
	ptestdistribution.cpp 
	resultsrepository.cpp 
	snpaligned.cpp  
	snpevalsuite.cpp 
	snppool.cpp 
	snpcontainer.cpp 
	snprecipient.cpp 
	snprepository.cpp 
	snprepostxtsorted.cpp 
	snpverificationmethod.cpp 
)

add_library (genetics ${GeneticsFiles})
SET_TARGET_PROPERTIES( genetics PROPERTIES DEBUG_POSTFIX d)
