//
// C++ Interface: gtevaluation
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONGTEVALUATION_H
#define GENETICS_EVALUATIONGTEVALUATION_H
#include "genotypedata.h"

namespace Genetics {

namespace Evaluation {


/**
@brief Functors used to evaluate an entry in a repository

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GtEvaluation{
public:
    GtEvaluation();
    virtual ~GtEvaluation();

	virtual uint GetGtCount();
	virtual void ResetGtCount();
	virtual bool Evaluate(GenotypeData *snp, uint position) = 0;	
protected:
	uint genotypeCount;

};

	

inline
GtEvaluation::GtEvaluation() : genotypeCount(0) {}

inline
GtEvaluation::~GtEvaluation() {}

inline
uint GtEvaluation::GetGtCount() {
	return genotypeCount;
}

inline
void GtEvaluation::ResetGtCount() {
	genotypeCount=0;
}
}
}
#endif
