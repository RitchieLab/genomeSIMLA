//
// C++ Interface: familyrepoevaluationbasic
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFAMILYREPOEVALUATIONBASIC_H
#define GENETICS_EVALUATIONFAMILYREPOEVALUATIONBASIC_H
#include "familyevaluation.h"
#include "familyrepoevaluation.h"

namespace Genetics {

namespace Evaluation {

/**
@brief Base class for evaluation methods that should iterate over individuals that are found within families inside a repository

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class FamilyRepoEvaluationBasic : public FamilyRepoEvaluation, public FamilyEvaluation {
public:
	FamilyRepoEvaluationBasic() {}
	~FamilyRepoEvaluationBasic() {}
	bool Evaluate(FamilyNode *family, uint position);

protected:
	FamilyNode *family;
};


inline
bool FamilyRepoEvaluationBasic::Evaluate(FamilyNode *family, uint position) {
	this->family = family;
	repoCount++;
	return family->PerformEvaluation(this);
}

}

}

#endif
