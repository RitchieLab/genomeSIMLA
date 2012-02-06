//
// C++ Interface: familyrepoevaluation
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFAMILYREPOEVALUATION_H
#define GENETICS_EVALUATIONFAMILYREPOEVALUATION_H
#include "utility/utility.h"
#include "familynode.h"

namespace Genetics {
namespace Evaluation {


class FamilyRepoEvaluation {
public:
	FamilyRepoEvaluation();
	virtual ~FamilyRepoEvaluation();
	/**
	 * @brief Returns the number of repositories have been iterated over (such as families)
	 */
	virtual uint GetRepoCount();
	virtual void ResetRepoCount();
	virtual bool Evaluate(FamilyNode *family, uint position)=0;	
protected:
	uint repoCount;
};


inline
FamilyRepoEvaluation::FamilyRepoEvaluation() : repoCount(0)
{
}

inline
FamilyRepoEvaluation::~FamilyRepoEvaluation()
{
}


inline
uint FamilyRepoEvaluation::GetRepoCount() {
	return repoCount;
}

inline
void FamilyRepoEvaluation::ResetRepoCount() {
	repoCount=0;
}

}

}

#endif
