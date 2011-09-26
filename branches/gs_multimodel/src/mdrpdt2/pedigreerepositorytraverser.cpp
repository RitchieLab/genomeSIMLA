//
// C++ Implementation: pedigreerepositorytraverser
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreerepositorytraverser.h"
#include <sstream>

using namespace std;

namespace MdrPDT {

namespace Validator {

PedigreeRepositoryTraverser::PedigreeRepositoryTraverser()
{
}


PedigreeRepositoryTraverser::~PedigreeRepositoryTraverser()
{
}


void PedigreeRepositoryTraverser::EvaluateRepository(PedigreeRepository *repo) {
	PedigreeRepository::Iterator itr = repo->GetIterator();
	Pedigree *ped = itr.GetNext();
	PrepEvaluation();

	stringstream summary;
	while (ped) {
		_Evaluate(ped);
		ped = itr.GetNext();
	}
	
	PostEvaluation();
}

void PedigreeRepositoryTraverser::_Evaluate(Pedigree *ped) {
	Pedigree::Iterator itr = ped->GetIterator();
	Individual *ind = itr.GetNext();
		
	Evaluate(ped);

	while (ind) {
		Evaluate(ped, ind);
		ind = itr.GetNext();
	}
}

}

}
