//
// C++ Interface: simwizard
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUISIMWIZARD_H
#define GENOMESIM_GUISIMWIZARD_H

#include "appcontroller.h"
#include "wxuser.h"
namespace GenomeSIM {

namespace GUI {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SimWizard: public WxUser {
public:
	SimWizard() : appController(NULL) { }
	virtual ~SimWizard() {}
	virtual void SetAppController(AppController *appController) {
		this->appController=appController;
	}
	//virtual void Initialize() = 0;

protected:
	AppController *appController;

};

}

}

#endif
