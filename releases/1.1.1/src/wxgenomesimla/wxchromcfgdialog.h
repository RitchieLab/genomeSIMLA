//
// C++ Interface: wxchromcfgdialog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXCHROMCFGDIALOG_H
#define GENOMESIM_GUIWXCHROMCFGDIALOG_H


#include <wx/notebook.h>
#include <sstream>
#include "simulation/chrompool.h"
#include "appinterface.h"

namespace GenomeSIM {
namespace GUI {

using namespace std;
using namespace Simulation;

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class wxChromCfgDialog : public WxUser {
public:
    wxChromCfgDialog();
    virtual ~wxChromCfgDialog();

	void SetID(int id);
	void UpdateLabel(const char *label);
	string Initialize(int id);
	virtual void InitLabel(const char *label, uint chromID) = 0;
	virtual void RefreshSize() = 0;
	virtual void Commit() = 0;
	virtual bool HasChanged() { return hasChanged; }

	void SetPool(ChromPool *chromosome);
	ChromPool *GetPool();
protected:
	bool hasChanged;
	int pageID;
	ChromPool *chromosome;
	
};

inline
void wxChromCfgDialog::SetPool(ChromPool *chromosome) {
	this->chromosome = chromosome;
}

inline
ChromPool *wxChromCfgDialog::GetPool() {
	return chromosome;
}

}
}

#endif
