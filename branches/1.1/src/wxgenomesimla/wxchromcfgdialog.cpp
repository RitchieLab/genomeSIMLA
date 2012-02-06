//
// C++ Implementation: wxchromcfgdialog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "wxchromcfgdialog.h"

namespace GenomeSIM {
namespace GUI {

	wxChromCfgDialog::wxChromCfgDialog() : hasChanged(false), chromosome(NULL) { }
	
	wxChromCfgDialog::~wxChromCfgDialog() { }
	
	string wxChromCfgDialog::Initialize(int id) {
		SetID(id);
		char label[2048];
		sprintf(label, "Chromosome-%d", id);
		InitLabel(label, id);
		return string(label);
	}
	
	void wxChromCfgDialog::SetID(int id) { pageID=id; }
	
	void wxChromCfgDialog::UpdateLabel(const char *label) {
		wxNotebook* mommy = (wxNotebook*)(((wxWindow*)this)->GetParent());
		if (mommy) 
			mommy->SetPageText(mommy->GetSelection(), wxT(label));
	}
}

}
