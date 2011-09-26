//
// C++ Interface: wxuser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXUSER_H
#define GENOMESIM_GUIWXUSER_H
#include "wx/textctrl.h"
namespace GenomeSIM {

namespace GUI {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class WxUser {
public:
	void UpdateTextField(wxTextCtrl *ctrl, int value);
	void UpdateTextField(wxTextCtrl *ctrl, double value);

	
	int ToInt(wxString s);
	int ExtractInteger(wxTextCtrl *ctrl);
	
	double ToDouble(wxString s);
	double ExtractDouble(wxTextCtrl *ctrl);
};

}

}

#endif
