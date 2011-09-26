//
// C++ Implementation: wxuser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "wxuser.h"

namespace GenomeSIM {

namespace GUI {


void WxUser::UpdateTextField(wxTextCtrl *ctrl, int val) {
	wxString value;
	value.Printf(wxT("%d"), val);
	ctrl->SetValue(value);
}

void WxUser::UpdateTextField(wxTextCtrl *ctrl, double val) {
	wxString value;
	value.Printf(wxT("%f"), val);
	ctrl->SetValue(value);
}

int WxUser::ToInt(wxString v) {
	long i;
	v.ToLong(&i);
	return i;
}

int WxUser::ExtractInteger(wxTextCtrl *ctrl) {
	long v;
	ctrl->GetValue().ToLong(&v);
	return v;
}

double WxUser::ToDouble(wxString v) {
	double i;
	v.ToDouble(&i);
	return i;
}

double WxUser::ExtractDouble(wxTextCtrl *ctrl) {
	double v;
	ctrl->GetValue().ToDouble(&v);
	return v;
}

}

}
