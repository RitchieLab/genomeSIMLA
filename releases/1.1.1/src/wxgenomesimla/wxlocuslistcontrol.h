//
// C++ Interface: wxlocuslistcontrol
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXLOCUSLISTCONTROL_H
#define GENOMESIM_GUIWXLOCUSLISTCONTROL_H

#include <wx/listctrl.h>
#include <vector>
#include "simulation/locus.h"

namespace GenomeSIM {

namespace GUI {

using namespace std;
using namespace Simulation;


/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class wxLocusListControl: public wxListCtrl
{
public:
    wxLocusListControl(wxWindow *parent, const wxWindowID id, const wxPoint& pos,
               const wxSize& size, long style);
    ~wxLocusListControl();
	void Initialize(vector<Locus *> *loci, int idx);

private:
    wxString OnGetItemText(long item, long column) const;
    int OnGetItemColumnImage(long item, long column) const;
    wxListItemAttr *OnGetItemAttr(long item) const;

protected:
	vector<Locus *> *loci;
	wxListItemAttr altRowAttr;

};


inline
wxLocusListControl::wxLocusListControl(wxWindow *parent, const wxWindowID id, 
	const wxPoint& pos, const wxSize& size, long style) 
	: wxListCtrl(parent, id, pos, size, style), 
	altRowAttr(*wxBLACK, *wxLIGHT_GREY, wxNullFont)       {  
	
}

inline
void wxLocusListControl::Initialize(vector<Locus *> *loci, int idx) {
	this->loci = loci;

	if (loci)
		SetItemCount(loci->size());
	InsertColumn(0, _T("Index"));
	InsertColumn(1, _T("Label"));
	InsertColumn(2, _T("Freq (A1)"));
	InsertColumn(3, _T("Position"));
	InsertColumn(4, _T("Description"));
	
}

}

}

#endif
