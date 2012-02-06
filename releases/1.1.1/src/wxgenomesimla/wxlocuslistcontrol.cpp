//
// C++ Implementation: wxlocuslistcontrol
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "wxlocuslistcontrol.h"

namespace GenomeSIM {

namespace GUI {


wxLocusListControl::~wxLocusListControl()
{
}


wxString wxLocusListControl::OnGetItemText(long item, long column) const {
	assert(item < loci->size());
	wxString textValue;
	switch (column) {
		case 0:
			textValue.Printf("%d", (int)item + 1);
			break;
		case 1:
			textValue = (*loci)[item]->GetLabel().c_str();	
			break;
		case 2:
			textValue.Printf("%f", (*loci)[item]->Freq1());
			break;
		case 3:
			textValue.Printf("%d", (*loci)[item]->GetLocation());
			break;
		case 4:
			textValue = wxT(((*loci)[item]->GetDescription().c_str()));
			break;
	}
	return textValue;
}

int wxLocusListControl::OnGetItemColumnImage(long item, long column) const {
	return 0;
}

wxListItemAttr *wxLocusListControl::OnGetItemAttr(long item) const {
	return NULL;
	return item % 2 ? NULL : (wxListItemAttr*)&altRowAttr;
}

}

}
