//
// C++ Implementation: wxgridloci
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "wxgridloci.h"
#include <fstream>
#include <sstream>
#include <iomanip>


namespace GenomeSIM {

namespace GUI {

int wxGridLoci::GetNumberCols() {
	return 6;
}

int wxGridLoci::GetNumberRows() {
	cout<<"GetNumberRows() : "<<loci->size()<<"\n";
	return loci->size();
}


void wxGridLoci::Reset() {
	vector<Locus *>::iterator itr=loci->begin();
	vector<Locus *>::iterator end=loci->end();
	
	while (itr != end) {
		delete *itr;
		itr++;
	}

	ClearTable();
	loci->clear();
}

void wxGridLoci::ClearTable(size_t lociToBeCleared) {
	if (lociToBeCleared == 0)
		lociToBeCleared = loci->size();

	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, lociToBeCleared);
		GetView()->ProcessTableMessage(msg);
	}
}	


void wxGridLoci::Refresh() {
	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_INSERTED, 0, loci->size());
		GetView()->ProcessTableMessage(msg);
	}
}

bool wxGridLoci::IsEmptyCell( int row, int col ) {
	return false;
}

wxString wxGridLoci::GetValue( int row, int col ) {	
	wxString rVal = wxT("");


	if ((size_t)row < loci->size()) {
		Locus *l = (*loci)[row];
	
//		cout<<"Fetching	"<<row<<"x"<<col<<" from SNP\n";

		switch (col) {
			case 0:
				rVal = wxT(l->GetLabel().c_str());
				break;
			case 1:			
				rVal.Printf("%f", l->Freq1());
				break;
			case 2:
				rVal.Printf("%f", l->Freq2());
				break;
			case 3:
				rVal.Printf("%f", l->MapDistance());
				break;
			case 4:
				rVal.Printf("%d", l->GetLocation());
				break;
			case 5:
				rVal = wxT(l->GetDescription().c_str());
				break;
			default:
				cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
				assert(0);
		}
	}
	return rVal;
}



void wxGridLoci::SetValue( int row, int col, const wxString& value ) {
	Locus *l = (*loci)[row];
	double temp;
	switch (col) {
		case 0:
			hasChanged = true;
			l->SetLabel( value.c_str() );
			break;
		case 1:
			if (value.ToDouble(&temp)) {
				hasChanged = true;
				l->SetFreq1(temp);
			}
			break;
		case 2:
			if (value.ToDouble(&temp)) {
				hasChanged = true;
				l->SetFreq2(temp);
			}
			break;
		case 3:
			if (value.ToDouble(&temp)) {
				hasChanged = true;
				l->SetMapDistance(temp);
			}
			break;
		case 4:
			long t;
			if (value.ToLong(&t)) {
				hasChanged =  true;
				l->SetLocation(t);
			}
			break;
		case 5:
			l->SetDescription( value.c_str());
			break;
		default:
			cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
			assert(0);
	}

}

wxString wxGridLoci::GetColLabelValue( int col ) {
	wxString rVal;
	switch (col) {
		case 0:
			rVal = wxT("Label");
			break;
		case 1:			
			rVal = wxT("Freq (Al1)");
			break;
		case 2:
			rVal = wxT("Freq (Al2)");
			break;
		case 3:
			rVal = wxT("Map Distance (cM)");
			break;
		case 4:
			rVal = wxT("Position Mb");
			break;
		case 5:
			rVal = wxT("Description");
			break;
		default:
			cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
			assert(0);
	}
	return rVal;

}

wxString wxGridLoci::GetTypeName( int row, int col ) {
	switch (col) {
		case 0:
		case 5:
			return wxGRID_VALUE_STRING;
		case 4:
			return wxGRID_VALUE_NUMBER;
		case 1:		
		case 2:
		case 3:
			return wxGRID_VALUE_FLOAT;
		default:
			cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
			assert(0);
	}
	return "";

}


/*wxGridCellAttr *wxGridLoci::GetAttr(int row, int col) {
	return attr;
}
*/
bool wxGridLoci::CanGetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)row >= loci->size())
		return false;

	if (typeName == wxGRID_VALUE_STRING)
		return true;
	else if (typeName == wxGRID_VALUE_NUMBER && col ==4)
		return true;
	else if (typeName == wxGRID_VALUE_FLOAT &&  col > 0 && col < 4)
		return true;
	else 
		return false;
}

bool wxGridLoci::CanSetValueAs( int row, int col, const wxString& typeName ) {
	if (col != 3 && (size_t)row < loci->size())
		return CanGetValueAs(row, col, typeName);
	else
		return false;
}

long wxGridLoci::GetValueAsLong( int row, int col ) {
	Locus *l = (*loci)[row];
	
	assert(col == 4);
	return l->GetLocation();
}


double wxGridLoci::GetValueAsDouble( int row, int col ) {
	Locus *l = (*loci)[row];
	
	switch (col) {
		case 1:
			return l->Freq1();
		case 2:
			return l->Freq2();
		case 3:
			return l->MapDistance();
	}
	return false;
}

void wxGridLoci::SetValueAsLong( int row, int col, long value ) {
	Locus *l = (*loci)[row];
	assert(col == 4);
	l->SetLocation(value);
}

void wxGridLoci::SetValueAsDouble( int row, int col, double value ) {
	Locus *l = (*loci)[row];
	
	switch (col) {
		case 1:
			l->SetFreq1(value);
			l->SetFreq2(1.0 - value);
			Refresh();
			break;
		case 2:
			l->SetFreq2(value);
			l->SetFreq1(1.0 - value);
			Refresh();
			break;
		case 3:
			l->SetMapDistance(value);
			break;
		default:
			cout<<"Invalid Index assignment: "<<row<<" x "<<col<<"\n";
			assert(0);
	}
}


}

}
