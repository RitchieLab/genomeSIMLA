/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigurepeddatasets.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 21 Mar 2008 12:39:38 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "wxdlgconfigurepeddatasets.h"
#include <sstream>
#include <iostream>

////@begin XPM images
#include "img/plus.xpm"
#include "img/minus.xpm"
////@end XPM images
namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxDlgConfigurePedDatasets type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgConfigurePedDatasets, wxDialog )


/*!
 * wxDlgConfigurePedDatasets event table definition
 */

BEGIN_EVENT_TABLE( wxDlgConfigurePedDatasets, wxDialog )

////@begin wxDlgConfigurePedDatasets event table entries
    EVT_SIZE( wxDlgConfigurePedDatasets::OnSize )

    EVT_MENU( ID_BTN_ADD_FAMILY_TYPE, wxDlgConfigurePedDatasets::OnBtnAddFamilyTypeClick )

    EVT_MENU( ID_BTN_DEL_FAMILY_TYPE, wxDlgConfigurePedDatasets::OnBtnDelFamilyTypeClick )

    EVT_BUTTON( wxID_CANCEL, wxDlgConfigurePedDatasets::OnCancelClick )

    EVT_BUTTON( wxID_OK, wxDlgConfigurePedDatasets::OnOkClick )

////@end wxDlgConfigurePedDatasets event table entries

END_EVENT_TABLE()


/*!
 * wxDlgConfigurePedDatasets constructors
 */

wxDlgConfigurePedDatasets::wxDlgConfigurePedDatasets() : gridMaster(NULL)
{
    Init();
}

wxDlgConfigurePedDatasets::wxDlgConfigurePedDatasets( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : gridMaster(NULL)
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgConfigurePedDatasets creator
 */

bool wxDlgConfigurePedDatasets::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgConfigurePedDatasets creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgConfigurePedDatasets creation
    return true;
}


/*!
 * wxDlgConfigurePedDatasets destructor
 */

wxDlgConfigurePedDatasets::~wxDlgConfigurePedDatasets()
{
////@begin wxDlgConfigurePedDatasets destruction
////@end wxDlgConfigurePedDatasets destruction
}


/*!
 * Member initialisation
 */

void wxDlgConfigurePedDatasets::Init()
{
////@begin wxDlgConfigurePedDatasets member initialisation
    txtLabel = NULL;
    txtGenotypeErr = NULL;
    txtPhenocopy = NULL;
    txtMissingData = NULL;
    grdFamTypes = NULL;
////@end wxDlgConfigurePedDatasets member initialisation
}


/*!
 * Control creation for wxDlgConfigurePedDatasets
 */

void wxDlgConfigurePedDatasets::CreateControls()
{    
////@begin wxDlgConfigurePedDatasets content construction
    wxDlgConfigurePedDatasets* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxFlexGridSizer* itemFlexGridSizer3 = new wxFlexGridSizer(3, 3, 1, 1);
    itemFlexGridSizer3->AddGrowableCol(1);
    itemFlexGridSizer3->AddGrowableCol(2);
    itemBoxSizer2->Add(itemFlexGridSizer3, 1, wxGROW|wxALL, 0);

    wxStaticText* itemStaticText4 = new wxStaticText( itemDialog1, wxID_STATIC, _("Label:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText4, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtLabel = new wxTextCtrl( itemDialog1, ID_DISEASE_LABEL, _("label"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_ALPHANUMERIC) );
    itemFlexGridSizer3->Add(txtLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText10 = new wxStaticText( itemDialog1, wxID_STATIC, _("Genotype error:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText10, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtGenotypeErr = new wxTextCtrl( itemDialog1, ID_DISEASE_GT_ERROR, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtGenotypeErr, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText13 = new wxStaticText( itemDialog1, wxID_STATIC, _("Phenocopy:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText13, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtPhenocopy = new wxTextCtrl( itemDialog1, ID_DISEASE_PHENO_ERR, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtPhenocopy, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText16 = new wxStaticText( itemDialog1, wxID_STATIC, _("Missing Data:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText16, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMissingData = new wxTextCtrl( itemDialog1, ID_DISEASE_MISSING, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtMissingData, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    grdFamTypes = new wxGrid( itemDialog1, ID_FAMTYPE_GRID, wxDefaultPosition, wxSize(400, 150), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdFamTypes->SetDefaultColSize(50);
    grdFamTypes->SetDefaultRowSize(25);
    grdFamTypes->SetColLabelSize(25);
    grdFamTypes->SetRowLabelSize(50);
    grdFamTypes->CreateGrid(5, 5, wxGrid::wxGridSelectCells);
    itemBoxSizer2->Add(grdFamTypes, 1, wxGROW|wxALL, 5);

    wxToolBar* itemToolBar20 = new wxToolBar( itemDialog1, wxID_ANY, wxDefaultPosition, wxSize(-1, 30), wxTB_FLAT|wxNO_BORDER );
    itemToolBar20->SetBackgroundColour(wxColour(255, 255, 255));
    itemToolBar20->SetFont(wxFont(6, wxDEFAULT, wxNORMAL, wxNORMAL, false, wxT("Sans")));
    itemToolBar20->SetToolPacking(2);
    itemToolBar20->SetToolBitmapSize(wxSize(16, 16));
    wxBitmap itemtool21Bitmap(itemDialog1->GetBitmapResource(wxT("img/plus.xpm")));
    wxBitmap itemtool21BitmapDisabled;
    itemToolBar20->AddTool(ID_BTN_ADD_FAMILY_TYPE, _("Add Block"), itemtool21Bitmap, itemtool21BitmapDisabled, wxITEM_NORMAL, _("Add a new block to the avable block list"), wxEmptyString);
    itemToolBar20->AddSeparator();
    wxBitmap itemtool23Bitmap(itemDialog1->GetBitmapResource(wxT("img/minus.xpm")));
    wxBitmap itemtool23BitmapDisabled;
    itemToolBar20->AddTool(ID_BTN_DEL_FAMILY_TYPE, _("Remove"), itemtool23Bitmap, itemtool23BitmapDisabled, wxITEM_NORMAL, _("Delete the block"), wxEmptyString);
    itemToolBar20->Realize();
    itemBoxSizer2->Add(itemToolBar20, 0, wxALIGN_RIGHT|wxALL, 0);

    wxBoxSizer* itemBoxSizer24 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer24, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton25 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer24->Add(itemButton25, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton26 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer24->Add(itemButton26, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgConfigurePedDatasets content construction

	if (gridMaster)
		delete gridMaster;
	gridMaster = new wxGridFamilyTypes();
	grdFamTypes->SetTable(gridMaster, true);

}


/*!
 * Should we show tooltips?
 */

bool wxDlgConfigurePedDatasets::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgConfigurePedDatasets::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgConfigurePedDatasets bitmap retrieval
    wxUnusedVar(name);
    if (name == _T("img/plus.xpm"))
    {
        wxBitmap bitmap( plus_xpm);
        return bitmap;
    }
    else if (name == _T("img/minus.xpm"))
    {
        wxBitmap bitmap( minus_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxDlgConfigurePedDatasets bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgConfigurePedDatasets::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgConfigurePedDatasets icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgConfigurePedDatasets icon retrieval
}



/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
 */

void wxDlgConfigurePedDatasets::OnCancelClick( wxCommandEvent& event )
{
	if ( IsModal() )
		EndModal(wxID_CANCEL);
	else {
		SetReturnCode(wxID_CANCEL);
		this->Show(false);
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
 */

void wxDlgConfigurePedDatasets::OnOkClick( wxCommandEvent& event )
{
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}

bool wxDlgConfigurePedDatasets::FamType::Parse(const char *line) {
	stringstream input(line);
	string val="", val2="";

	input>>val>>val2;
	bool success = val=="DATASET" && val2 == "FAMTYPE";

	if (success) {
		input>>affected>>unaffected>>extrasibs>>famcount;
		//cout<<"Loaded: "<<affected<<" "<<unaffected<<" "<<extrasibs<<" "<<famcount<<"\n";
	}

	return success;
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::ConfigureModelDetails(istream &modelDetails) {
	while (!modelDetails.eof()) {
		FamType newType;
		char line[2048];
		modelDetails.getline(line, 2048);
		std::cout<<"New Line: "<<line<<"\n";
		if (newType.Parse(line)) 
			familyTypes.push_back(newType);
	}

	if (GetView()){
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, familyTypes.size());
		GetView()->ProcessTableMessage(msg);
	}

}

string wxDlgConfigurePedDatasets::GetModelLabel() {
	return txtLabel->GetLineText(0).c_str();
}

string wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetSummary() {
	stringstream summary;
	int famCount = familyTypes.size();
	int totalIndividuals = 0;
	for (int i=0; i<famCount; i++) {
		int indCount =2;
		indCount += familyTypes[i].affected;
		indCount += familyTypes[i].unaffected;
		indCount += (familyTypes[i].extrasibs / 2);
		totalIndividuals+=(indCount*familyTypes[i].famcount);
	}
	
	summary<<famCount<<" Families\t~"<<totalIndividuals<<" Individuals";
	return summary.str().c_str();
}

void wxDlgConfigurePedDatasets::ConfigureModelDetails(const char *details) {
	stringstream modelDetails(details);

	string val;	
	//Eat the first keywords, DATASET PED
	modelDetails>>val;
	modelDetails>>val;


	modelDetails>>val;
	txtLabel->SetValue(val.c_str());
		
	modelDetails>>val;
	txtGenotypeErr->SetValue(val.c_str());

	modelDetails>>val;
	txtPhenocopy->SetValue(val.c_str());
	
	modelDetails>>val;
	txtMissingData->SetValue(val.c_str());

	gridMaster->ConfigureModelDetails(modelDetails);
	
}

string wxDlgConfigurePedDatasets::GetModelSummary() {
	stringstream summary("PED\t");
	summary<<gridMaster->GetSummary()<<"\n";
	return summary.str();	
}

string wxDlgConfigurePedDatasets::GetModelDetails() {
	stringstream modelDetails;
	modelDetails<<"DATASET PED "<<txtLabel->GetLineText(0)<<" "
		<<txtGenotypeErr->GetValue()<<" "<<txtPhenocopy->GetValue()<<" "<<txtMissingData->GetValue()<<"\n";
	int familyTypeCount = gridMaster->GetNumberRows();
	
	while (familyTypeCount-- > 0) 
		modelDetails<<gridMaster->familyTypes[familyTypeCount].GetDetails()<<"\n";
	return modelDetails.str().c_str();
}



bool wxDlgConfigurePedDatasets::wxGridFamilyTypes::IsEmptyCell(int row, int col) {
	return row<(int)familyTypes.size() && col<4;
}

/*
	int affected;
	int unaffected;
	int extrasibs;
	int famcount;
*/

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::SetValue(int row, int col, const wxString& value) {
	assert(row<(int)familyTypes.size());
	
	FamType &entry = familyTypes[row];
	switch (col) {
		case 0:
			entry.affected = atoi(value.c_str());			
			break;
		case 1:
			entry.unaffected = atoi(value.c_str());
			break;
		case 2:
			entry.extrasibs = atoi(value.c_str());
			break;
		case 3:
			entry.famcount = atoi(value.c_str());
			break;
	}	
}

wxString wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetValue(int row, int col) {
	assert(row<(int)familyTypes.size());
	FamType &entry = familyTypes[row];
	wxString val;
	switch (col) { 
		case 0:
			val.Printf("%d", entry.affected);
			break;
		case 1:
			val.Printf("%d", entry.unaffected);
			break;
		case 2:	
			val.Printf("%d", entry.extrasibs);
			break;
		 case 3:
			val.Printf("%d", entry.famcount);
			break;
	}
	return val;
}

long wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetValueAsLong(int row, int col) {
	assert((size_t)row<familyTypes.size());
	FamType &entry = familyTypes[row];
	int value;

	switch (col) {
		case 0:
			return entry.affected;
		case 1:	
			return entry.unaffected;
		case 2:
			return entry.extrasibs;
		case 3:
			return entry.famcount;
	}
	return -1;
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::SetValueAsLong(int row, int col, long value) {
	assert((size_t)row<familyTypes.size());
	FamType &entry = familyTypes[row];
	
	switch(col) {
		case 0:
			entry.affected = value;
			break;
		case 1:
			entry.unaffected = value;
			 break;
		case 2:
			entry.extrasibs = value;
			break;
		case 3:
			entry.famcount = value;
			break;
	}
	return;
}
	
bool wxDlgConfigurePedDatasets::wxGridFamilyTypes::CanGetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)row >= familyTypes.size())
		return false;
	if (typeName == wxGRID_VALUE_NUMBER)
		return true;
	else 
		return false;
}

bool wxDlgConfigurePedDatasets::wxGridFamilyTypes::CanSetValueAs( int row, int col, const wxString& typeName ) {
	return CanGetValueAs(row, col, typeName);
}

wxString wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetColLabelValue(int col) {
	assert(col<4);

	switch(col) {
		case 0:
			return "Affected";
		case 1:
			return "Unaffected";
		case 2:
			return "Extra";
		case 3:
			return "Family Count";
	}
	return "";
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::AddFamilyType() {
	FamType newFamilyType;
	std::cout<<"Adding Family: "<<newFamilyType.affected<<" "<<newFamilyType.unaffected<<" "<<newFamilyType.famcount<<"\n";
	familyTypes.push_back(newFamilyType);
	std::cout<<"Family Type Count: "<<familyTypes.size()<<"\n";

	if (GetView()){
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 1);//, penetrances.size());
		GetView()->ProcessTableMessage(msg);
	}
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::DelFamilyType(int idx) {
	if (idx >= 0  && idx < familyTypes.size()) {
		vector<FamType>::iterator itr = familyTypes.begin();
		vector<FamType>::iterator end = familyTypes.end();
		int i = idx;
		while (i-- > 0) 
			itr++;
		familyTypes.erase(itr);
	if (GetView()){
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 1, idx);//, penetrances.size());
		GetView()->ProcessTableMessage(msg);
	}
	}
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::AddFamilyType(const char *data) {
	stringstream ss(data);
	int aff=0, unaff=0, extra=0, famCount=0;
	ss>>aff>>unaff>>extra>>famCount;
	familyTypes.push_back(FamType(aff, unaff, extra, famCount));
	Refresh();
}

wxString wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetRowLabelValue(int row) {
	wxString rVal;
	rVal.Printf("%d", row);
	return rVal;
}


wxString wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetTypeName( int row, int col ) {
	return wxGRID_VALUE_NUMBER;
}

string wxDlgConfigurePedDatasets::wxGridFamilyTypes::GetDetails() {	
	return "Working on it";
}

string wxDlgConfigurePedDatasets::FamType::GetDetails() {
	stringstream ss;
	ss<<"DATASET FAMTYPE "<<affected<<" "<<unaffected<<" "<<extrasibs<<" "<<famcount;
	return ss.str();
}


/*!
 * wxEVT_SIZE event handler for ID_WXDLGCONFIGUREPEDDATASETS
 */
void wxDlgConfigurePedDatasets::OnSize( wxSizeEvent& event )
{
	int width, height;
	grdFamTypes->GetSize(&width, &height);
	width= (width - grdFamTypes->GetRowLabelSize() - 20)/4;

	for (int i=0; i<4; i++) 
		grdFamTypes->SetColSize(i, width);
	
	event.Skip();
}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_TOOL1
 */

void wxDlgConfigurePedDatasets::OnBtnDelFamilyTypeClick( wxCommandEvent& event )
{
	wxArrayInt selections = grdFamTypes->GetSelectedRows();
	if (selections.GetCount() > 0) {
		gridMaster->DelFamilyType(selections[0]);
	}
}

void wxDlgConfigurePedDatasets::wxGridFamilyTypes::Refresh() {
	if (GetView()){
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 0);//, penetrances.size());
		GetView()->ProcessTableMessage(msg);
	}
}




/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_BTN_ADD_FAMILY_TYPE
 */

void wxDlgConfigurePedDatasets::OnBtnAddFamilyTypeClick( wxCommandEvent& event )
{
	gridMaster->AddFamilyType();
}

}

}



