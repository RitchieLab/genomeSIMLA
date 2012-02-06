/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigureccdatasets.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 21 Mar 2008 12:07:48 CDT
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

#include "wxdlgconfigureccdatasets.h"
#include <sstream>

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxDlgConfigureCCDatasets type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgConfigureCCDatasets, wxDialog )


/*!
 * wxDlgConfigureCCDatasets event table definition
 */

BEGIN_EVENT_TABLE( wxDlgConfigureCCDatasets, wxDialog )

////@begin wxDlgConfigureCCDatasets event table entries
    EVT_BUTTON( wxID_CANCEL, wxDlgConfigureCCDatasets::OnCancelClick )

    EVT_BUTTON( wxID_OK, wxDlgConfigureCCDatasets::OnOkClick )

////@end wxDlgConfigureCCDatasets event table entries

END_EVENT_TABLE()


/*!
 * wxDlgConfigureCCDatasets constructors
 */

wxDlgConfigureCCDatasets::wxDlgConfigureCCDatasets()
{
    Init();
}

wxDlgConfigureCCDatasets::wxDlgConfigureCCDatasets( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgConfigureCCDatasets creator
 */

bool wxDlgConfigureCCDatasets::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgConfigureCCDatasets creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgConfigureCCDatasets creation
    return true;
}


/*!
 * wxDlgConfigureCCDatasets destructor
 */

wxDlgConfigureCCDatasets::~wxDlgConfigureCCDatasets()
{
////@begin wxDlgConfigureCCDatasets destruction
////@end wxDlgConfigureCCDatasets destruction
}


/*!
 * Member initialisation
 */

void wxDlgConfigureCCDatasets::Init()
{
////@begin wxDlgConfigureCCDatasets member initialisation
    txtLabel = NULL;
    txtGenotypeErr = NULL;
    txtPhenocopy = NULL;
    txtMissingData = NULL;
    txtCases = NULL;
    txtControls = NULL;
////@end wxDlgConfigureCCDatasets member initialisation
}


/*!
 * Control creation for wxDlgConfigureCCDatasets
 */

void wxDlgConfigureCCDatasets::CreateControls()
{    
////@begin wxDlgConfigureCCDatasets content construction
    wxDlgConfigureCCDatasets* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxFlexGridSizer* itemFlexGridSizer3 = new wxFlexGridSizer(3, 3, 1, 1);
    itemFlexGridSizer3->AddGrowableCol(1);
    itemFlexGridSizer3->AddGrowableCol(2);
    itemBoxSizer2->Add(itemFlexGridSizer3, 1, wxGROW|wxALL, 0);

    wxStaticText* itemStaticText4 = new wxStaticText( itemDialog1, wxID_STATIC, _("Label:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText4, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtLabel = new wxTextCtrl( itemDialog1, ID_TEXTCTRL1, _("label"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_ALPHANUMERIC) );
    itemFlexGridSizer3->Add(txtLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText10 = new wxStaticText( itemDialog1, wxID_STATIC, _("Genotype error:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText10, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtGenotypeErr = new wxTextCtrl( itemDialog1, ID_TEXTCTRL4, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtGenotypeErr, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText13 = new wxStaticText( itemDialog1, wxID_STATIC, _("Phenocopy:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText13, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtPhenocopy = new wxTextCtrl( itemDialog1, ID_TEXTCTRL5, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtPhenocopy, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    wxStaticText* itemStaticText16 = new wxStaticText( itemDialog1, wxID_STATIC, _("Missing Data:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemFlexGridSizer3->Add(itemStaticText16, 0, wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMissingData = new wxTextCtrl( itemDialog1, ID_TEXTCTRL6, _("0.0"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemFlexGridSizer3->Add(txtMissingData, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemFlexGridSizer3->Add(5, 5, 1, wxALIGN_CENTER_HORIZONTAL|wxGROW|wxALL, 5);

    itemBoxSizer2->Add(25, 25, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer20 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer20, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer21 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer20->Add(itemBoxSizer21, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText22 = new wxStaticText( itemDialog1, wxID_STATIC, _("Cases:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer21->Add(itemStaticText22, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtCases = new wxTextCtrl( itemDialog1, ID_TEXTCTRL7, _("500"), wxDefaultPosition, wxDefaultSize, 0 , wxTextValidator(wxFILTER_NUMERIC));
    itemBoxSizer21->Add(txtCases, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer24 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer20->Add(itemBoxSizer24, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText25 = new wxStaticText( itemDialog1, wxID_STATIC, _("Controls:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer24->Add(itemStaticText25, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtControls = new wxTextCtrl( itemDialog1, ID_TEXTCTRL8, _("500"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer24->Add(txtControls, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer27 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer27, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton28 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer27->Add(itemButton28, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton29 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer27->Add(itemButton29, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgConfigureCCDatasets content construction
}


/*!
 * Should we show tooltips?
 */

bool wxDlgConfigureCCDatasets::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgConfigureCCDatasets::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgConfigureCCDatasets bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgConfigureCCDatasets bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgConfigureCCDatasets::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgConfigureCCDatasets icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgConfigureCCDatasets icon retrieval
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
 */

void wxDlgConfigureCCDatasets::OnOkClick( wxCommandEvent& event )
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


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
 */

void wxDlgConfigureCCDatasets::OnCancelClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL in wxDlgConfigureCCDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL in wxDlgConfigureCCDatasets. 
}

string wxDlgConfigureCCDatasets::GetModelDetails() {
	stringstream modelDetails;
	modelDetails<<"DATASET CC "<<txtLabel->GetLineText(0)<<" "<<txtCases->GetValue()<<" "<<txtControls->GetValue()<<" "
		<<txtGenotypeErr->GetValue()<<" "<<txtPhenocopy->GetValue()<<" "<<txtMissingData->GetValue();
	return modelDetails.str();
}

string wxDlgConfigureCCDatasets::GetModelSummary() {
	stringstream summary;

	int count = ExtractInteger(txtCases) + ExtractInteger(txtControls);

	summary<<count<<" Individuals\n";
	return summary.str();
}

string wxDlgConfigureCCDatasets::GetModelLabel() {
	return txtLabel->GetLineText(0).c_str();
}
void wxDlgConfigureCCDatasets::ConfigureModelDetails(const char *details) {
	stringstream modelDetails(details);

	string val;
	//Eat the DATASET keyword
	modelDetails>>val;
	//Eat the CC keyword
	modelDetails>>val;

	modelDetails>>val;
	txtLabel->SetValue(val.c_str());
		
	modelDetails>>val;
	txtCases->SetValue(val.c_str());
		
	modelDetails>>val;
	txtControls->SetValue(val.c_str());
		
	modelDetails>>val;
	txtGenotypeErr->SetValue(val.c_str());
	
	modelDetails>>val;
	txtPhenocopy->SetValue(val.c_str());
	
	modelDetails>>val;
	txtMissingData->SetValue(val.c_str());
}


}

}

