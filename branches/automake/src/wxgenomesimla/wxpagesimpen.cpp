/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagesimpen.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 26 Feb 2008 13:23:04 CST
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

#include "wxpagesimpen.h"
#include "appcontroller.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {


/*!
 * wxPageSimPEN type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageSimPEN, wxDiseaseConfiguration )


/*!
 * wxPageSimPEN event table definition
 */

BEGIN_EVENT_TABLE( wxPageSimPEN, wxDiseaseConfiguration )

////@begin wxPageSimPEN event table entries
    EVT_CHOICE( ID_CMB_GA_SETTINGS, wxPageSimPEN::OnCmbGaSettingsSelected )

    EVT_BUTTON( ID_CMD_CONFIGURE_GA, wxPageSimPEN::OnCmdConfigureGaClick )

////@end wxPageSimPEN event table entries

END_EVENT_TABLE()


/*!
 * wxPageSimPEN constructors
 */

wxPageSimPEN::wxPageSimPEN()
{
    Init();
}

wxPageSimPEN::wxPageSimPEN( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPagesimPEN creator
 */

bool wxPageSimPEN::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageSimPEN creation
    wxDiseaseConfiguration::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageSimPEN creation
    return true;
}


/*!
 * wxPageSimPEN destructor
 */

wxPageSimPEN::~wxPageSimPEN()
{
////@begin wxPageSimPEN destruction
////@end wxPageSimPEN destruction
}


/*!
 * Member initialisation
 */

void wxPageSimPEN::Init()
{
////@begin wxPageSimPEN member initialisation
    txtPrevalence = NULL;
    txtHeritTarget = NULL;
    txtHeritWeight = NULL;
    txtORTarget = NULL;
    txtORWeight = NULL;
    txtVarTarget = NULL;
    txtVarianceWeight = NULL;
    cmbGASettings = NULL;
////@end wxPageSimPEN member initialisation
}


/*!
 * Control creation for wxPagesimPEN
 */

void wxPageSimPEN::CreateControls()
{    
////@begin wxPageSimPEN content construction
    wxPageSimPEN* itemDiseaseConfiguration1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDiseaseConfiguration1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 0, wxALIGN_LEFT|wxALL, 5);

    wxStaticText* itemStaticText4 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Target Prevalence:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer3->Add(itemStaticText4, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtPrevalence = new wxTextCtrl( itemDiseaseConfiguration1, ID_TEXTCTRL, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer3->Add(txtPrevalence, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer6 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer6, 0, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer6->Add(itemBoxSizer7, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer8Static = new wxStaticBox(itemDiseaseConfiguration1, wxID_ANY, _("Heritability"));
    wxStaticBoxSizer* itemStaticBoxSizer8 = new wxStaticBoxSizer(itemStaticBoxSizer8Static, wxVERTICAL);
    itemBoxSizer7->Add(itemStaticBoxSizer8, 0, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer9 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer8->Add(itemBoxSizer9, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText10 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Target:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer9->Add(itemStaticText10, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer9->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtHeritTarget = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_HERIT_TARGET, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer9->Add(txtHeritTarget, 0, wxGROW|wxALL, 0);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer8->Add(itemBoxSizer13, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText14 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Weight:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemStaticText14, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer13->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtHeritWeight = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_HERIT_WEIGHT, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer13->Add(txtHeritWeight, 0, wxGROW|wxALL, 0);

    wxStaticBox* itemStaticBoxSizer17Static = new wxStaticBox(itemDiseaseConfiguration1, wxID_ANY, _("Odds Ratio"));
    wxStaticBoxSizer* itemStaticBoxSizer17 = new wxStaticBoxSizer(itemStaticBoxSizer17Static, wxVERTICAL);
    itemBoxSizer7->Add(itemStaticBoxSizer17, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer18 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer17->Add(itemBoxSizer18, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticText* itemStaticText19 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Target:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer18->Add(itemStaticText19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer18->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtORTarget = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_OR_TARGET, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer18->Add(txtORTarget, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer22 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer17->Add(itemBoxSizer22, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticText* itemStaticText23 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Weight:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer22->Add(itemStaticText23, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer22->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtORWeight = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_OR_WEIGHT, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer22->Add(txtORWeight, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxStaticBox* itemStaticBoxSizer26Static = new wxStaticBox(itemDiseaseConfiguration1, wxID_ANY, _("Marginal Variance"));
    wxStaticBoxSizer* itemStaticBoxSizer26 = new wxStaticBoxSizer(itemStaticBoxSizer26Static, wxVERTICAL);
    itemBoxSizer7->Add(itemStaticBoxSizer26, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer27 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer26->Add(itemBoxSizer27, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText28 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Target:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer27->Add(itemStaticText28, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer27->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtVarTarget = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_VARIANCE_TARGET, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer27->Add(txtVarTarget, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer31 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer26->Add(itemBoxSizer31, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText32 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("Weight:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer31->Add(itemStaticText32, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer31->Add(5, 5, 1, wxGROW|wxALL, 5);

    txtVarianceWeight = new wxTextCtrl( itemDiseaseConfiguration1, ID_TXT_VARIANCE_WEIGHT, _T(""), wxDefaultPosition, wxSize(125, -1), 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer31->Add(txtVarianceWeight, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer35 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer35, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText36 = new wxStaticText( itemDiseaseConfiguration1, wxID_STATIC, _("GA Settings:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer35->Add(itemStaticText36, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxArrayString cmbGASettingsStrings;
    cmbGASettings = new wxChoice( itemDiseaseConfiguration1, ID_CMB_GA_SETTINGS, wxDefaultPosition, wxDefaultSize, cmbGASettingsStrings, 0 );
    itemBoxSizer35->Add(cmbGASettings, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton38 = new wxButton( itemDiseaseConfiguration1, ID_CMD_CONFIGURE_GA, _("Configure"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer35->Add(itemButton38, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxPageSimPEN content construction
	int settingCount = AppController::gaSettings.size();
	while (settingCount-- > 0) 
		cmbGASettings->Insert(_(AppController::gaSettings[settingCount].c_str()), settingCount);
//	cmbGASettings->Insert(_T("tight.ga"), 0);
	cmbGASettings->SetSelection(0);
}


/*!
 * wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_GA_SETTINGS
 */

void wxPageSimPEN::OnCmbGaSettingsSelected( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_GA_SETTINGS in wxPagesimPEN.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_GA_SETTINGS in wxPagesimPEN. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIGURE_GA
 */

void wxPageSimPEN::OnCmdConfigureGaClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIGURE_GA in wxPagesimPEN.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIGURE_GA in wxPagesimPEN. 
}


/*!
 * Should we show tooltips?
 */

bool wxPageSimPEN::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageSimPEN::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageSimPEN bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPageSimPEN bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageSimPEN::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageSimPEN icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageSimPEN icon retrieval
}

bool wxPageSimPEN::Save(const char *filename, const char *description) {
	bool success = false;
	int selection = cmbGASettings->GetSelection();

	//Possibly verify the settings
	ofstream file(filename);
	file<<"HERITWEIGHT\t"<<ExtractDouble(txtHeritWeight)<<"\n";
	file<<"HERIT\t"<<ExtractDouble(txtHeritTarget)<<"\n";
	file<<"MARGWEIGHT\t"<<ExtractDouble(txtVarianceWeight)<<"\n";
	file<<"MARGVAR\t"<<ExtractDouble(txtVarTarget)<<"\n";
	file<<"ODDSWEIGHT\t"<<ExtractDouble(txtORWeight)<<"\n";
	file<<"ODDSRATIO\t"<<ExtractDouble(txtORTarget)<<"\n";
	file<<"PENTARGET\t"<<ExtractDouble(txtPrevalence)<<"\n";
	file<<"\n#GA Settings (filename containing the GA details\n";
	file<<"GA_SETTINGS "<<cmbGASettings->GetString(selection).c_str()<<"\n";


	return true;
	
}
bool wxPageSimPEN::Import(const char *filename, string &desc) {
	bool success = false;
	int selection = cmbGASettings->GetSelection();
	if (selection >= 0) {
		success = fileContents.Parse(cmbGASettings->GetString(selection).c_str());

		success = success && fileContents.Parse(filename);
		if (success) {
			string comment = fileContents.comments;
	
			cout<<"Original Comment: "<<comment<<"\n";
			if (comment.length() > 0) {
				size_t loc = comment.rfind("\n#", comment.length() - 1);
					
				while (loc != string::npos) {
					cout<<"Found at: "<<loc<<"\n";
					comment = comment.erase(loc + 1, 1);
					loc = comment.rfind("\n#", loc - 1);	
				}
			}		
	
			if (comment[0] == '#')
				comment=comment.erase(0,1);
			cout<<"Loaded Comment:   "<<comment<<"\n";
	
			string value;
			txtHeritWeight->SetValue(fileContents.GetLine("HERITWEIGHT").c_str());		
			txtHeritTarget->SetValue(fileContents.GetLine("HERIT").c_str());
			txtVarianceWeight->SetValue(fileContents.GetLine("MARGWEIGHT").c_str());
			txtVarTarget->SetValue(fileContents.GetLine("MARGVAR").c_str());
			txtORWeight->SetValue(fileContents.GetLine("ODDSWEIGHT").c_str());
			txtORTarget->SetValue(fileContents.GetLine("ODDSRATIO").c_str());
			txtPrevalence->SetValue(fileContents.GetLine("PENTARGET").c_str());	
			
		}
	}

}



}

}
