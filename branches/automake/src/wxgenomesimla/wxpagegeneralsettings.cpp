/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagegeneralsettings.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 20 Dec 2007 15:26:39 CST
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

#include "wxpagegeneralsettings.h"
#include "growthrateconfig.h"
#include "dlggrowthreview.h"


////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxPageGeneralSettings type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageGeneralSettings, wxPanel )


/*!
 * wxPageGeneralSettings event table definition
 */

BEGIN_EVENT_TABLE( wxPageGeneralSettings, wxPanel )

////@begin wxPageGeneralSettings event table entries
    EVT_TEXT( ID_TXT_RANDOM_SEED, wxPageGeneralSettings::OnTxtRandomSeedUpdated )
    EVT_TEXT_ENTER( ID_TXT_RANDOM_SEED, wxPageGeneralSettings::OnTxtRandomSeedEnter )

    EVT_BUTTON( ID_CMD_CONFIG_GROWTH, wxPageGeneralSettings::OnCmdConfigGrowthClick )

////@end wxPageGeneralSettings event table entries
	EVT_WIZARD_PAGE_CHANGED(-1,wxPageGeneralSettings::OnWizConfirmationPageChanged)
END_EVENT_TABLE()


/*!
 * wxPageGeneralSettings constructors
 */

	wxPageGeneralSettings::wxPageGeneralSettings(): hasChanged(false)
{
    Init();
}

void wxPageGeneralSettings::OnWizConfirmationPageChanged( wxWizardEvent& event )
{
	imgGrowthChart->RefreshImage();
}

wxPageGeneralSettings::wxPageGeneralSettings( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPageGeneralSettings creator
 */

bool wxPageGeneralSettings::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageGeneralSettings creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageGeneralSettings creation
    return true;
}


/*!
 * wxPageGeneralSettings destructor
 */

wxPageGeneralSettings::~wxPageGeneralSettings()
{
////@begin wxPageGeneralSettings destruction
////@end wxPageGeneralSettings destruction
}


/*!
 * Member initialisation
 */

void wxPageGeneralSettings::Init()
{
////@begin wxPageGeneralSettings member initialisation
    txtRandomSeed = NULL;
    txtSimultaneousChrom = NULL;
    txtThreadsPerChrom = NULL;
    imgGrowthChart = NULL;
    txtGrowthParams = NULL;
////@end wxPageGeneralSettings member initialisation
}




void wxPageGeneralSettings::Init(AppController *ctl) {
	appController = ctl;
	RefreshSettings();
	hasChanged = false;
}
/**
 * @brief This is called prior to saving the configuration
 */
void wxPageGeneralSettings::Commit() {
	appController->parameters.SetRandomSeed(ExtractInteger(txtRandomSeed));
	appController->parameters.SetSimultaneousChroms(ExtractInteger(txtSimultaneousChrom));
	appController->parameters.SetThreadsPerChrom(ExtractInteger(txtThreadsPerChrom));
	hasChanged = false;
}


bool wxPageGeneralSettings::HasChanged() {
	return hasChanged;
}

/**
 * @brief This is called just after a configuration has been loaded
 */
void wxPageGeneralSettings::RefreshSettings() {
	UpdateTextField(txtRandomSeed, (int)appController->parameters.GetRandomSeed());
	UpdateTextField(txtSimultaneousChrom, (int)appController->parameters.GetSimultaneousChroms());
	UpdateTextField(txtThreadsPerChrom, (int)appController->parameters.GetThreadsPerChrom());

	
	if (imgGrowthChart) {
		
		uint estimatedGenerations = 1500;
		imgGrowthChart->SetGraphPoints( 1, estimatedGenerations, (uint)(estimatedGenerations * 0.05));
		imgGrowthChart->InitAppController(appController);


		stringstream ss;
		GrowthRate *growth = appController->parameters.GetGrowthRate();
		if (growth) {
			growth->GenerateReport( ss, 30);
			txtGrowthParams->ChangeValue(wxT(ss.str().c_str()));
		}
	}
	hasChanged = false;
}


/*!
 * Control creation for wxPageGeneralSettings
 */

void wxPageGeneralSettings::CreateControls()
{    
////@begin wxPageGeneralSettings content construction
    wxPageGeneralSettings* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer3->Add(itemBoxSizer4, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer4->Add(itemBoxSizer5, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer6 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer5->Add(itemBoxSizer6, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer6->Add(itemBoxSizer7, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticText* itemStaticText8 = new wxStaticText( itemPanel1, wxID_STATIC, _("Radom Seed:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    itemBoxSizer7->Add(itemStaticText8, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtRandomSeed = new wxTextCtrl( itemPanel1, ID_TXT_RANDOM_SEED, _("1397"), wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER, wxTextValidator(wxFILTER_NUMERIC) );
    if (wxPageGeneralSettings::ShowToolTips())
        txtRandomSeed->SetToolTip(_("The seed to be used for the pseudo random number generator"));
    itemBoxSizer7->Add(txtRandomSeed, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer5->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer11Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Threading"));
    wxStaticBoxSizer* itemStaticBoxSizer11 = new wxStaticBoxSizer(itemStaticBoxSizer11Static, wxHORIZONTAL);
    itemBoxSizer5->Add(itemStaticBoxSizer11, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer12 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer11->Add(itemBoxSizer12, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText13 = new wxStaticText( itemPanel1, wxID_STATIC, _("Simultaneous Chromosomes:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer12->Add(itemStaticText13, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtSimultaneousChrom = new wxTextCtrl( itemPanel1, ID_TXT_SIMULTANEOUS_CHROMS, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer12->Add(txtSimultaneousChrom, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer15 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer11->Add(itemBoxSizer15, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText16 = new wxStaticText( itemPanel1, wxID_STATIC, _("Threads Per Chromosome:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer15->Add(itemStaticText16, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtThreadsPerChrom = new wxTextCtrl( itemPanel1, IS_TXT_THR_PER_CHROM, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer15->Add(txtThreadsPerChrom, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer5->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer19 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer5->Add(itemBoxSizer19, 1, wxGROW|wxALL, 0);

    wxBoxSizer* itemBoxSizer20 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer19->Add(itemBoxSizer20, 1, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer21Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Population Growth"));
    wxStaticBoxSizer* itemStaticBoxSizer21 = new wxStaticBoxSizer(itemStaticBoxSizer21Static, wxHORIZONTAL);
    itemBoxSizer20->Add(itemStaticBoxSizer21, 1, wxGROW|wxALL, 5);

    imgGrowthChart = new wxImgGrowthChart( itemPanel1, ID_IMGGROWTHCHART, wxDefaultPosition, wxSize(350, 125), wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
    if (wxPageGeneralSettings::ShowToolTips())
        imgGrowthChart->SetToolTip(_("Shows the Growth Rate for current model"));
    itemStaticBoxSizer21->Add(imgGrowthChart, 2, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer23 = new wxBoxSizer(wxVERTICAL);
    itemStaticBoxSizer21->Add(itemBoxSizer23, 1, wxGROW|wxALL, 0);

    txtGrowthParams = new wxTextCtrl( itemPanel1, ID_TXT_GROWTH_CFG, _T(""), wxDefaultPosition, wxSize(200, 115), wxTE_MULTILINE|wxTE_READONLY|wxTE_LEFT|wxHSCROLL );
    if (wxPageGeneralSettings::ShowToolTips())
        txtGrowthParams->SetToolTip(_("This lists the details of the currently configured growth rate."));
    txtGrowthParams->SetFont(wxFont(10, wxTELETYPE, wxNORMAL, wxNORMAL, false, wxT("Courier")));
    itemBoxSizer23->Add(txtGrowthParams, 1, wxGROW|wxALL, 5);

    wxButton* itemButton25 = new wxButton( itemPanel1, ID_CMD_CONFIG_GROWTH, _("Configure"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer23->Add(itemButton25, 0, wxALIGN_RIGHT|wxALL, 5);

    // Connect events and objects
    imgGrowthChart->Connect(ID_IMGGROWTHCHART, wxEVT_SIZE, wxSizeEventHandler(wxPageGeneralSettings::OnSize), NULL, this);
    imgGrowthChart->Connect(ID_IMGGROWTHCHART, wxEVT_LEFT_DOWN, wxMouseEventHandler(wxPageGeneralSettings::OnLeftDown), NULL, this);
////@end wxPageGeneralSettings content construction
}


/*!
 * Should we show tooltips?
 */

bool wxPageGeneralSettings::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageGeneralSettings::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageGeneralSettings bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPageGeneralSettings bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageGeneralSettings::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageGeneralSettings icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageGeneralSettings icon retrieval
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_RANDOM_SEED
 */

void wxPageGeneralSettings::OnTxtRandomSeedUpdated( wxCommandEvent& event )
{
	hasChanged = true;
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_RANDOM_SEED in wxPageGeneralSettings.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_RANDOM_SEED in wxPageGeneralSettings. 
}


/*!
 * wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_RANDOM_SEED
 */

void wxPageGeneralSettings::OnTxtRandomSeedEnter( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_RANDOM_SEED in wxPageGeneralSettings.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_RANDOM_SEED in wxPageGeneralSettings. 
}



/*!
 * wxEVT_SIZE event handler for ID_IMGGROWTHCHART
 */

void wxPageGeneralSettings::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_IMGGROWTHCHART in wxPageGeneralSettings.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_IMGGROWTHCHART in wxPageGeneralSettings. 
}


/*!
 * wxEVT_LEFT_DOWN event handler for ID_IMGGROWTHCHART
 */

void wxPageGeneralSettings::OnLeftDown( wxMouseEvent& event )
{
////@begin wxEVT_LEFT_DOWN event handler for ID_IMGGROWTHCHART in wxPageGeneralSettings.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_LEFT_DOWN event handler for ID_IMGGROWTHCHART in wxPageGeneralSettings. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON
 */

void wxPageGeneralSettings::OnCmdConfigGrowthClick( wxCommandEvent& event )
{
	hasChanged = true;
	GrowthRate *model = appController->parameters.GetGrowthRate();
	GrowthRateType modelType = model->GetModelType();
    // Before editing this code, remove the block markers.
    GrowthCfg window(NULL, ID_DIALOG_GROWTH_RATE, _("Configure Growth Rate"));
	window.SetCurrentGrowthModel(model, modelType);
//	window.SetForceDropPoints(ChromPool::targetPop);

    if (window.ShowModal() == wxID_OK) {
		stringstream ss;
		GrowthRate *growth = window.GetGrowthModel();
		growth->GenerateReport( ss, 30);
		imgGrowthChart->SetModel( growth->GenerateCfgString().c_str() );
		txtGrowthParams->ChangeValue(wxT(ss.str().c_str()));
//		ChromPool::targetPop = window.GetForceDropPoints();
	}

////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON in wxPageGeneralSettings.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON in wxPageGeneralSettings. 
}




}

}
