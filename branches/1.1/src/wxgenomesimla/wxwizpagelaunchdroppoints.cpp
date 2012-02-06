/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchdroppoints.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:11 CST
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

#include "wxwizpagelaunchdroppoints.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageLaunchDropPoints type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageLaunchDropPoints, wxWizardPageSimple )


/*!
 * wxWizPageLaunchDropPoints event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageLaunchDropPoints, wxWizardPageSimple )

////@begin wxWizPageLaunchDropPoints event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageLaunchDropPoints::OnWizDropPointsPageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageLaunchDropPoints::OnWizDropPointsPageChanging )
    EVT_WIZARD_CANCEL( -1, wxWizPageLaunchDropPoints::OnWizDropPointsCancel )
    EVT_WIZARD_FINISHED( ID_WIZ_DROP_POINTS, wxWizPageLaunchDropPoints::OnWizDropPointsFinished )

    EVT_TEXT( ID_TXT_DROP_FIRST, wxWizPageLaunchDropPoints::OnTxtDropFirstUpdated )
    EVT_TEXT_ENTER( ID_TXT_DROP_FIRST, wxWizPageLaunchDropPoints::OnTxtDropFirstEnter )

    EVT_TEXT( ID_TXT_DROP_FREQ, wxWizPageLaunchDropPoints::OnTxtDropFreqUpdated )
    EVT_TEXT_ENTER( ID_TXT_DROP_FREQ, wxWizPageLaunchDropPoints::OnTxtDropFreqEnter )

    EVT_TEXT( ID_TXT_DROP_COUNT, wxWizPageLaunchDropPoints::OnTxtDropCountUpdated )
    EVT_TEXT_ENTER( ID_TXT_DROP_COUNT, wxWizPageLaunchDropPoints::OnTxtDropCountEnter )

////@end wxWizPageLaunchDropPoints event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageLaunchDropPoints constructors
 */

wxWizPageLaunchDropPoints::wxWizPageLaunchDropPoints()
{
    Init();
}

wxWizPageLaunchDropPoints::wxWizPageLaunchDropPoints( wxWizard* parent )
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageLaunchDropPoints creator
 */

bool wxWizPageLaunchDropPoints::Create( wxWizard* parent )
{
////@begin wxWizPageLaunchDropPoints creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageLaunchDropPoints creation
    return true;
}


/*!
 * wxWizPageLaunchDropPoints destructor
 */

wxWizPageLaunchDropPoints::~wxWizPageLaunchDropPoints()
{
////@begin wxWizPageLaunchDropPoints destruction
////@end wxWizPageLaunchDropPoints destruction
}


/*!
 * Member initialisation
 */

void wxWizPageLaunchDropPoints::Init()
{
////@begin wxWizPageLaunchDropPoints member initialisation
    txtDropInitial = NULL;
    txtDropFrequency = NULL;
    txtDropCount = NULL;
////@end wxWizPageLaunchDropPoints member initialisation
}


/*!
 * Control creation for wxWizPageLaunchDropPoints
 */

void wxWizPageLaunchDropPoints::CreateControls()
{    
////@begin wxWizPageLaunchDropPoints content construction
    wxWizPageLaunchDropPoints* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxHORIZONTAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    wxStaticBox* itemStaticBoxSizer3Static = new wxStaticBox(itemWizardPageSimple1, wxID_ANY, _("Pool Drop Points:"));
    wxStaticBoxSizer* itemStaticBoxSizer3 = new wxStaticBoxSizer(itemStaticBoxSizer3Static, wxVERTICAL);
    itemBoxSizer2->Add(itemStaticBoxSizer3, 1, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText4 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Describe the generations at which genomeSIMLA will perform pool \"drops\". These drops will be available for detailed analysis as well as sources for dataset productions. "), wxDefaultPosition, wxSize(400, 100), wxST_NO_AUTORESIZE );
    itemStaticBoxSizer3->Add(itemStaticText4, 0, wxGROW|wxALL, 5);

    itemStaticBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticLine* itemStaticLine6 = new wxStaticLine( itemWizardPageSimple1, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    itemStaticBoxSizer3->Add(itemStaticLine6, 0, wxGROW|wxALL, 5);

    itemStaticBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer8 = new wxBoxSizer(wxVERTICAL);
    itemStaticBoxSizer3->Add(itemBoxSizer8, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer9 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer8->Add(itemBoxSizer9, 0, wxALIGN_LEFT|wxALL, 5);

    wxStaticText* itemStaticText10 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("First Drop: This represents the first generation at which the pool is written to file. "), wxDefaultPosition, wxSize(200, -1), 0 );
    itemStaticText10->Wrap(200);
    itemBoxSizer9->Add(itemStaticText10, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticLine* itemStaticLine11 = new wxStaticLine( itemWizardPageSimple1, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_VERTICAL );
    itemBoxSizer9->Add(itemStaticLine11, 0, wxGROW|wxALL, 10);

    wxBoxSizer* itemBoxSizer12 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer9->Add(itemBoxSizer12, 1, wxALIGN_TOP|wxALL, 5);

    wxStaticText* itemStaticText13 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("First Drop:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    itemBoxSizer12->Add(itemStaticText13, 0, wxALIGN_LEFT|wxALL|wxADJUST_MINSIZE, 5);

    txtDropInitial = new wxTextCtrl( itemWizardPageSimple1, ID_TXT_DROP_FIRST, _("250"), wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER, wxTextValidator(wxFILTER_NUMERIC) );
    if (wxWizPageLaunchDropPoints::ShowToolTips())
        txtDropInitial->SetToolTip(_("The generation at which genomeSIMLA should drop and analyze the \"gene pool\". This pool can later be used for generating datasets or to perform additional analysis and/or advancements."));
    itemBoxSizer12->Add(txtDropInitial, 0, wxALIGN_LEFT|wxALL, 5);

    wxBoxSizer* itemBoxSizer15 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer8->Add(itemBoxSizer15, 0, wxALIGN_LEFT|wxALL, 5);

    wxStaticText* itemStaticText16 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Drop Frequency: Indicates how many generations between each drop will be iterated over. "), wxDefaultPosition, wxSize(200, -1), 0 );
    itemStaticText16->Wrap(200);
    itemBoxSizer15->Add(itemStaticText16, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticLine* itemStaticLine17 = new wxStaticLine( itemWizardPageSimple1, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_VERTICAL );
    itemBoxSizer15->Add(itemStaticLine17, 0, wxGROW|wxALL, 10);

    wxBoxSizer* itemBoxSizer18 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer15->Add(itemBoxSizer18, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText19 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Drop Frequency:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    itemBoxSizer18->Add(itemStaticText19, 0, wxALIGN_LEFT|wxALL|wxADJUST_MINSIZE, 5);

    txtDropFrequency = new wxTextCtrl( itemWizardPageSimple1, ID_TXT_DROP_FREQ, _("75"), wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER, wxTextValidator(wxFILTER_NUMERIC) );
    if (wxWizPageLaunchDropPoints::ShowToolTips())
        txtDropFrequency->SetToolTip(_("How many generations are to be advanced between drop points. "));
    itemBoxSizer18->Add(txtDropFrequency, 0, wxALIGN_LEFT|wxALL, 5);

    wxBoxSizer* itemBoxSizer21 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer8->Add(itemBoxSizer21, 0, wxALIGN_LEFT|wxALL, 5);

    wxStaticText* itemStaticText22 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Drop Count: Indicates how many total drops are to be made. "), wxDefaultPosition, wxSize(200, -1), 0 );
    itemStaticText22->Wrap(200);
    itemBoxSizer21->Add(itemStaticText22, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticLine* itemStaticLine23 = new wxStaticLine( itemWizardPageSimple1, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_VERTICAL );
    itemBoxSizer21->Add(itemStaticLine23, 0, wxGROW|wxALL, 10);

    wxBoxSizer* itemBoxSizer24 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer21->Add(itemBoxSizer24, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText25 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Drop Count:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
    itemBoxSizer24->Add(itemStaticText25, 0, wxALIGN_LEFT|wxALL|wxADJUST_MINSIZE, 5);

    txtDropCount = new wxTextCtrl( itemWizardPageSimple1, ID_TXT_DROP_COUNT, _("3"), wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER, wxTextValidator(wxFILTER_NUMERIC) );
    if (wxWizPageLaunchDropPoints::ShowToolTips())
        txtDropCount->SetToolTip(_("How many drops in total should be made (including initial)"));
    itemBoxSizer24->Add(txtDropCount, 0, wxALIGN_LEFT|wxALL, 5);

    itemBoxSizer8->Add(5, 5, 2, wxGROW|wxALL, 5);

////@end wxWizPageLaunchDropPoints content construction
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_DROP_POINTS
 */

void wxWizPageLaunchDropPoints::OnWizDropPointsPageChanged( wxWizardEvent& event )
{
    if (event.GetDirection()) {
		UpdateTextField(txtDropInitial, (int)(appController->parameters.GetStartingGeneration() + appController->parameters.GetDropFreq()));
		UpdateTextField(txtDropFrequency, (int)(appController->parameters.GetDropFreq()));
		UpdateTextField(txtDropCount, (int)appController->parameters.GetDropCount());
	}

}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_DROP_POINTS
 */

void wxWizPageLaunchDropPoints::OnWizDropPointsPageChanging( wxWizardEvent& event )
{

	if (event.GetDirection()) {
		appController->parameters.SetDropCount(ExtractInteger(txtDropCount));

		appController->parameters.SetDropFreq(ExtractInteger(txtDropFrequency));

		appController->parameters.SetDropInit(ExtractInteger(txtDropInitial));
	}
}


/*!
 * wxEVT_WIZARD_CANCEL event handler for ID_WIZ_DROP_POINTS
 */

void wxWizPageLaunchDropPoints::OnWizDropPointsCancel( wxWizardEvent& event )
{
////@begin wxEVT_WIZARD_CANCEL event handler for ID_WIZ_DROP_POINTS in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_WIZARD_CANCEL event handler for ID_WIZ_DROP_POINTS in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_WIZARD_FINISHED event handler for ID_WIZ_DROP_POINTS
 */

void wxWizPageLaunchDropPoints::OnWizDropPointsFinished( wxWizardEvent& event )
{
	appController->parameters.SetDropCount(ExtractInteger(txtDropCount));
	appController->parameters.SetDropFreq(ExtractInteger(txtDropFrequency));
	appController->parameters.SetDropInit(ExtractInteger(txtDropInitial));

	event.Skip();
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FIRST
 */

void wxWizPageLaunchDropPoints::OnTxtDropFirstUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FIRST in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FIRST in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FIRST
 */

void wxWizPageLaunchDropPoints::OnTxtDropFirstEnter( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FIRST in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FIRST in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FREQ
 */

void wxWizPageLaunchDropPoints::OnTxtDropFreqUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FREQ in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FREQ in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FREQ
 */

void wxWizPageLaunchDropPoints::OnTxtDropFreqEnter( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FREQ in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FREQ in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_COUNT
 */

void wxWizPageLaunchDropPoints::OnTxtDropCountUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_COUNT in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_COUNT in wxWizPageLaunchDropPoints. 
}


/*!
 * wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_COUNT
 */

void wxWizPageLaunchDropPoints::OnTxtDropCountEnter( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_COUNT in wxWizPageLaunchDropPoints.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_COUNT in wxWizPageLaunchDropPoints. 
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageLaunchDropPoints::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageLaunchDropPoints::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageLaunchDropPoints bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageLaunchDropPoints bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageLaunchDropPoints::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageLaunchDropPoints icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageLaunchDropPoints icon retrieval
}

wxWizardPage *wxWizPageLaunchDropPoints::GetPrev() const {
	wxWizardPage *prev = wxWizardPageSimple::GetPrev();
/*	if (appController->parameters.GetStartingGeneration() == 0 && prev)
		return prev->GetPrev();
	else
*/		return prev;
}

}

}
