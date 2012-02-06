/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchconfirmation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:19 CST
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

#include "wxwizpagelaunchconfirmation.h"
#include "wxbeginsimulation.h"
////@begin XPM images
////@end XPM images
#include "wxwizpagereview.h"

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageLaunchConfirmation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageLaunchConfirmation, wxWizardPageSimple )


/*!
 * wxWizPageLaunchConfirmation event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageLaunchConfirmation, wxWizardPageSimple )

////@begin wxWizPageLaunchConfirmation event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageLaunchConfirmation::OnWizConfirmationPageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageLaunchConfirmation::OnWizConfirmationPageChanging )

////@end wxWizPageLaunchConfirmation event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageLaunchConfirmation constructors
 */

wxWizPageLaunchConfirmation::wxWizPageLaunchConfirmation()
{
    Init();
}

wxWizPageLaunchConfirmation::wxWizPageLaunchConfirmation( wxWizard* parent )
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageLaunchConfirmation creator
 */

bool wxWizPageLaunchConfirmation::Create( wxWizard* parent )
{
////@begin wxWizPageLaunchConfirmation creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageLaunchConfirmation creation
    return true;
}


/*!
 * wxWizPageLaunchConfirmation destructor
 */

wxWizPageLaunchConfirmation::~wxWizPageLaunchConfirmation()
{
////@begin wxWizPageLaunchConfirmation destruction
////@end wxWizPageLaunchConfirmation destruction
}


/*!
 * Member initialisation
 */

void wxWizPageLaunchConfirmation::Init()
{
////@begin wxWizPageLaunchConfirmation member initialisation
    txtSimDetails = NULL;
////@end wxWizPageLaunchConfirmation member initialisation
}


/*!
 * Control creation for wxWizPageLaunchConfirmation
 */

void wxWizPageLaunchConfirmation::CreateControls()
{    
////@begin wxWizPageLaunchConfirmation content construction
    wxWizPageLaunchConfirmation* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Simulation Details:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer2->Add(itemStaticText3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    txtSimDetails = new wxRichTextCtrl( itemWizardPageSimple1, ID_TXT_SIM_DETAILS, _T(""), wxDefaultPosition, wxSize(100, 100), wxTE_READONLY|wxWANTS_CHARS );
    txtSimDetails->SetFont(wxFont(10, wxTELETYPE, wxNORMAL, wxNORMAL, false, wxT("Courier")));
    itemBoxSizer2->Add(txtSimDetails, 1, wxGROW|wxALL, 5);

////@end wxWizPageLaunchConfirmation content construction
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_CONFIRMATION
 */

void wxWizPageLaunchConfirmation::OnWizConfirmationPageChanged( wxWizardEvent& event )
{
	txtSimDetails->Clear();
	stringstream ss;
	appController->SummarizeSimulation(ss);
	txtSimDetails->WriteText(_(ss.str().c_str()));

}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_CONFIRMATION
 */

void wxWizPageLaunchConfirmation::OnWizConfirmationPageChanging( wxWizardEvent& event )
{	
	if (event.GetDirection()) {
		wxBeginSimulation dlg(this, wxID_ANY, _("Run Simulation"));
		dlg.Initialize(appController);
		if (dlg.ShowModal())
			((wxWizPageReview*)GetNext())->LoadResults();
	}
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageLaunchConfirmation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageLaunchConfirmation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageLaunchConfirmation bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageLaunchConfirmation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageLaunchConfirmation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageLaunchConfirmation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageLaunchConfirmation icon retrieval
}

}
}
