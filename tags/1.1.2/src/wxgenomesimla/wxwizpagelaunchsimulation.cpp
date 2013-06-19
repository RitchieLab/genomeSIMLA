/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchsimulation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:26 CST
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

#include "wxwizpagelaunchsimulation.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxWizPageLaunchSimulation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageLaunchSimulation, wxWizardPageSimple )


/*!
 * wxWizPageLaunchSimulation event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageLaunchSimulation, wxWizardPageSimple )

////@begin wxWizPageLaunchSimulation event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageLaunchSimulation::OnWizMonitorRunPageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageLaunchSimulation::OnWizMonitorRunPageChanging )
    EVT_WIZARD_CANCEL( -1, wxWizPageLaunchSimulation::OnWizMonitorRunCancel )

////@end wxWizPageLaunchSimulation event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageLaunchSimulation constructors
 */

wxWizPageLaunchSimulation::wxWizPageLaunchSimulation()
{
    Init();
}

wxWizPageLaunchSimulation::wxWizPageLaunchSimulation( wxWizard* parent )
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageLaunchSimulation creator
 */

bool wxWizPageLaunchSimulation::Create( wxWizard* parent )
{
////@begin wxWizPageLaunchSimulation creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageLaunchSimulation creation
    return true;
}


/*!
 * wxWizPageLaunchSimulation destructor
 */

wxWizPageLaunchSimulation::~wxWizPageLaunchSimulation()
{
////@begin wxWizPageLaunchSimulation destruction
////@end wxWizPageLaunchSimulation destruction
}


/*!
 * Member initialisation
 */

void wxWizPageLaunchSimulation::Init()
{
////@begin wxWizPageLaunchSimulation member initialisation
    txtSimOutput = NULL;
    guageSimCompletion = NULL;
////@end wxWizPageLaunchSimulation member initialisation
}


/*!
 * Control creation for wxWizPageLaunchSimulation
 */

void wxWizPageLaunchSimulation::CreateControls()
{    
////@begin wxWizPageLaunchSimulation content construction
    wxWizPageLaunchSimulation* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("genomeSIMLA is now running. Progress is reported in the log window below. When the simulation has completed, you will be able to review the results. "), wxDefaultPosition, wxDefaultSize, 0 );
    itemStaticText3->Wrap(400);
    itemBoxSizer2->Add(itemStaticText3, 0, wxALIGN_LEFT|wxALL, 5);

    txtSimOutput = new wxTextCtrl( itemWizardPageSimple1, ID_TXT_SIM_LOG, _T(""), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE||wxTE_READONLY );
    itemBoxSizer2->Add(txtSimOutput, 1, wxGROW|wxALL, 5);

    guageSimCompletion = new wxGauge( itemWizardPageSimple1, ID_SIM_COMPLETION, 100, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
    guageSimCompletion->SetValue(1);
    itemBoxSizer2->Add(guageSimCompletion, 0, wxGROW|wxALL, 5);

////@end wxWizPageLaunchSimulation content construction
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_MONITOR_RUN
 */

void wxWizPageLaunchSimulation::OnWizMonitorRunPageChanged( wxWizardEvent& event )
{
	event.Skip();
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_MONITOR_RUN
 */

void wxWizPageLaunchSimulation::OnWizMonitorRunPageChanging( wxWizardEvent& event )
{
////@begin wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_MONITOR_RUN in wxWizPageLaunchSimulation.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_MONITOR_RUN in wxWizPageLaunchSimulation. 
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageLaunchSimulation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageLaunchSimulation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageLaunchSimulation bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageLaunchSimulation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageLaunchSimulation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageLaunchSimulation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageLaunchSimulation icon retrieval
}

/*!
 * wxEVT_WIZARD_CANCEL event handler for ID_WIZ_MONITOR_RUN
 */

void wxWizPageLaunchSimulation::OnWizMonitorRunCancel( wxWizardEvent& event )
{
////@begin wxEVT_WIZARD_CANCEL event handler for ID_WIZ_MONITOR_RUN in wxWizPageLaunchSimulation.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_WIZARD_CANCEL event handler for ID_WIZ_MONITOR_RUN in wxWizPageLaunchSimulation. 
}



}

}


