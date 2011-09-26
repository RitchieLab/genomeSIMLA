/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpageselectgeneration.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 18 Feb 2008 13:55:24 CST
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

#include "wxwizpageselectgeneration.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxWizPageSelectGeneration type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageSelectGeneration, wxWizardPageSimple )


/*!
 * wxWizPageSelectGeneration event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageSelectGeneration, wxWizardPageSimple )

////@begin wxWizPageSelectGeneration event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageSelectGeneration::OnWizardpage1PageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageSelectGeneration::OnWizardpage1PageChanging )

////@end wxWizPageSelectGeneration event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageSelectGeneration constructors
 */

wxWizPageSelectGeneration::wxWizPageSelectGeneration()
{
    Init();
}

wxWizPageSelectGeneration::wxWizPageSelectGeneration( wxWizard* parent, AppController *ctrl )
{
	appController = ctrl;
    Init();
    Create( parent );
}


/*!
 * wxWizPageSelectGeneration creator
 */

bool wxWizPageSelectGeneration::Create( wxWizard* parent )
{
////@begin wxWizPageSelectGeneration creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageSelectGeneration creation
    return true;
}


/*!
 * wxWizPageSelectGeneration destructor
 */

wxWizPageSelectGeneration::~wxWizPageSelectGeneration()
{
////@begin wxWizPageSelectGeneration destruction
////@end wxWizPageSelectGeneration destruction
}


/*!
 * Member initialisation
 */

void wxWizPageSelectGeneration::Init()
{
////@begin wxWizPageSelectGeneration member initialisation
    treeGenerationSelection = NULL;
////@end wxWizPageSelectGeneration member initialisation
}


/*!
 * Control creation for wxWizPageSelectGeneration
 */

void wxWizPageSelectGeneration::CreateControls()
{    
////@begin wxWizPageSelectGeneration content construction
    wxWizPageSelectGeneration* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    treeGenerationSelection = new wxPageSelectGeneration( itemWizardPageSimple1, appController, ID_PANEL1, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
    itemBoxSizer2->Add(treeGenerationSelection, 1, wxGROW|wxALL, 5);

    wxCheckBox* itemCheckBox4 = new wxCheckBox( itemWizardPageSimple1, ID_SKIP_ANALYSIS, _("Skip Analysis"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    itemCheckBox4->SetValue(false);
    if (wxWizPageSelectGeneration::ShowToolTips())
        itemCheckBox4->SetToolTip(_("Skip detailed analysis and just rescan for the locus reports"));
    itemCheckBox4->SetName(_("chkSkipAnalysis"));
    itemBoxSizer2->Add(itemCheckBox4, 0, wxALIGN_RIGHT|wxALL, 5);

////@end wxWizPageSelectGeneration content construction
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZARDPAGE1
 */

void wxWizPageSelectGeneration::OnWizardpage1PageChanged( wxWizardEvent& event )
{
	appController->ClearCurrentLogEntries();
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZARDPAGE1
 */

void wxWizPageSelectGeneration::OnWizardpage1PageChanging( wxWizardEvent& event )
{
	if (GetNext()) {
		ExecutionLog::LogEntry logEntry;
		if (treeGenerationSelection->GetSelectedGeneration(logEntry))
			appController->AppendActiveEntry(logEntry);
	}
}

/*!
 * Should we show tooltips?
 */

bool wxWizPageSelectGeneration::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageSelectGeneration::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageSelectGeneration bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageSelectGeneration bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageSelectGeneration::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageSelectGeneration icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageSelectGeneration icon retrieval
}

}

}
