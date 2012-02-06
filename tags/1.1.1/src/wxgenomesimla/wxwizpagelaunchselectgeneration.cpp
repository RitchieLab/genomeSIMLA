/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchselectgeneration.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:02 CST
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

#include "wxwizpagelaunchselectgeneration.h"

////@begin XPM images
////@end XPM images
namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageLaunchSelectGeneration type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageLaunchSelectGeneration, wxWizardPageSimple )


/*!
 * wxWizPageLaunchSelectGeneration event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageLaunchSelectGeneration, wxWizardPageSimple )

////@begin wxWizPageLaunchSelectGeneration event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunPageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunPageChanging )
    EVT_WIZARD_FINISHED( ID_WIZ_SELECT_PREVIOUS_RUN, wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunFinished )
    EVT_SIZE( wxWizPageLaunchSelectGeneration::OnSize )

////@end wxWizPageLaunchSelectGeneration event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageLaunchSelectGeneration constructors
 */

wxWizPageLaunchSelectGeneration::wxWizPageLaunchSelectGeneration()
{
    Init();
}

wxWizPageLaunchSelectGeneration::wxWizPageLaunchSelectGeneration( wxWizard* parent, AppController *ctrl )
{
	appController = ctrl;
    Init();
    Create( parent );
}


/*!
 * wxWizPageLaunchSelectGeneration creator
 */

bool wxWizPageLaunchSelectGeneration::Create( wxWizard* parent )
{
////@begin wxWizPageLaunchSelectGeneration creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageLaunchSelectGeneration creation
    return true;
}


/*!
 * wxWizPageLaunchSelectGeneration destructor
 */

wxWizPageLaunchSelectGeneration::~wxWizPageLaunchSelectGeneration()
{
////@begin wxWizPageLaunchSelectGeneration destruction
////@end wxWizPageLaunchSelectGeneration destruction
}


/*!
 * Member initialisation
 */

void wxWizPageLaunchSelectGeneration::Init()
{
////@begin wxWizPageLaunchSelectGeneration member initialisation
    lblDescription = NULL;
    listGenerations = NULL;
////@end wxWizPageLaunchSelectGeneration member initialisation
}


/*!
 * Control creation for wxWizPageLaunchSelectGeneration
 */

void wxWizPageLaunchSelectGeneration::CreateControls()
{    
////@begin wxWizPageLaunchSelectGeneration content construction
    wxWizPageLaunchSelectGeneration* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxHORIZONTAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer2->Add(itemBoxSizer3, 2, wxGROW|wxALL, 5);

    lblDescription = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Starting Generation:\n\nIf this isn't the first time the application has run under this particular project name, genomeSIMLA can start over again at a previous generation. From there, the user can advance further in time using new drop configuration. \n\nSelect the generation from which you want to start. If you choose none (or 0), genomeSIMLA will start all over again. "), wxDefaultPosition, wxSize(400, 150), wxST_NO_AUTORESIZE|wxSUNKEN_BORDER );
    itemBoxSizer3->Add(lblDescription, 0, wxGROW|wxALL, 5);

    listGenerations = new wxListCtrl( itemWizardPageSimple1, ID_LST_SELECT_GENERATION, wxDefaultPosition, wxSize(100, 300), wxLC_REPORT|wxSUNKEN_BORDER );
    itemBoxSizer3->Add(listGenerations, 3, wxGROW|wxALL, 5);

////@end wxWizPageLaunchSelectGeneration content construction

	InitializeColumns();
}


void wxWizPageLaunchSelectGeneration::InitializeColumns() {
	listGenerations->ClearAll();

	listGenerations->InsertColumn(0, _("Gen."));
	listGenerations->InsertColumn(1, _("Pop. \nSize"));
	listGenerations->InsertColumn(2, _("Exp. Count"));
	listGenerations->InsertColumn(3, _("Date Completed"));
}


void wxWizPageLaunchSelectGeneration::AddEntry(ExecutionLog::LogEntry *entry) {
	int count = listGenerations->GetItemCount();
	if (entry) {
		time_t theTime = entry->endTime;
		tm *timeinfo = localtime( &theTime );
		listGenerations->InsertItem(count, _(wxString::Format("%d", (int)entry->currentGeneration)));
		listGenerations->SetItem(count, 1, _("xxxx"));
		listGenerations->SetItem(count, 2, wxString::Format("%d", (int)entry->expressionCount));
		listGenerations->SetItem(count, 3, _(asctime(timeinfo)));
		listGenerations->SetItemData(count, (int)entry->currentGeneration);
	}
	else {
		listGenerations->InsertItem(count, _("0"), 0);
		listGenerations->SetItem(count, 1, _("-"));
		listGenerations->SetItem(count, 2, _("-"));
		listGenerations->SetItem(count, 3, _("-"));
		listGenerations->SetItemData(count, 0);
	}
		
}
/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_SELECT_PREVIOUS_RUN
 */

void wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunPageChanged( wxWizardEvent& event )
{
	if (event.GetDirection()) {
		string projectName = appController->parameters.GetProjectName();
		ExecutionLog::RunType *entries	= appController->GetProjectEntries();
		InitializeColumns();

		if (entries) {
			ExecutionLog::RunType::iterator itr = entries->begin();
			ExecutionLog::RunType::iterator end = entries->end();
	
			cout<<"Preparing the Combo box: "<<entries->size()<<" drops should be there\n";
			AddEntry(NULL);
			while (itr != end) {
				ExecutionLog::LogEntry &entry = itr->second;
				AddEntry(&entry);
				itr++;
			}
			if (entries->size() > 0)
				listGenerations->FindItem(0, _("0"));
		}
	}
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_SELECT_PREVIOUS_RUN
 */
void wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunPageChanging( wxWizardEvent& event )
{
	if (event.GetDirection()) {
		int gen = 0;
	
		long selected = listGenerations->GetNextItem(-1, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
		if (selected >= 0 ){
			gen = listGenerations->GetItemData(selected);
			cout<<"Setting up generation: "<<gen<<"\n";
		}
		appController->parameters.SetStartingGeneration(gen);
	}

}

/*!
 * wxEVT_SIZE event handler for ID_WIZ_SELECT_PREVIOUS_RUN
 */

void wxWizPageLaunchSelectGeneration::OnSize( wxSizeEvent& event )
{
	if (listGenerations) {
		wxRect listSize = GetRect();
		int width = (int)(listSize.GetWidth() * 0.95);
		int perc10 = (int)(0.15 * (float)width);
		listGenerations->SetColumnWidth(0, perc10);
		listGenerations->SetColumnWidth(1, perc10);
		listGenerations->SetColumnWidth(2, perc10);
		listGenerations->SetColumnWidth(3, width - (3 * perc10));
	}
	//We need this so that the rest of the sizing stuff is performed (not by ourselves, of course)
	event.Skip();
}


/*!
 * wxEVT_WIZARD_FINISHED event handler for ID_WIZ_SELECT_PREVIOUS_RUN
 */

void wxWizPageLaunchSelectGeneration::OnWizSelectPreviousRunFinished( wxWizardEvent& event )
{
////@begin wxEVT_WIZARD_FINISHED event handler for ID_WIZ_SELECT_PREVIOUS_RUN in wxWizPageLaunchSelectGeneration.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_WIZARD_FINISHED event handler for ID_WIZ_SELECT_PREVIOUS_RUN in wxWizPageLaunchSelectGeneration. 
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageLaunchSelectGeneration::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageLaunchSelectGeneration::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageLaunchSelectGeneration bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageLaunchSelectGeneration bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageLaunchSelectGeneration::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageLaunchSelectGeneration icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageLaunchSelectGeneration icon retrieval
}


}

}


