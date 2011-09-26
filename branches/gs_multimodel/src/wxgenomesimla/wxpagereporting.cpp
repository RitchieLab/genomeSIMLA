/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagereporting.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 21 Feb 2008 10:15:24 CST
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

#include "wxpagereporting.h"
#include "wxwizrunsimulation.h"
#include "wxdlgprocessmonitor.h"
#include "wxdlgconfigurediseasemodel.h"
////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxPageReporting type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageReporting, wxPanel )


/*!
 * wxPageReporting event table definition
 */

BEGIN_EVENT_TABLE( wxPageReporting, wxPanel )

////@begin wxPageReporting event table entries
    EVT_INIT_DIALOG( wxPageReporting::OnInitDialog )
//    EVT_SIZE( wxPageReporting::OnSize )

/*	EVT_MENU(ID_CONTINUE_SIMULATION, wxPageReporting::OnContinueSimulation)
	EVT_MENU(ID_DETAILED_ANALYSIS, wxPageReporting::OnDetailedAnalysis)
	EVT_MENU(ID_OPEN_REPORT, wxPageReporting::OpenReport)
	EVT_MENU(ID_EXTRACT_DATASETS, wxPageReporting::ExtractData)
*/
////@end wxPageReporting event table entries

END_EVENT_TABLE()


/*!
 * wxPageReporting constructors
 */

wxPageReporting::wxPageReporting()
{
    Init();
}

wxPageReporting::wxPageReporting( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPageReporting creator
 */

bool wxPageReporting::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageReporting creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageReporting creation
    return true;
}


/*!
 * wxPageReporting destructor
 */

wxPageReporting::~wxPageReporting()
{
////@begin wxPageReporting destruction
////@end wxPageReporting destruction
}


/*!
 * Member initialisation
 */

void wxPageReporting::Init()
{
////@begin wxPageReporting member initialisation
    treeGenerations = NULL;
////@end wxPageReporting member initialisation
}


/*!
 * Control creation for wxPageReporting
 */

void wxPageReporting::CreateControls()
{    
////@begin wxPageReporting content construction

////@end wxPageReporting content construction
}


/*!
 * Should we show tooltips?
 */

bool wxPageReporting::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageReporting::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageReporting bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPageReporting bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageReporting::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageReporting icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageReporting icon retrieval
}



/*!
 * wxEVT_INIT_DIALOG event handler for ID_WXPAGEREPORTING
 */

void wxPageReporting::OnInitDialog( wxInitDialogEvent& event )
{
////@begin wxEVT_INIT_DIALOG event handler for ID_WXPAGEREPORTING in wxPageReporting.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_INIT_DIALOG event handler for ID_WXPAGEREPORTING in wxPageReporting. 
}


/*!
 * wxEVT_SIZE event handler for ID_WXPAGEREPORTING
 */

void wxPageReporting::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_WXPAGEREPORTING in wxPageReporting.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_WXPAGEREPORTING in wxPageReporting. 
}

void wxPageReporting::InitAppController(AppController *ctrl) {
	appController = ctrl;   

	wxPageReporting* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);
	
//	wxTreeCtrl *theTree = new wxTreeCtrl( itemPanel1, wxID_ANY, wxDefaultPosition, wxSize(100, 100), wxTR_SINGLE );
 //   itemBoxSizer2->Add(theTree, 1, wxGROW|wxALL, 5);


    treeGenerations = new wxPanelReports( itemPanel1, ctrl, ID_PAGESELECTGENERATION, wxDefaultPosition, wxSize(100, 100), wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
	treeGenerations->InitAppController(ctrl);
    itemBoxSizer2->Add(treeGenerations, 1, wxGROW|wxALL, 5);

}



/**
 * @brief This is called prior to saving the configuration
 */
void wxPageReporting::Commit() {
	return;
}

/**
 * @brief This is called just after a configuration has been loaded
 */
void wxPageReporting::RefreshSettings() {
	cout<<"We are refreshing again!\n";
	treeGenerations->Refresh();
}

bool wxPageReporting::HasChanged() {
	return false;
}



void wxPageReporting::OnContinueSimulation(wxCommandEvent& event) {

	LogEntryData *data = treeGenerations->GetSelectedLogEntyData();
	appController->parameters.SetProjectName(data->entry.projectName.c_str());
	appController->parameters.SetStartingGeneration(data->entry.currentGeneration);

	wxWizRunSimulation *window = new wxWizRunSimulation(this, true, appController, ID_WXWIZRUNSIMULATION);
	window->Run();
	window->Destroy();
	RefreshSettings();

}
void wxPageReporting::OnDetailedAnalysis(wxCommandEvent& event) {
	LogEntryData *data = treeGenerations->GetSelectedLogEntyData();
	appController->AppendActiveEntry(data->entry);

	wxDlgProcessMonitor anl(this, appController);
	cout<<"Showing Analysis dialog\n";
	anl.ShowModal();

	Refresh();
}
void wxPageReporting::OpenReport(wxCommandEvent& event) {
	LogEntryData *data = treeGenerations->GetSelectedLogEntyData();
	if (!data->Open())
		event.Skip();

}
void wxPageReporting::ExtractData(wxCommandEvent& event) {
	LogEntryData *data = treeGenerations->GetSelectedLogEntyData();

	cerr<<"ExtractData\n";

	wxDlgConfigureDiseaseModel dlg(this, wxID_ANY, _("Model Configuration"));

	if (dlg.ShowModal()) {
		cout<<"You said OK\n";
	}

	cout<<"Extracting Datasets from generation "<<data->entry.currentGeneration<<"\n";
}

}
}

