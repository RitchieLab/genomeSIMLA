/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagereview.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 04 Feb 2008 16:58:16 CST
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

#include "wxwizpagereview.h"
#include "wxpageinfo.h"
#include <wx/html/htmlwin.h>
////@begin XPM images
////@end XPM images

#include "treelistctrl.h"

namespace GenomeSIM {

namespace GUI {

class ReportData : public wxTreeItemData {
public:
	ReportData(const char *filename) : reportFilename(filename) {}
	
	void Open() { 
		if (reportFilename.length() == 0)
			return;
		wxFileName filename(_(reportFilename.c_str()));
		filename.MakeAbsolute();

		wxString fn;
		fn=wxString(_(FILE_PREDICATE))+filename.GetFullPath();
		fn.Replace(" ", "%20", true);

		wxLaunchDefaultBrowser(fn);
			return;
	}
	string reportFilename;
};


/*!
 * wxWizPageReview type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageReview, wxWizardPageSimple )


/*!
 * wxWizPageReview event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageReview, wxWizardPageSimple )

////@begin wxWizPageReview event table entries
    EVT_TREE_ITEM_ACTIVATED( ID_TREECTRL1, wxWizPageReview::OnTreectrl1ItemActivated )

////@end wxWizPageReview event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageReview constructors
 */

wxWizPageReview::wxWizPageReview() : doLoadReports(false)	{
    Init();
}

wxWizPageReview::wxWizPageReview( wxWizard* parent ) : doLoadReports(false)	{
    Init();
    Create( parent );
}


void wxWizPageReview::LoadResults() {
		InitSimResults();
}

/*!
 * wxWizPageReivew creator
 */

bool wxWizPageReview::Create( wxWizard* parent )
{
////@begin wxWizPageReview creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageReview creation


    return true;
}


/*!
 * wxWizPageReview destructor
 */

wxWizPageReview::~wxWizPageReview()
{
////@begin wxWizPageReview destruction
////@end wxWizPageReview destruction
}


/*!
 * Member initialisation
 */

void wxWizPageReview::Init()
{
////@begin wxWizPageReview member initialisation
    reportTree = NULL;
////@end wxWizPageReview member initialisation
}


/*!
 * Control creation for wxWizPageReivew
 */

void wxWizPageReview::CreateControls()
{    
////@begin wxWizPageReview content construction
    wxWizPageReview* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxHORIZONTAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    reportTree = new wxTreeListCtrl( itemWizardPageSimple1, ID_TREECTRL1, wxDefaultPosition, wxSize(100, 100), wxTR_SINGLE );
    itemBoxSizer2->Add(reportTree, 1, wxGROW|wxALL, 5);

////@end wxWizPageReview content construction
	reportTree->AddColumn(_("Generation"));
	reportTree->SetColumnEditable(0, false);
	reportTree->AddColumn(_("Pool Size"));
	reportTree->SetColumnEditable(1, false);	
	reportTree->AddColumn(_("Report File"));
	reportTree->SetColumnEditable(2, false);
	root = reportTree->AddRoot(_("Reports"));
}

void wxWizPageReview::ClearReports() {
}	

void wxWizPageReview::AddNewEntry(ExecutionLog::LogEntry &entry) {
	wxString filename(_T(entry.reportFilename.c_str()));
	ReportData *e = new ReportData(_T(entry.reportFilename.c_str()));
	//ExecutionLog::LogEntry *e = new ExecutionLog::LogEntry(entry);
	map<string,string>::iterator itr = entry.reports.begin();
	map<string,string>::iterator end = entry.reports.end();
	wxTreeItemId gen = reportTree->AppendItem(root, wxString::Format(_T("Gen. %d"), uint(entry.currentGeneration)), -1, -1, e);
	reportTree->SetItemText(gen, 1, wxString::Format(_T("%d"), (uint)entry.expressionCount));
	reportTree->SetItemText(gen, 2, entry.reportFilename.c_str());
	while (itr != end) {
		ReportData *report = new ReportData(_T(itr->second.c_str()));
		wxTreeItemId rep = reportTree->AppendItem(gen, _T(itr->first.c_str()),-1, -1, report);
		reportTree->SetItemText(rep, 1, wxString::Format(_T("%d"), (uint)entry.expressionCount));
		//reportTree->SetItemText(gen, 2, _T(entry.reportFilename.c_str()));
		reportTree->SetItemText(rep, 2, _T(itr->second.c_str()));
		reportTree->EnsureVisible(rep);
		itr++;
	}
	reportTree->ExpandAll(root);

}

void wxWizPageReview::InitSimResults() {
	vector<ExecutionLog::LogEntry> results;
	appController->GetSessionEntries(results);
	
	ClearReports();

	size_t count = results.size();
	for (size_t i=0; i<count; i++) {
		AddNewEntry(results[i]);
	}
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageReview::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageReview::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageReview bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageReview bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageReview::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageReview icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageReview icon retrieval
}




/*!
 * wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREECTRL1
 */

void wxWizPageReview::OnTreectrl1ItemActivated( wxTreeEvent& event )
{
	ReportData *report = (ReportData*)reportTree->GetItemData(event.GetItem());
	report->Open();
}





}

}


