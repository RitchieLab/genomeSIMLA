/////////////////////////////////////////////////////////////////////////////
// Name:        wxpageselectgeneration.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 15 Feb 2008 11:08:28 CST
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

#include "wxpanelreports.h"
#include <wx/imaglist.h>
#include <sstream>
////@begin XPM images
////@end XPM images
#include "img/LD-Plot.xpm"
#include "img/RedFolderOpen.xpm"
#include "img/LD-PlotBW.xpm"
#include "img/form.xpm"
#include "img/formu.xpm"
#include "img/popc.xpm"
#include "img/popu.xpm"

namespace GenomeSIM {

namespace GUI {

/*!
 * wxPanelReports type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPanelReports, wxPanel )


/*!
 * wxPanelReports event table definition
 */

BEGIN_EVENT_TABLE( wxPanelReports, wxPanel )

////@begin wxPanelReports event table entries
    EVT_SIZE( wxPanelReports::OnSize )

	EVT_TREE_ITEM_ACTIVATED( ID_REPORT_TREE_GEN_SELECTION, wxPanelReports::OnTreeOpenReport )
//	EVT_CONTEXT_MENU(wxPanelReports::OnContextMenu)

////@end wxPanelReports event table entries
	EVT_MENU(ID_CONTINUE_SIMULATION, wxPanelReports::OnContinueSimulation)
	EVT_MENU(ID_DETAILED_ANALYSIS, wxPanelReports::OnDetailedAnalysis)
	EVT_MENU(ID_OPEN_REPORT, wxPanelReports::OpenReport)
	EVT_MENU(ID_EXTRACT_DATASETS, wxPanelReports::ExtractData)

    EVT_TREE_SEL_CHANGED( ID_REPORT_TREE_GEN_SELECTION, wxPanelReports::OnReportTreeSelChanged )
    EVT_TREE_ITEM_ACTIVATED( ID_REPORT_TREE_GEN_SELECTION, wxPanelReports::OnReportTreeItemActivated )


END_EVENT_TABLE()


void wxPanelReports::OnContinueSimulation(wxCommandEvent& event) {
	cerr<<"wxPanelReports::OnContinueSimulation()\n";
	wxMessageBox(_("wxPanelReports::OnContinueSimulation()"), _T("Continue Sim"), wxICON_WARNING|wxOK, this);
}

void wxPanelReports::OnDetailedAnalysis(wxCommandEvent& event) {
	cerr<<"wxPanelReports::OnDetailedAnalysis()\n";
	wxMessageBox(_("wxPanelReports::OnDetailedAnalysis()"), _T("Detailed Analysis"), wxICON_WARNING|wxOK, this);
}

void wxPanelReports::OpenReport(wxCommandEvent& event) {
	cerr<<"wxPAnelReports::OpenReport()\n";
	wxMessageBox(_("wxPanelReports::OpenReport()"), _T("Open Report"), wxICON_WARNING|wxOK, this);
}

void wxPanelReports::ExtractData(wxCommandEvent& event) {
	wxMessageBox(_("wxPanelReports::ExtractData()"), _T("Extract Data"), wxICON_WARNING|wxOK, this);
}
/*!
 * wxPanelReports constructors
 */

wxPanelReports::wxPanelReports() : entry(NULL)
{
    Init();
}

wxPanelReports::wxPanelReports( wxWindow* parent, AppController *ctrl, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : entry(NULL)
{
	appController = ctrl;
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPanelReports creator
 */

bool wxPanelReports::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPanelReports creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPanelReports creation
    return true;
}


/*!
 * wxPanelReports destructor
 */

wxPanelReports::~wxPanelReports()
{
////@begin wxPanelReports destruction
////@end wxPanelReports destruction
}


/*!
 * Member initialisation
 */

void wxPanelReports::Init()
{
////@begin wxPanelReports member initialisation
    treeGenerationSelection = NULL;
////@end wxPanelReports member initialisation
}


/*!
 * Control creation for wxPanelReports
 */

void wxPanelReports::CreateControls()
{    
	wxImageList *imageList = new wxImageList(16, 16);
	imageList->Add(wxIcon(popc_xpm));
	imageList->Add(wxIcon(popc_xpm));
	imageList->Add(wxIcon(popu_xpm));
	imageList->Add(wxIcon(form_xpm));
	imageList->Add(wxIcon(formu_xpm));
	imageList->Add(wxIcon(LD_PlotBW_xpm));
	imageList->Add(wxIcon(LD_Plot_xpm));
	

////@begin wxPanelReports content construction
    wxPanelReports* itemPanel1 = this;






    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);


    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer7, 0, wxALIGN_RIGHT|wxALL, 5);

    cmdOpen = new wxButton( this, ID_CMD_OPEN_REPORT, _("Open Report"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdOpen, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	cmdLaunch = new wxButton(this, ID_MENU_LAUNCH, _("Launch"), wxDefaultPosition, wxDefaultSize, 0);
	itemBoxSizer7->Add(cmdLaunch, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdContinue = new wxButton( this, ID_CMD_CONTINUE_SIMULATION, _("Continue"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdContinue, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdExtract = new wxButton( this, ID_CMD_EXTRACT_DATASETS, _("Extract Datasets"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdExtract, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdDetailed = new wxButton( this, ID_CMD_DETAILED_ANALYSIS, _("Detailed Analysis"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdDetailed, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	cmdExtractFromZero = new wxButton(this, ID_MENU_EXTRACT_DATASETS_GEN_0, _("Extract Datasets (Gen. 0)"), wxDefaultPosition, wxDefaultSize, 0);
	itemBoxSizer7->Add(cmdExtractFromZero, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
	


	cmdOpen->Enable(false);
	cmdContinue->Enable(false);
	cmdExtract->Enable(false);
	cmdDetailed->Enable(false);



//    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
 //   itemPanel1->SetSizer(itemBoxSizer2);

//    treeGenerationSelection = new wxTreeCtrl( itemPanel1, ID_TREE_GEN_SELECTION, wxDefaultPosition, wxSize(400,400), wxTR_SINGLE );
//    itemBoxSizer2->Add(treeGenerationSelection, 1, wxGROW|wxALL, 5);
	
    treeGenerationSelection = new wxTreeCtrl( itemPanel1, ID_REPORT_TREE_GEN_SELECTION, wxDefaultPosition, wxSize(100, 100), wxTR_SINGLE );
    itemBoxSizer2->Add(treeGenerationSelection, 1, wxGROW|wxALL, 5);


////@end wxPanelReports content construction
	treeGenerationSelection->AssignImageList(imageList);
}

void wxPanelReports::InitAppController(AppController *ctrl) {
	appController = ctrl;
	RefreshSettings();
}

bool wxPanelReports::HasChanged() {
	return false;
}

void wxPanelReports::Commit() {
	return;
}

/*!
 * Should we show tooltips?
 */

bool wxPanelReports::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPanelReports::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPanelReports bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPanelReports bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPanelReports::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPanelReports icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPanelReports icon retrieval
}




void wxPanelReports::AddProjectEntries(const char *project, ExecutionLog::RunType *entries) {
	wxTreeItemId pEntry = root; 
	
	if (project) {
		pEntry = treeGenerationSelection->AppendItem(root, project , 1, 1);
//		treeGenerationSelection->SetItemText(pEntry, 1, _T("-"));
//		treeGenerationSelection->SetItemText(pEntry, 2, _T("-"));
	}

	if (entries) {
		ExecutionLog::RunType::iterator cur = entries->begin();
		ExecutionLog::RunType::iterator lastEntry = entries->end();

		while (cur != lastEntry) {
			ExecutionLog::LogEntry &entry = (cur++)->second;
			
			time_t theTime = entry.endTime;
			tm *timeinfo = localtime( &theTime );
			
			LogEntryData *data = new LogEntryData(entry, "");
			stringstream ss;
			ss<<"Gen. "<<setw(8)<<left<<entry.currentGeneration<<" pop. "<<setw(8)<<left<<entry.expressionCount<<"  "<<asctime(timeinfo);
			wxTreeItemId newEntry = treeGenerationSelection->AppendItem(pEntry, wxString::Format(_T((ss.str().substr(0,ss.str().length()-1)).c_str())), 1, 2, data);
			

			map<string, string>::iterator itr = entry.reports.begin();
			map<string, string>::iterator end = entry.reports.end();

			while (itr != end) {
				string key = itr->first;
				string reportFilename = itr->second;
					LogEntryData *report = new LogEntryData(entry, reportFilename.c_str());
					int pic1 = 5, pic2 = 6;
					if (key.c_str() != "Index") {
						pic1 = 3;
						pic2 = 4;
					}
					wxTreeItemId rep = treeGenerationSelection->AppendItem(newEntry, _T(key.c_str()), pic1, pic2, report);

					treeGenerationSelection->EnsureVisible(rep);
				itr++;
			}
//				listGenerations->SetItemData(count, data);
//			treeGenerationSelection->EnsureVisible(newEntry);
		}
	}			

}

bool wxPanelReports::ValidSelection() {
	wxTreeItemId selection = treeGenerationSelection->GetSelection();
	return selection.IsOk();
}

bool wxPanelReports::GetSelectedGeneration(ExecutionLog::LogEntry &entry) {
	bool success = ValidSelection();

	if (success) {
		wxTreeItemId selection = treeGenerationSelection->GetSelection();	
		cout<<"Selected item: "<<treeGenerationSelection->GetItemText(selection)<<"\n";
		LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);
		entry = data->entry;
	}
	return success;
}

void wxPanelReports::RefreshSettings() {
	if (root.IsOk())
		treeGenerationSelection->Delete(root);
	
	root = treeGenerationSelection->AddRoot(_T("Reports"));

	vector<string> projectNames;
	appController->GetProjectList(projectNames);

	vector<string>::iterator itr = projectNames.begin();
	vector<string>::iterator end = projectNames.end();

	while (itr != end) {
		string project = *(itr++);
//	reportTree->SetItemText(gen, 1, wxString::Format(_T("%d"), (uint)entry.expressionCount));
		ExecutionLog::RunType *entries	= appController->GetProjectEntries(project.c_str());
		AddProjectEntries(project.c_str(), entries);	
	}

}

bool wxPanelReports::TransferDataFromWindow() {
	wxTreeItemId selection = treeGenerationSelection->GetSelection();
	LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);
	return data != NULL;
}

void wxPanelReports::OnTreeOpenReport( wxTreeEvent& event ){
	LogEntryData *report = (LogEntryData*)treeGenerationSelection->GetItemData(event.GetItem());
	if (report == NULL || !report->Open())
		event.Skip();
}

void wxPanelReports::OnContextMenu(wxContextMenuEvent& event) {
	wxPoint point = event.GetPosition();

	if (point.x == 0.1 && point.y == -1) {
		wxSize size = GetSize();
		point.x = size.x / 2;
		point.y = size.y / 2;
	}	else 
		point = ScreenToClient(point);

	ShowContextMenu(point);
}

void wxPanelReports::ShowContextMenu(const wxPoint& pos) {
	wxMenu menu;
		
	wxTreeItemId selection = treeGenerationSelection->GetSelection();	
	cout<<"Selected item: "<<treeGenerationSelection->GetItemText(selection)<<"\n";
	LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);
	if (data) {
		data->GetMenu(menu);
		PopupMenu(&menu, pos.x, pos.y);
	}
}

LogEntryData *wxPanelReports::GetSelectedLogEntyData() {
	wxTreeItemId selection = treeGenerationSelection->GetSelection();	
	LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);

	return data;
}


/*!
 * wxEVT_COMMAND_TREE_SEL_CHANGED event handler for ID_REPORT_TREE

	wxButton* cmdOpen;
	wxButton* cmdContinue;
	wxButton* cmdExtract;
	wxButton* cmdDetailed;

 */

void wxPanelReports::OnReportTreeSelChanged( wxTreeEvent& event )
{
	LogEntryData *data = GetSelectedLogEntyData();
	cmdContinue->Enable(data && data->CanContinue());
	cmdOpen->Enable(data && data->CanOpenReport());
	cmdDetailed->Enable(data && data->CanPerfromDetailed());
	cmdExtract->Enable(data && data->CanExtractDatasets());
//	event.Skip();
}



/*!
 * wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_REPORT_TREE
 */

void wxPanelReports::OnReportTreeItemActivated( wxTreeEvent& event )
{
	cerr<<"ReportTreeitemActivated\n";
	event.Skip();
}

}

}
