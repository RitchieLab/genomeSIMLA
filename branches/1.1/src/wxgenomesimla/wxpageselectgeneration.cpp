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

#include "wxpageselectgeneration.h"
#include <wx/imaglist.h>
#include <sstream>
////@begin XPM images
////@end XPM images
#include "img/LD-Plot.xpm"
#include "img/RedFolderOpen.xpm"
#include "img/LD-PlotBW.xpm"
#include <wx/artprov.h>

namespace GenomeSIM {

namespace GUI {

/*!
 * wxPageSelectGeneration type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageSelectGeneration, wxPanel )


/*!
 * wxPageSelectGeneration event table definition
 */

BEGIN_EVENT_TABLE( wxPageSelectGeneration, wxPanel )

////@begin wxPageSelectGeneration event table entries
    EVT_SIZE( wxPageSelectGeneration::OnSize )

    EVT_TREE_SEL_CHANGED( ID_TREE_GEN_SELECTION, wxPageSelectGeneration::OnTreeGenSelectionSelChanged )
    EVT_TREE_ITEM_ACTIVATED( ID_TREE_GEN_SELECTION, wxPageSelectGeneration::OnTreeGenSelectionItemActivated )
    EVT_TREE_ITEM_EXPANDED( ID_TREE_GEN_SELECTION, wxPageSelectGeneration::OnTreeGenSelectionItemExpanded )

////@end wxPageSelectGeneration event table entries

END_EVENT_TABLE()


	

/*!
 * wxPageSelectGeneration constructors
 */

wxPageSelectGeneration::wxPageSelectGeneration() : entry(NULL)
{
    Init();
}

wxPageSelectGeneration::wxPageSelectGeneration( wxWindow* parent, AppController *ctrl, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : entry(NULL)
{
	appController = ctrl;
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPageSelectGeneration creator
 */

bool wxPageSelectGeneration::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageSelectGeneration creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageSelectGeneration creation
    return true;
}


/*!
 * wxPageSelectGeneration destructor
 */

wxPageSelectGeneration::~wxPageSelectGeneration()
{
////@begin wxPageSelectGeneration destruction
////@end wxPageSelectGeneration destruction
}


/*!
 * Member initialisation
 */

void wxPageSelectGeneration::Init()
{
////@begin wxPageSelectGeneration member initialisation
    treeGenerationSelection = NULL;
////@end wxPageSelectGeneration member initialisation
}


/*!
 * Control creation for wxPageSelectGeneration
 */

void wxPageSelectGeneration::CreateControls()
{    
	wxImageList *imageList = new wxImageList(16, 16);
	imageList->Add(wxIcon(wxART_HELP_FOLDER));
	imageList->Add(wxIcon(LD_PlotBW_xpm));
	imageList->Add(wxIcon(LD_Plot_xpm));

////@begin wxPageSelectGeneration content construction
    wxPageSelectGeneration* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemPanel1, wxID_STATIC, _("Select the generation from the project(s) listed below."), wxDefaultPosition, wxDefaultSize);
    itemBoxSizer2->Add(itemStaticText3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    treeGenerationSelection = new wxTreeCtrl( itemPanel1, ID_TREE_GEN_SELECTION, wxDefaultPosition, wxSize(100, 100), wxTR_HAS_BUTTONS|wxTR_HIDE_ROOT|wxTR_EXTENDED );
    itemBoxSizer2->Add(treeGenerationSelection, 1, wxGROW|wxALL, 5);

////@end wxPageSelectGeneration content construction
	treeGenerationSelection->AssignImageList(imageList);
	InitializeColumns();
}


/*!
 * Should we show tooltips?
 */

bool wxPageSelectGeneration::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageSelectGeneration::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageSelectGeneration bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPageSelectGeneration bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageSelectGeneration::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageSelectGeneration icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageSelectGeneration icon retrieval
}


/*!
 * wxEVT_SIZE event handler for ID_WXPAGESELECTGENERATION
 */

void wxPageSelectGeneration::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_WXPAGESELECTGENERATION in wxPageSelectGeneration.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_WXPAGESELECTGENERATION in wxPageSelectGeneration. 
}


/*!
 * wxEVT_COMMAND_TREE_SEL_CHANGED event handler for ID_TREE_GEN_SELECTION
 */

void wxPageSelectGeneration::OnTreeGenSelectionSelChanged( wxTreeEvent& event )
{
////@begin wxEVT_COMMAND_TREE_SEL_CHANGED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TREE_SEL_CHANGED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration. 
}


/*!
 * wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREE_GEN_SELECTION
 */

void wxPageSelectGeneration::OnTreeGenSelectionItemActivated( wxTreeEvent& event )
{
////@begin wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration. 
}


/*!
 * wxEVT_COMMAND_TREE_ITEM_EXPANDED event handler for ID_TREE_GEN_SELECTION
 */

void wxPageSelectGeneration::OnTreeGenSelectionItemExpanded( wxTreeEvent& event )
{
////@begin wxEVT_COMMAND_TREE_ITEM_EXPANDED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TREE_ITEM_EXPANDED event handler for ID_TREE_GEN_SELECTION in wxPageSelectGeneration. 
}


void wxPageSelectGeneration::AddProjectEntries(const char *project, ExecutionLog::RunType *entries) {
	wxTreeItemId pEntry = root; 
	
	if (project) {
		pEntry = treeGenerationSelection->AppendItem(root, project , 0, 0);
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
			
			LogEntryData *data = new LogEntryData(entry);
			stringstream ss;
			ss<<"Gen. "<<entry.currentGeneration<<"   \tpop. "<<entry.expressionCount<<"  \t"<<asctime(timeinfo);
			wxTreeItemId newEntry = treeGenerationSelection->AppendItem(pEntry, wxString::Format(_T((ss.str().substr(0,ss.str().length()-1)).c_str())), 1, 2, data);

//				listGenerations->SetItemData(count, data);
//			treeGenerationSelection->EnsureVisible(newEntry);
		}
	}			

}

bool wxPageSelectGeneration::ValidSelection() {
	wxTreeItemId selection = treeGenerationSelection->GetSelection();
	return selection.IsOk();
}

bool wxPageSelectGeneration::GetSelectedGeneration(ExecutionLog::LogEntry &entry) {
	bool success = ValidSelection();

	if (success) {
		wxTreeItemId selection = treeGenerationSelection->GetSelection();	
		cout<<"Selected item: "<<treeGenerationSelection->GetItemText(selection)<<"\n";
		LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);
		entry = data->entry;
	}
	return success;
}

void wxPageSelectGeneration::InitializeColumns() {
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

bool wxPageSelectGeneration::TransferDataFromWindow() {
	wxTreeItemId selection = treeGenerationSelection->GetSelection();
	LogEntryData *data = (LogEntryData*)treeGenerationSelection->GetItemData(selection);
	return data != NULL;
}


}

}
