/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagelocusreporting.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 11 Dec 2007 09:21:06 CST
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

#include "wxpagelocusreporting.h"
#include "wxdlgeditlocusselector.h"

////@begin XPM images
#include "img/plus.xpm"
#include "img/minus.xpm"
////@end XPM images
namespace GenomeSIM {

namespace GUI {

/*!
 * wxPageLocusReporting type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageLocusReporting, wxPanel )


/*!
 * wxPageLocusReporting event table definition
 */

BEGIN_EVENT_TABLE( wxPageLocusReporting, wxPanel )

////@begin wxPageLocusReporting event table entries
    EVT_LIST_ITEM_SELECTED( ID_LST_SELECTORS, wxPageLocusReporting::OnLstSelectorsSelected )
    EVT_LIST_ITEM_DESELECTED( ID_LST_SELECTORS, wxPageLocusReporting::OnLstSelectorsDeselected )
	EVT_LIST_ITEM_ACTIVATED( ID_LST_SELECTORS, wxPageLocusReporting::OnLeftDouble )

    EVT_MENU( ID_CMD_ADD, wxPageLocusReporting::OnCmdAddClick )

    EVT_MENU( ID_CMD_REMOVE, wxPageLocusReporting::OnCmdRemoveClick )

////@end wxPageLocusReporting event table entries

END_EVENT_TABLE()


/*!
 * wxPageLocusReporting constructors
 */

wxPageLocusReporting::wxPageLocusReporting() : selectedRow(-1)	{
    Init();
}

wxPageLocusReporting::wxPageLocusReporting( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : selectedRow(-1) {
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPageLocusReporting creator
 */

bool wxPageLocusReporting::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageLocusReporting creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageLocusReporting creation
    return true;
}


/*!
 * wxPageLocusReporting destructor
 */

wxPageLocusReporting::~wxPageLocusReporting()
{
////@begin wxPageLocusReporting destruction
////@end wxPageLocusReporting destruction
}

/*
void wxPageLocusReporting::RefreshSize() {

	int width, height;
	gridRegions->GetSize(&width, &height);
	width= (int)((width - (gridRegions->GetRowLabelSize() + 20)) * 0.5);

	
}*/




/**
 * @brief This is called prior to saving the configuration
 */
void wxPageLocusReporting::Commit() {
	appController->parameters.ClearLocusSelections();
	vector<LocusSelection>::iterator itr = selections.begin();
	vector<LocusSelection>::iterator end = selections.end();

	while (itr != end) {
		appController->parameters.AddLocusSelector(*itr);
		itr++;	
	}
	
}


void wxPageLocusReporting::ClearSelections() {
	selections.clear();

	size_t count = lstLocusSelectors->GetItemCount();
	while (count) {
		lstLocusSelectors->DeleteItem(--count);
	}

}

/**
 * @brief This is called just after a configuration has been loaded
 */
void wxPageLocusReporting::RefreshSettings() {
	ClearSelections();
	LocusSelectionMap locSel = appController->parameters.GetLocusSelections();
	
	LocusSelectionMap::iterator itr = locSel.begin();
	LocusSelectionMap::iterator end = locSel.end();

	while (itr != end) {
		AddSelector(itr->second);
		itr++;
	}
}


void wxPageLocusReporting::InitList() {

	lstLocusSelectors->ClearAll();

	lstLocusSelectors->InsertColumn(0, wxT("Selector ID"));
	lstLocusSelectors->InsertColumn(1, wxT("Bl Size Criterion"));
	lstLocusSelectors->InsertColumn(2, wxT("MAF Criterion"));
	lstLocusSelectors->InsertColumn(3, wxT("Region(s)"));
	lstLocusSelectors->InsertColumn(4, wxT("Description"));

//	list->Show();
	RefreshSize();
	
}


void wxPageLocusReporting::RefreshSize() {
	int width, height;
	float widths[] = { 0.15, 0.12, 0.18, 0.2, 0.35 };
	lstLocusSelectors->GetSize(&width, &height);
//	width-=(lstLocusSelectors->GetColumnWidth(0) + 20);

	for (int i=0; i<5; i++) 
		lstLocusSelectors->SetColumnWidth(i, (uint)(widths[i] * width));

}

bool wxPageLocusReporting::HasChanged() {
	bool hasChanged = false;

	return hasChanged;
}

/*!
 * Member initialisation
 */

void wxPageLocusReporting::Init()
{
////@begin wxPageLocusReporting member initialisation
    lstLocusSelectors = NULL;
////@end wxPageLocusReporting member initialisation
}


/*!
 * Control creation for wxPageLocusReporting
 */

void wxPageLocusReporting::CreateControls()
{    
////@begin wxPageLocusReporting content construction
    wxPageLocusReporting* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    lstLocusSelectors = new wxListCtrl( itemPanel1, ID_LST_SELECTORS, wxDefaultPosition, wxSize(500, 300), wxLC_REPORT|wxLC_EDIT_LABELS|wxLC_HRULES|wxLC_VRULES );
    itemBoxSizer3->Add(lstLocusSelectors, 1, wxGROW|wxALL, 5);

    wxToolBar* itemToolBar5 = new wxToolBar( itemPanel1, wxID_ANY, wxDefaultPosition, wxSize(-1, 30), wxTB_FLAT|wxNO_BORDER );
    itemToolBar5->SetBackgroundColour(wxColour(255, 255, 255));
    itemToolBar5->SetFont(wxFont(6, wxDEFAULT, wxNORMAL, wxNORMAL, false, wxT("Sans")));
    itemToolBar5->SetToolPacking(2);
    itemToolBar5->SetToolBitmapSize(wxSize(16, 16));
    wxBitmap itemtool6Bitmap(itemPanel1->GetBitmapResource(wxT("img/plus.xpm")));
    wxBitmap itemtool6BitmapDisabled;
    itemToolBar5->AddTool(ID_CMD_ADD, _("Add Region"), itemtool6Bitmap, itemtool6BitmapDisabled, wxITEM_NORMAL, _("Add a new region to the list"), wxEmptyString);
    itemToolBar5->AddSeparator();
    wxBitmap itemtool8Bitmap(itemPanel1->GetBitmapResource(wxT("img/minus.xpm")));
    wxBitmap itemtool8BitmapDisabled;
    itemToolBar5->AddTool(ID_CMD_REMOVE, _("Remove"), itemtool8Bitmap, itemtool8BitmapDisabled, wxITEM_NORMAL, _("Delete the currently selected region from the list"), wxEmptyString);
    itemToolBar5->Realize();
    itemBoxSizer3->Add(itemToolBar5, 0, wxALIGN_RIGHT|wxALL, 0);

    // Connect events and objects
    lstLocusSelectors->Connect(ID_LST_SELECTORS, wxEVT_SIZE, wxSizeEventHandler(wxPageLocusReporting::OnSize), NULL, this);
    lstLocusSelectors->Connect(ID_LST_SELECTORS, wxEVT_LEFT_DCLICK, wxMouseEventHandler(wxPageLocusReporting::OnLeftDClick), NULL, this);
////@end wxPageLocusReporting content construction

	InitList();
}


/*!
 * Should we show tooltips?
 */

bool wxPageLocusReporting::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageLocusReporting::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageLocusReporting bitmap retrieval
    wxUnusedVar(name);
    if (name == _T("img/plus.xpm"))
    {
        wxBitmap bitmap( plus_xpm);
        return bitmap;
    }
    else if (name == _T("img/minus.xpm"))
    {
        wxBitmap bitmap( minus_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxPageLocusReporting bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageLocusReporting::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageLocusReporting icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageLocusReporting icon retrieval
}

void wxPageLocusReporting::InitAppController( AppController *ctrl) {
	appController = ctrl;

}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD
 */

void wxPageLocusReporting::OnCmdAddClick( wxCommandEvent& event )
{
	EditLocusSelector dlg(NULL, wxID_ANY, _T("Configure Locus Selection"));
	dlg.Init(chromosomes);
	if (dlg.ShowModal() == wxID_OK) {
		AddSelector(dlg.locSelection);
	}
	
////@begin wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
 //   event.Skip();
////@end wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD in wxPageLocusReporting. 
}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE
 */

void wxPageLocusReporting::OnCmdRemoveClick( wxCommandEvent& event )
{

	long item = -1;

	for (;;) {
		item = lstLocusSelectors->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
		
		if (item == -1)
			break;

		lstLocusSelectors->DeleteItem(item);
	}
/*
	wxArrayInt rows = lstLocusSelectors->GetSelectedRows();
	size_t count = rows.Count();

	for (size_t i=0; i<count; i++) 
		lstLocusSelectors->DeleteRows(rows.Item(i));
*/
////@begin wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
 //   event.Skip();
////@end wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE in wxPageLocusReporting. 
}



void wxPageLocusReporting::SetSelector(long idx, LocusSelection &sel) {
	assert(idx < (long)selections.size());

	lstLocusSelectors->SetItem(idx, 0, wxT(sel.GetLabel().c_str()));
	lstLocusSelectors->SetItem(idx, 1, wxT(sel.GetDescBlockSize().c_str()));
	lstLocusSelectors->SetItem(idx, 2, wxT(sel.GetDescMAF().c_str()));
	lstLocusSelectors->SetItem(idx, 3, wxT(sel.GetDescRanges().c_str()));
	lstLocusSelectors->SetItem(idx, 4, wxT(sel.GetDescription().c_str()));
	
}

void wxPageLocusReporting::AddSelector(LocusSelection &sel) {
	long idx = lstLocusSelectors->InsertItem(selections.size(), wxT(sel.GetLabel().c_str()));
	selections.push_back(sel);
	if (idx >= 0) 
		SetSelector(idx, sel);
}


/*!
 * wxEVT_SIZE event handler for ID_LST_SELECTORS
 */

void wxPageLocusReporting::OnSize( wxSizeEvent& event )
{
	RefreshSize();
////@begin wxEVT_SIZE event handler for ID_LST_SELECTORS in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_LST_SELECTORS in wxPageLocusReporting. 
}




void wxPageLocusReporting::OnLeftDouble( wxListEvent &event) {

	cout<<"We got a double click!\n";
	EditLocusSelector dlg(NULL, wxID_ANY, _T("Configure Locus Selection"));
	dlg.locSelection = selections[selectedRow];
	dlg.Init(chromosomes);
	if (dlg.ShowModal() == wxID_OK) {
		selections[selectedRow] = dlg.locSelection;
		SetSelector( selectedRow, dlg.locSelection);
		cout<<"Woohoo!\n";
	}

//	event.Skip();
}
/*!
 * wxEVT_LEFT_DCLICK event handler for ID_LST_SELECTORS
 */

void wxPageLocusReporting::OnLeftDClick( wxMouseEvent& event )	{

	cout<<"We got a double click!\n";
	EditLocusSelector dlg(NULL, wxID_ANY, _T("Configure Locus Selection"));
	dlg.locSelection = selections[selectedRow];
	dlg.Init(chromosomes);
	if (dlg.ShowModal() == wxID_OK) {
		selections[selectedRow] = dlg.locSelection;
		SetSelector( selectedRow, dlg.locSelection);
		cout<<"Woohoo!\n";
	}
	
////@begin wxEVT_LEFT_DCLICK event handler for ID_LST_SELECTORS in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
//    event.Skip();
////@end wxEVT_LEFT_DCLICK event handler for ID_LST_SELECTORS in wxPageLocusReporting. 
}



/*!
 * wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_SELECTORS
 */

void wxPageLocusReporting::OnLstSelectorsSelected( wxListEvent& event )
{
	selectedRow = event.GetIndex();
////@begin wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_SELECTORS in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_SELECTORS in wxPageLocusReporting. 
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_SELECTORS
 */

void wxPageLocusReporting::OnLstSelectorsDeselected( wxListEvent& event )
{
	selectedRow = -1;
////@begin wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_SELECTORS in wxPageLocusReporting.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_SELECTORS in wxPageLocusReporting. 
}

}

}




