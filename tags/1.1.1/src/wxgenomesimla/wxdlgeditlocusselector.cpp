/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgeditlocusselector.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 13 Dec 2007 14:03:55 CST
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

#include "wxdlgeditlocusselector.h"
#include "wxgridloci.h"
#include "wxdlgselectlocusrange.h"

////@begin XPM images
#include "img/plus.xpm"
#include "img/minus.xpm"
////@end XPM images

using namespace Simulation;

namespace GenomeSIM {

namespace GUI {

/*!
 * EditLocusSelector type definition
 */

IMPLEMENT_DYNAMIC_CLASS( EditLocusSelector, wxDialog )

	
/*!
 * EditLocusSelector event table definition
 */

BEGIN_EVENT_TABLE( EditLocusSelector, wxDialog )
	EVT_COMBOBOX(wxID_ANY, EditLocusSelector::OnCombo)
////@begin EditLocusSelector event table entries
    EVT_SIZE( EditLocusSelector::OnSize )

    EVT_TEXT( ID_TXT_LABEL, EditLocusSelector::OnTxtLabelTextUpdated )

    EVT_TEXT( ID_TXT_DESCRIPTION, EditLocusSelector::OnTxtDescriptionTextUpdated )

    EVT_TEXT( ID_TXT_BL_TARGET, EditLocusSelector::OnTxtBlTargetTextUpdated )

    EVT_TEXT( ID_TXT_BL_MIN, EditLocusSelector::OnTxtBlMinTextUpdated )

    EVT_TEXT( ID_TXT_BL_MAX, EditLocusSelector::OnTxtBlMaxTextUpdated )

    EVT_TEXT( ID_TXT_MAC_TARGET, EditLocusSelector::OnTxtMacTargetTextUpdated )

    EVT_TEXT( ID_TXT_MAF_MIN, EditLocusSelector::OnTxtMafMinTextUpdated )

    EVT_TEXT( ID_TXT_MAF_MAX, EditLocusSelector::OnTxtMafMaxTextUpdated )

    EVT_GRID_CELL_LEFT_DCLICK( EditLocusSelector::OnLeftDClick )
    EVT_GRID_CELL_CHANGE( EditLocusSelector::OnCellChange )
    EVT_GRID_EDITOR_HIDDEN( EditLocusSelector::OnEditorHidden )
    EVT_GRID_EDITOR_SHOWN( EditLocusSelector::OnEditorShown )

    EVT_MENU( ID_CMD_ADD_REGION, EditLocusSelector::OnCmdAddRegionClick )

    EVT_MENU( ID_CMD_REMOVE_REGION, EditLocusSelector::OnCmdRemoveRegionClick )

    EVT_BUTTON( wxID_OK, EditLocusSelector::OnOkClick )

////@end EditLocusSelector event table entries

END_EVENT_TABLE()


/*!
 * EditLocusSelector constructors
 */

EditLocusSelector::EditLocusSelector() : activeRow(-1), initializationComplete(false)
{
    Init();
}

EditLocusSelector::EditLocusSelector( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ): activeRow(-1), initializationComplete(false)
{
	
    Init();
    Create(parent, id, caption, pos, size, style);
}

void EditLocusSelector::Init(vector<FileBasedChromosome *> *chrom) {

	vector<FileBasedChromosome *>::iterator itr = chrom->begin();
	vector<FileBasedChromosome *>::iterator end = chrom->end();

	while (itr != end) {
		FileBasedChromosome *cur = *itr++;
		chromosomes[cur->label] = cur;
	}
		
	initializationComplete = false;

	Init();
	CreateControls();

	wxString val;

	LinearRange<float> &mafRange = locSelection.GetMafRange();
	val.Printf("%f", mafRange.target);
	txtMafTarget->SetValue(val);
	
	val.Printf("%f", mafRange.min);
	txtMafMin->SetValue(val);

	val.Printf("%f", mafRange.max);
	txtMafMax->SetValue(val);

	LinearRange<uint> &blRange = locSelection.GetBlockSizeRange();;
	val.Printf("%d", blRange.target);
	txtBlSizeTarget->SetValue(val);

	val.Printf("%d", blRange.min);
	txtBlSizeMin->SetValue(val);

	val.Printf("%d", blRange.max);
	txtBlSizeMax->SetValue(val);
	
	size_t count = locSelection.GetRangeCount();
	if (gridRegions->InsertRows(0, count)) {
		for (size_t i=0; i<count; i++) {
			SimpleRange<Locus> &range = locSelection.GetRange(i);
			RenderRegion(i, &(range.min), &(range.max));
		}
	}
	

    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();

	initializationComplete = true;
}


/*!
 * EditLocusSelector creator
 */

bool EditLocusSelector::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );
    return true;
}


/*!
 * EditLocusSelector destructor
 */

EditLocusSelector::~EditLocusSelector()
{
	delete[] labels;
////@begin EditLocusSelector destruction
////@end EditLocusSelector destruction
}


/*!
 * Member initialisation
 */

void EditLocusSelector::Init()
{
////@begin EditLocusSelector member initialisation
    txtLabel = NULL;
    txtDescription = NULL;
    txtBlSizeTarget = NULL;
    txtBlSizeMin = NULL;
    txtBlSizeMax = NULL;
    txtMafTarget = NULL;
    txtMafMin = NULL;
    txtMafMax = NULL;
    gridRegions = NULL;
////@end EditLocusSelector member initialisation
}


/*!
 * Control creation for EditLocusSelector
 */

void EditLocusSelector::CreateControls()
{    
////@begin EditLocusSelector content construction
    EditLocusSelector* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer4, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer4->Add(itemBoxSizer5, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText6 = new wxStaticText( itemDialog1, wxID_STATIC, _("Selection Label:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer5->Add(itemStaticText6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtLabel = new wxTextCtrl( itemDialog1, ID_TXT_LABEL, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_ALPHANUMERIC) );
    itemBoxSizer5->Add(txtLabel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer4->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 10);

    wxBoxSizer* itemBoxSizer9 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer4->Add(itemBoxSizer9, 0, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText10 = new wxStaticText( itemDialog1, wxID_STATIC, _("Description:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer9->Add(itemStaticText10, 0, wxALIGN_LEFT|wxALL, 5);

    txtDescription = new wxTextCtrl( itemDialog1, ID_TXT_DESCRIPTION, _T(""), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer9->Add(txtDescription, 1, wxGROW|wxALL, 5);

    itemBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 10);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer13, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer14Static = new wxStaticBox(itemDialog1, wxID_ANY, _("Block Size Criterion"));
    wxStaticBoxSizer* itemStaticBoxSizer14 = new wxStaticBoxSizer(itemStaticBoxSizer14Static, wxHORIZONTAL);
    itemBoxSizer13->Add(itemStaticBoxSizer14, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer15 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer14->Add(itemBoxSizer15, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText16 = new wxStaticText( itemDialog1, wxID_STATIC, _("Target:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer15->Add(itemStaticText16, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtBlSizeTarget = new wxTextCtrl( itemDialog1, ID_TXT_BL_TARGET, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer15->Add(txtBlSizeTarget, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer18 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer14->Add(itemBoxSizer18, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText19 = new wxStaticText( itemDialog1, wxID_STATIC, _("Min:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer18->Add(itemStaticText19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtBlSizeMin = new wxTextCtrl( itemDialog1, ID_TXT_BL_MIN, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer18->Add(txtBlSizeMin, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer21 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer14->Add(itemBoxSizer21, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText22 = new wxStaticText( itemDialog1, wxID_STATIC, _("Max:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer21->Add(itemStaticText22, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtBlSizeMax = new wxTextCtrl( itemDialog1, ID_TXT_BL_MAX, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer21->Add(txtBlSizeMax, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer13->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer25Static = new wxStaticBox(itemDialog1, wxID_ANY, _("MAF Criterion"));
    wxStaticBoxSizer* itemStaticBoxSizer25 = new wxStaticBoxSizer(itemStaticBoxSizer25Static, wxHORIZONTAL);
    itemBoxSizer13->Add(itemStaticBoxSizer25, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer26 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer25->Add(itemBoxSizer26, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText27 = new wxStaticText( itemDialog1, wxID_STATIC, _("Target:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer26->Add(itemStaticText27, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMafTarget = new wxTextCtrl( itemDialog1, ID_TXT_MAC_TARGET, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer26->Add(txtMafTarget, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer29 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer25->Add(itemBoxSizer29, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText30 = new wxStaticText( itemDialog1, wxID_STATIC, _("Min:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer29->Add(itemStaticText30, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMafMin = new wxTextCtrl( itemDialog1, ID_TXT_MAF_MIN, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer29->Add(txtMafMin, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer32 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer25->Add(itemBoxSizer32, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText33 = new wxStaticText( itemDialog1, wxID_STATIC, _("Max:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer32->Add(itemStaticText33, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMafMax = new wxTextCtrl( itemDialog1, ID_TXT_MAF_MAX, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer32->Add(txtMafMax, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 10);

    wxStaticText* itemStaticText36 = new wxStaticText( itemDialog1, wxID_STATIC, _("Region:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer3->Add(itemStaticText36, 0, wxALIGN_RIGHT|wxALL, 0);

    gridRegions = new wxGrid( itemDialog1, ID_GRD_REGIONS, wxDefaultPosition, wxSize(500, 300), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    gridRegions->SetDefaultColSize(50);
    gridRegions->SetDefaultRowSize(25);
    gridRegions->SetColLabelSize(25);
    gridRegions->SetRowLabelSize(50);
    itemBoxSizer3->Add(gridRegions, 1, wxGROW|wxALL, 0);

    wxToolBar* itemToolBar38 = new wxToolBar( itemDialog1, wxID_ANY, wxDefaultPosition, wxSize(-1, 30), wxTB_FLAT|wxNO_BORDER );
    itemToolBar38->SetBackgroundColour(wxColour(255, 255, 255));
    itemToolBar38->SetFont(wxFont(6, wxDEFAULT, wxNORMAL, wxNORMAL, false, wxT("Sans")));
    itemToolBar38->SetToolPacking(2);
    itemToolBar38->SetToolBitmapSize(wxSize(16, 16));
    wxBitmap itemtool39Bitmap(itemDialog1->GetBitmapResource(wxT("img/plus.xpm")));
    wxBitmap itemtool39BitmapDisabled;
    itemToolBar38->AddTool(ID_CMD_ADD_REGION, _("Add Region"), itemtool39Bitmap, itemtool39BitmapDisabled, wxITEM_NORMAL, _("Add a new region to the list"), wxEmptyString);
    itemToolBar38->AddSeparator();
    wxBitmap itemtool41Bitmap(itemDialog1->GetBitmapResource(wxT("img/minus.xpm")));
    wxBitmap itemtool41BitmapDisabled;
    itemToolBar38->AddTool(ID_CMD_REMOVE_REGION, _("Remove"), itemtool41Bitmap, itemtool41BitmapDisabled, wxITEM_NORMAL, _("Delete the currently selected region from the list"), wxEmptyString);
    itemToolBar38->Realize();
    itemBoxSizer3->Add(itemToolBar38, 0, wxALIGN_RIGHT|wxALL, 0);

    itemBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer43 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer3->Add(itemBoxSizer43, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton44 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer43->Add(itemButton44, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton45 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer43->Add(itemButton45, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // Connect events and objects
    gridRegions->Connect(ID_GRD_REGIONS, wxEVT_SIZE, wxSizeEventHandler(EditLocusSelector::OnSize), NULL, this);
////@end EditLocusSelector content construction
	InitList();

	txtLabel->SetValue((locSelection.GetLabel().c_str()));
	txtDescription->SetValue((locSelection.GetDescription().c_str()));



}


void EditLocusSelector::RenderRegion(int pos, Locus *start, Locus *stop) {
	gridRegions->SetCellValue(pos, 0, wxT(labels[start->GetChromID()]));
	gridRegions->SetCellValue(pos, 1, wxT(start->GetLabel().c_str()));	
	gridRegions->SetReadOnly(pos, 1);
	gridRegions->SetCellValue(pos, 2, wxT(stop->GetLabel().c_str()));
	gridRegions->SetReadOnly(pos, 2);
}

void EditLocusSelector::AddNewRange() {
	//We'll start with the first chromosome, and the first/last chromosome in it 
	FileBasedChromosome *chr = chromosomes.begin()->second;
	Locus *l1 = chr->GetLocus(0);
	Locus *l2 = chr->GetLocus(chr->GetLocusCount() - 1);
	gridRegions->InsertRows(locSelection.GetRangeCount(), 1);
	RenderRegion(locSelection.GetRangeCount(), l1, l2);
	locSelection.AddRange( *l1, *l2);
}


void EditLocusSelector::OnCombo(wxCommandEvent &event) {
	wxString selection = event.GetString();
	cout<<"This was the string value: "<<event.GetString()<<"\n";
	if (activeRow >= 0) {
		
		FileBasedChromosome *chr = chromosomes[selection.c_str()];
		Locus *l1 = chr->GetLocus(0);
		Locus *l2 = chr->GetLocus(chr->GetLocusCount() - 1);
		RenderRegion(activeRow, l1, l2);
		locSelection.SetRange(activeRow, *l1, *l2);
	}
}

/*!
 * wxEVT_SIZE event handler for ID_EDITLOCUSSELECTOR
 */

void EditLocusSelector::OnSize( wxSizeEvent& event )	{
	RefreshSize();
////@begin wxEVT_SIZE event handler for ID_EDITLOCUSSELECTOR in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_EDITLOCUSSELECTOR in EditLocusSelector. 
}


/*!
 * wxEVT_GRID_CELL_LEFT_DCLICK event handler for ID_GRD_REGIONS
 */

void EditLocusSelector::OnLeftDClick( wxGridEvent& event )		{
	int row = event.GetRow();
	int col = event.GetCol() - 1;
	wxString selection = gridRegions->GetCellValue(row, 0);

	if (row >= 0 && col >= 0 && col < 2) {
		Locus *loc[2];
		FileBasedChromosome *chr = chromosomes[selection.c_str()];

		loc[0] = chr->GetLocus(gridRegions->GetCellValue(row, 1).c_str());
		loc[1] = chr->GetLocus(gridRegions->GetCellValue(row, 2).c_str());
		
		wxDlgSelectLocusRange dlg(NULL, wxID_ANY);
		dlg.Initialize(chr, loc[col]);
		
		if (dlg.ShowModal() == wxID_OK) {
			loc[col] = dlg.GetSelection();
			if (loc[col]) {
				gridRegions->SetCellValue(row, col + 1, loc[col]->GetLabel().c_str());
				locSelection.SetRange(row, *(loc[0]), *(loc[1]));				
			}
		}

	}	
	
}


/*!
 * wxEVT_GRID_CELL_CHANGE event handler for ID_GRD_REGIONS
 */

void EditLocusSelector::OnCellChange( wxGridEvent& event )
{
	cout<<"OnCellChange()\n";
////@begin wxEVT_GRID_CELL_CHANGE event handler for ID_GRD_REGIONS in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_CELL_CHANGE event handler for ID_GRD_REGIONS in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD_REGION
 */

void EditLocusSelector::OnCmdAddRegionClick( wxCommandEvent& event )
{

	if (chromosomes.size() > 0) {
		AddNewRange();
	}
	else {
		wxMessageDialog dlg(NULL, wxT("There are no file based chromosomes to use. "), wxT("No File Based Chromosomes"), wxOK);
		dlg.ShowModal();
	}
////@begin wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD_REGION in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD_REGION in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE_REGION
 */

void EditLocusSelector::OnCmdRemoveRegionClick( wxCommandEvent& event )
{
	wxArrayInt rows = gridRegions->GetSelectedRows();
	size_t count = rows.Count();
	
	cout<<"Deleting: "<<count<<" columns\n";
	for (size_t i=0; i<count; i++) {
		locSelection.RemoveRange( rows.Item(i) );
		gridRegions->DeleteRows(rows.Item(i));
	}
	
	
////@begin wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE_REGION in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE_REGION in EditLocusSelector. 
}



void EditLocusSelector::RefreshSize() {

	float widths[] = { 0.4, 0.3, 0.3 };
	int width, height;
	gridRegions->GetSize(&width, &height);
	width-= (gridRegions->GetRowLabelSize() + 20);

	gridRegions->SetColSize(0, (int)(width * widths[0]));
	gridRegions->SetColSize(1, (int)(width * widths[1]));
	gridRegions->SetColSize(2, (int)(width * widths[2]));
	
}


void EditLocusSelector::InitList() {
	uint chromCount = chromosomes.size();
	map<string, FileBasedChromosome *>::iterator itr = chromosomes.begin();
	map<string, FileBasedChromosome *>::iterator end = chromosomes.end();
	labels = new wxString[chromCount];

	int i=0;
	while (itr != end) {
		labels[i++] = wxT(itr->second->label.c_str());
		itr++;
	}

	gridRegions->ClearGrid();
	gridRegions->CreateGrid(0, 3);

	gridRegions->SetRowLabelSize(20);
	//gridRegions->SetColLabelSize(40);

	gridRegions->SetColLabelValue(0, wxT("Chromosome"));
	gridRegions->SetColLabelValue(1, wxT("Start"));
	gridRegions->SetColLabelValue(2, wxT("Stop"));

	wxGridCellAttr *attrCombo  = new wxGridCellAttr;
	
	attrCombo->SetEditor(new wxGridCellChoiceEditor(chromCount,
                                                    labels));
	gridRegions->SetColAttr(0, attrCombo);

	RefreshSize();
	
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
 */

void EditLocusSelector::OnOkClick( wxCommandEvent& event )
{
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}


/*!
 * Should we show tooltips?
 */

bool EditLocusSelector::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap EditLocusSelector::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin EditLocusSelector bitmap retrieval
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
////@end EditLocusSelector bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon EditLocusSelector::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin EditLocusSelector icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end EditLocusSelector icon retrieval
}



/*!
 * wxEVT_GRID_EDITOR_HIDDEN event handler for ID_GRD_REGIONS
 */

void EditLocusSelector::OnEditorHidden( wxGridEvent& event )
{	
	activeRow = -1;
////@begin wxEVT_GRID_EDITOR_HIDDEN event handler for ID_GRD_REGIONS in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_EDITOR_HIDDEN event handler for ID_GRD_REGIONS in EditLocusSelector. 
}


/*!
 * wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_REGIONS
 */

void EditLocusSelector::OnEditorShown( wxGridEvent& event )
{
	activeRow = event.GetRow();
////@begin wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_REGIONS in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_REGIONS in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_LABEL
 */

void EditLocusSelector::OnTxtLabelTextUpdated( wxCommandEvent& event )
{

	if (txtLabel)
		locSelection.SetLabel( txtLabel->GetLineText(0).c_str());
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_LABEL in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_LABEL in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DESCRIPTION
 */

void EditLocusSelector::OnTxtDescriptionTextUpdated( wxCommandEvent& event )
{
	if (txtDescription) 
		locSelection.SetDescription(txtDescription->GetLineText(0).c_str());
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DESCRIPTION in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DESCRIPTION in EditLocusSelector. 
}



/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_TARGET
 */

void EditLocusSelector::OnTxtBlTargetTextUpdated( wxCommandEvent& event )
{
	UpdateBlockSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_TARGET in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_TARGET in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MIN
 */

void EditLocusSelector::OnTxtBlMinTextUpdated( wxCommandEvent& event )
{
	UpdateBlockSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MIN in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MIN in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MAX
 */

void EditLocusSelector::OnTxtBlMaxTextUpdated( wxCommandEvent& event )
{
	UpdateBlockSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MAX in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MAX in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAC_TARGET
 */

void EditLocusSelector::OnTxtMacTargetTextUpdated( wxCommandEvent& event )
{
	UpdateMAFSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAC_TARGET in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAC_TARGET in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MIN
 */

void EditLocusSelector::OnTxtMafMinTextUpdated( wxCommandEvent& event )
{
	UpdateMAFSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MIN in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MIN in EditLocusSelector. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MAX
 */

void EditLocusSelector::OnTxtMafMaxTextUpdated( wxCommandEvent& event )
{
	UpdateMAFSelection();	
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MAX in EditLocusSelector.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MAX in EditLocusSelector. 
}

void EditLocusSelector::UpdateBlockSelection() {
	if (txtBlSizeTarget && txtBlSizeMin && txtBlSizeMax) {
		if (initializationComplete) {
			LinearRange<uint> range(0,0,0);
			long value;
			txtBlSizeTarget->GetLineText(0).ToLong(&value);
			range.target = value;
			txtBlSizeMin->GetLineText(0).ToLong(&value);
			range.min = value;
			txtBlSizeMax->GetLineText(0).ToLong(&value);
			range.max = value;
			locSelection.SetBlockSizeRange( range );
		}
	}
}

void EditLocusSelector::UpdateMAFSelection() {

	if (txtMafMax && txtMafMin && txtMafTarget) {
		if (initializationComplete) {
			LinearRange<float> range(0.0, 0.0, 0.0);
			double value;
			
			txtMafMax->GetLineText(0).ToDouble(&value);
			range.max = value;
			txtMafMin->GetLineText(0).ToDouble(&value);
			range.min = value;
			txtMafTarget->GetLineText(0).ToDouble(&value);
			range.target = value;
			locSelection.SetMAFRange( range );
		}
	}
}


}

}



