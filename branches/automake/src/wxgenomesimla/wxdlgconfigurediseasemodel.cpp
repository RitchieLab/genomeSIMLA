/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigurediseasemodel.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 12:56:57 CST
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

#include "wxdlgconfigurediseasemodel.h"
#include <wx/filename.h>



////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgConfigureDiseaseModel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgConfigureDiseaseModel, wxDialog )


/*!
 * wxDlgConfigureDiseaseModel event table definition
 */

BEGIN_EVENT_TABLE( wxDlgConfigureDiseaseModel, wxDialog )

////@begin wxDlgConfigureDiseaseModel event table entries
    EVT_NOTEBOOK_PAGE_CHANGED( ID_NB_DiseaseCfg, wxDlgConfigureDiseaseModel::OnNBDiseaseCfgPageChanged )

    EVT_BUTTON( ID_CMD_EVALUATE, wxDlgConfigureDiseaseModel::OnCmdEvaluateClick )

    EVT_BUTTON( wxID_SAVE, wxDlgConfigureDiseaseModel::OnSaveClick )

    EVT_BUTTON( wxID_SAVEAS, wxDlgConfigureDiseaseModel::OnSaveasClick )

    EVT_BUTTON( wxID_OPEN, wxDlgConfigureDiseaseModel::OnOpenClick )

    EVT_BUTTON( wxID_CANCEL, wxDlgConfigureDiseaseModel::OnCancelClick )

////@end wxDlgConfigureDiseaseModel event table entries

END_EVENT_TABLE()


/*!
 * wxDlgConfigureDiseaseModel constructors
 */

wxDlgConfigureDiseaseModel::wxDlgConfigureDiseaseModel() : newFile(false)
{
    Init();
}

wxDlgConfigureDiseaseModel::wxDlgConfigureDiseaseModel( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : newFile(false)
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgConfigureDiseaseModel creator
 */

bool wxDlgConfigureDiseaseModel::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgConfigureDiseaseModel creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgConfigureDiseaseModel creation
    return true;
}


/*!
 * wxDlgConfigureDiseaseModel destructor
 */

wxDlgConfigureDiseaseModel::~wxDlgConfigureDiseaseModel()
{
////@begin wxDlgConfigureDiseaseModel destruction
////@end wxDlgConfigureDiseaseModel destruction
}


/*!
 * Member initialisation
 */

void wxDlgConfigureDiseaseModel::Init()
{
////@begin wxDlgConfigureDiseaseModel member initialisation
    txtDescription = NULL;
    nbDiseaseConfiguration = NULL;
    cmdEvaluate = NULL;
    cmdSave = NULL;
    cmdSaveAs = NULL;
    cmdOpen = NULL;
    cmdCancel = NULL;
////@end wxDlgConfigureDiseaseModel member initialisation
}


/*!
 * Control creation for wxDlgConfigureDiseaseModel
 */

void wxDlgConfigureDiseaseModel::CreateControls()
{    
////@begin wxDlgConfigureDiseaseModel content construction
    wxDlgConfigureDiseaseModel* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    itemBoxSizer2->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticText* itemStaticText4 = new wxStaticText( itemDialog1, wxID_STATIC, _("Comment:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer2->Add(itemStaticText4, 0, wxALIGN_LEFT|wxALL, 5);

    txtDescription = new wxTextCtrl( itemDialog1, ID_TEXTCTRL_DESCRIPTION, _T(""), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE );
    itemBoxSizer2->Add(txtDescription, 0, wxGROW|wxALL, 5);

    nbDiseaseConfiguration = new wxNotebook( itemDialog1, ID_NB_DiseaseCfg, wxDefaultPosition, wxDefaultSize, wxNB_DEFAULT );

    itemBoxSizer2->Add(nbDiseaseConfiguration, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer7, 0, wxALIGN_RIGHT|wxALL, 5);

    cmdEvaluate = new wxButton( itemDialog1, ID_CMD_EVALUATE, _("Evaluate"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdEvaluate, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdSave = new wxButton( itemDialog1, wxID_SAVE, _("&Save"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdSave, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdSaveAs = new wxButton( itemDialog1, wxID_SAVEAS, _("Save &As..."), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdSaveAs, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdOpen = new wxButton( itemDialog1, wxID_OPEN, _("&Open"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdOpen, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdCancel = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer7->Add(cmdCancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgConfigureDiseaseModel content construction
	pagePenetrance = new wxPagePenetranceModel(nbDiseaseConfiguration, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER | wxTAB_TRAVERSAL);
	nbDiseaseConfiguration->AddPage(pagePenetrance, _T("Penetrance"));

	pageSimpen = new wxPageSimPEN(nbDiseaseConfiguration, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER | wxTAB_TRAVERSAL);
	nbDiseaseConfiguration->AddPage(pageSimpen, _T("simPEN"));

	pageSimla = new wxDlgSIMLA(nbDiseaseConfiguration, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER | wxTAB_TRAVERSAL);
	nbDiseaseConfiguration->AddPage(pageSimla, _T("SIMLA"));

}


/*!
 * Should we show tooltips?
 */

bool wxDlgConfigureDiseaseModel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgConfigureDiseaseModel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgConfigureDiseaseModel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgConfigureDiseaseModel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgConfigureDiseaseModel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgConfigureDiseaseModel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgConfigureDiseaseModel icon retrieval
}


/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for wxID_SAVE
 */

void wxDlgConfigureDiseaseModel::OnSaveClick( wxCommandEvent& event )
{
	int sel = nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);

	if (cfg) {
		if (filename.length() == 0) {
			string fileTypes = cfg->GetImportFileTypes();
			wxFileDialog filenameSelect(this, wxT("Save As"), wxT(""), wxT(""), wxT(fileTypes.c_str()), wxSAVE|wxOVERWRITE_PROMPT);
			if (filenameSelect.ShowModal() == wxID_OK) {
				wxFileName newFilename(filenameSelect.GetPath());
				if (newFilename.GetExt().Len() == 0)
					newFilename.SetExt(_(cfg->GetPreferredExtension().c_str()));
			
				filename=newFilename.GetFullPath().c_str();
				newFile=true;			//Let the calling function know that we created a new source file
			}
			else
				return;
		}
		//If we have an error, we don't actually want to close the dialog
		if (!cfg->Save(filename.c_str(), txtDescription->GetValue().c_str()))
			return;
	}
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}

int wxDlgConfigureDiseaseModel::GetLocusCount() {
	int locusCount = 0;
	int sel = nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);
	if (cfg)
		locusCount = cfg->GetLocusCount();
	return locusCount;
}

/*!
 * wxEVT_COMMAND_MENU_SELECTED event handler for wxID_SAVEAS
 */

void wxDlgConfigureDiseaseModel::OnSaveasClick( wxCommandEvent& event )
{
	newFile = true;
	int sel = nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);
	if (cfg) {
		string fileTypes = cfg->GetImportFileTypes();

		wxFileDialog filenameSelect(this, wxT("Save As"), wxT(""), wxT(""), wxT(fileTypes.c_str()), wxSAVE|wxOVERWRITE_PROMPT);
		if (filenameSelect.ShowModal() == wxID_OK) {
			wxFileName newFilename(filenameSelect.GetPath());
			if (newFilename.GetExt().Len() == 0)
				newFilename.SetExt(_(cfg->GetPreferredExtension().c_str()));
			filename=newFilename.GetFullPath().c_str();
			cfg->Save(filename.c_str(), txtDescription->GetValue().c_str());
		}		

		cout<<"SaveAsing\n";
	}
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}



string wxDlgConfigureDiseaseModel::GetLabel() {
	wxFileName label(filename.c_str());
	return label.GetFullName().c_str();
}

bool wxDlgConfigureDiseaseModel::LoadModel(const char *filetype, const char *filename, bool editable) {
	this->filename=filename;
	bool success=false;
	newFile = false;
	int sel = -1;
	if (strcmp(filetype, "PENTABLE") == 0) 
		sel = 0;
	else if (strcmp(filetype, "SIMPEN") == 0)
		sel = 1;
	else if (strcmp(filetype, "SIMLA") == 0)
		sel = 2;
	if (sel > -1) {
		nbDiseaseConfiguration->SetSelection(sel);
		wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);

		cmdSave->Enable(editable);
//		toolbar->EnableTool(wxID_SAVE, editable);

		if (cfg) {
			string desc;
			cfg->Import(filename, desc);
			success= true;
		}
	}
	return success;
}



/*!
 * wxEVT_COMMAND_NOTEBOOK_PAGE_CHANGED event handler for ID_NB_DiseaseCfg
 */

void wxDlgConfigureDiseaseModel::OnNBDiseaseCfgPageChanged( wxNotebookEvent& event )
{
	static string types[]={"PENTABLE", "SIMPEN", "SIMLA"};
	int sel =  nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);
	cout<<"Page changing "<<sel<<":\n";

	if (cfg) {
		filetype=types[sel];
		cmdSave->Enable(cfg->CanSave());
		cmdSaveAs->Enable(cfg->CanSaveAs());
		cmdOpen->Enable(cfg->CanImport());
		cmdEvaluate->Enable(cfg->CanEvaluate());

	}
}




bool wxGridDiseaseAssignment::IsEmptyCell(int row, int col) {
	bool isEmpty=true;
	if ((size_t)row < loci.size() ) 
		isEmpty = !(loci[row].l.Valid());

	return isEmpty;
}

wxString wxGridDiseaseAssignment::GetValue(int row, int col) { 
	wxString rVal = wxT("");

	if ((size_t)row < loci.size()) {
		switch (col) {
			case 0:
				rVal = loci[row].l.GetLabel().c_str();
				break;
			case 1:
				rVal.Printf("%f", loci[row].l.Freq1());
				break;
			case 2:
				rVal.Printf("%f", loci[row].l.Freq2());
				break;
			case 3:
				rVal = loci[row].details.c_str();
				break;
		}
	}

	return rVal;
};

void wxGridDiseaseAssignment::SetValue(int row, int col, const wxString& value) {
	wxString rVal;
	assert((size_t)row<loci.size());
	
	wxGridDiseaseAssignment::LocusEntry &entry = loci[row];
	switch (col) {
		case 0:
			
			rVal = loci[row].l.GetLabel().c_str();
			break;
		case 1:
			rVal.Printf("%f", loci[row].l.Freq1());
			break;
		case 2:
			rVal.Printf("%f", loci[row].l.Freq2());
			break;
		case 3:
			rVal = loci[row].details.c_str();
			break;
	}
	
}
wxString wxGridDiseaseAssignment::GetColLabelValue( int col ) {
	assert(0);
}


wxString wxGridDiseaseAssignment::GetRowLabelValue( int row ) {
	assert(0);
}

wxString wxGridDiseaseAssignment::GetTypeName( int row, int col ) {
	assert(0);
}




/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
 */

void wxDlgConfigureDiseaseModel::OnCancelClick( wxCommandEvent& event )
{
	if ( IsModal() )
		EndModal(wxID_CANCEL);
	else {
		SetReturnCode(wxID_CANCEL);
		this->Show(false);
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OPEN
 */

void wxDlgConfigureDiseaseModel::OnOpenClick( wxCommandEvent& event )
{
	int sel = nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);
	newFile = true;
	if (cfg) {
		string fileTypes = cfg->GetImportFileTypes();
		cout<<"fileTypes("<<sel<<"): "<<fileTypes<<"\n";
		wxFileDialog filenameSelect(this, wxT("Import"), wxT(""), wxT(""), wxT(fileTypes.c_str()), wxOPEN);
		if (filenameSelect.ShowModal() == wxID_OK) {
			wxString path = filenameSelect.GetPath();
			string desc;
			this->filename=path.c_str();
			cfg->Import(path.c_str(), desc);
		}

	}
}



/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_EVALUATE
 */

void wxDlgConfigureDiseaseModel::OnCmdEvaluateClick( wxCommandEvent& event )
{
	int sel = nbDiseaseConfiguration->GetSelection();
	wxDiseaseConfiguration *cfg = (wxDiseaseConfiguration*)nbDiseaseConfiguration->GetPage(sel);
	if (cfg)
		cfg->Evaluate();
}


}

}




