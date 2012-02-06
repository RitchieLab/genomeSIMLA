/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsimulatedatasets.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 10:18:34 CDT
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

#include "wxdlgsimulatedatasets.h"

////@begin XPM images

////@end XPM images


/*!
 * wxDlgSimulateDatasets type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgSimulateDatasets, wxDialog )


/*!
 * wxDlgSimulateDatasets event table definition
 */

BEGIN_EVENT_TABLE( wxDlgSimulateDatasets, wxDialog )

////@begin wxDlgSimulateDatasets event table entries
    EVT_SIZE( wxDlgSimulateDatasets::OnSize )

    EVT_LIST_ITEM_SELECTED( ID_LIST_DATASETS, wxDlgSimulateDatasets::OnListDatasetsSelected )
    EVT_LIST_ITEM_DESELECTED( ID_LIST_DATASETS, wxDlgSimulateDatasets::OnListDatasetsDeselected )
    EVT_LIST_ITEM_ACTIVATED( ID_LIST_DATASETS, wxDlgSimulateDatasets::OnListDatasetsItemActivated )

    EVT_BUTTON( ID_CMD_ADD_CC, wxDlgSimulateDatasets::OnCmdAddCcClick )

    EVT_BUTTON( ID_CMD_ADD_PED, wxDlgSimulateDatasets::OnCmdAddPedClick )

    EVT_BUTTON( ID_CMD_DEL_DATASET, wxDlgSimulateDatasets::OnCmdDelDatasetClick )

    EVT_CHOICE( ID_CHOICE_DISEASE, wxDlgSimulateDatasets::OnChoiceDiseaseSelected )

    EVT_BUTTON( ID_CMD_CONFIG_MODEL, wxDlgSimulateDatasets::OnCmdConfigModelClick )

////@end wxDlgSimulateDatasets event table entries

END_EVENT_TABLE()


/*!
 * wxDlgSimulateDatasets constructors
 */

wxDlgSimulateDatasets::wxDlgSimulateDatasets()
{
    Init();
}

wxDlgSimulateDatasets::wxDlgSimulateDatasets( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgSimulateDatasets creator
 */

bool wxDlgSimulateDatasets::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgSimulateDatasets creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgSimulateDatasets creation
    return true;
}


/*!
 * wxDlgSimulateDatasets destructor
 */

wxDlgSimulateDatasets::~wxDlgSimulateDatasets()
{
////@begin wxDlgSimulateDatasets destruction
////@end wxDlgSimulateDatasets destruction
}


/*!
 * Member initialisation
 */

void wxDlgSimulateDatasets::Init()
{
////@begin wxDlgSimulateDatasets member initialisation
    txtDatasetCount = NULL;
    lstDatasets = NULL;
    chkStandardPedigreeHeader = NULL;
    lstDiseaseModels = NULL;
////@end wxDlgSimulateDatasets member initialisation
}


/*!
 * Control creation for wxDlgSimulateDatasets
 */

void wxDlgSimulateDatasets::CreateControls()
{    
////@begin wxDlgSimulateDatasets content construction
    wxDlgSimulateDatasets* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer4, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer4->Add(itemBoxSizer5, 0, wxALIGN_RIGHT|wxALL, 0);

    wxStaticText* itemStaticText6 = new wxStaticText( itemDialog1, wxID_STATIC, _("Number of Data Sets:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    itemBoxSizer5->Add(itemStaticText6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtDatasetCount = new wxTextCtrl( itemDialog1, ID_TEXTCTRL, _("1000"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer5->Add(txtDatasetCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    lstDatasets = new wxListCtrl( itemDialog1, ID_LIST_DATASETS, wxDefaultPosition, wxSize(300, 100), wxLC_REPORT );
    itemBoxSizer4->Add(lstDatasets, 1, wxGROW|wxALL, 5);

    chkStandardPedigreeHeader = new wxCheckBox( itemDialog1, ID_CHK_PED_HEADER, _("Use Standard Pedigree Header"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    chkStandardPedigreeHeader->SetValue(true);
    itemBoxSizer4->Add(chkStandardPedigreeHeader, 0, wxALIGN_RIGHT|wxALL, 5);

    wxCheckBox* itemCheckBox10 = new wxCheckBox( itemDialog1, ID_CHK_BINARY_DATASET, _("Write Binary Datasets"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    itemCheckBox10->SetValue(false);
    itemCheckBox10->SetHelpText(_("Save space with binary datasets, but only if your application can read them"));
    if (wxDlgSimulateDatasets::ShowToolTips())
        itemCheckBox10->SetToolTip(_("Save space with binary datasets, but only if your application can read them"));
    itemBoxSizer4->Add(itemCheckBox10, 0, wxALIGN_RIGHT|wxALL, 5);

    wxBoxSizer* itemBoxSizer11 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer4->Add(itemBoxSizer11, 0, wxGROW|wxALL, 0);

    itemBoxSizer4->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer13, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton14 = new wxButton( itemDialog1, ID_CMD_ADD_CC, _("Add &Case/Control"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton14, 0, wxGROW|wxALL, 5);

    wxButton* itemButton15 = new wxButton( itemDialog1, ID_CMD_ADD_PED, _("Add &Pedigree"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton15, 0, wxGROW|wxALL, 5);

    wxButton* itemButton16 = new wxButton( itemDialog1, ID_CMD_DEL_DATASET, _("&Delete Dataset"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton16, 0, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer17Static = new wxStaticBox(itemDialog1, wxID_ANY, _("Status Model"));
    wxStaticBoxSizer* itemStaticBoxSizer17 = new wxStaticBoxSizer(itemStaticBoxSizer17Static, wxHORIZONTAL);
    itemBoxSizer2->Add(itemStaticBoxSizer17, 0, wxGROW|wxALL, 5);

    wxArrayString lstDiseaseModelsStrings;
    lstDiseaseModels = new wxChoice( itemDialog1, ID_CHOICE_DISEASE, wxDefaultPosition, wxDefaultSize, lstDiseaseModelsStrings, 0 );
    itemStaticBoxSizer17->Add(lstDiseaseModels, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton19 = new wxButton( itemDialog1, ID_CMD_CONFIG_MODEL, _("Confi&gure"), wxDefaultPosition, wxDefaultSize, 0 );
    itemStaticBoxSizer17->Add(itemButton19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer20 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer20, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton21 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer20->Add(itemButton21, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton22 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer20->Add(itemButton22, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);


    // Connect events and objects
    lstDatasets->Connect(ID_LIST_DATASETS, wxEVT_LEFT_DCLICK, wxMouseEventHandler(wxDlgSimulateDatasets::OnLeftDClick), NULL, this);
////@end wxDlgSimulateDatasets content construction
}


/*!
 * wxEVT_SIZE event handler for ID_WXDLGSIMULATEDATASETS
 */

void wxDlgSimulateDatasets::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_WXDLGSIMULATEDATASETS in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_WXDLGSIMULATEDATASETS in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LIST_DATASETS
 */

void wxDlgSimulateDatasets::OnListDatasetsSelected( wxListEvent& event )
{
////@begin wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LIST_DATASETS
 */

void wxDlgSimulateDatasets::OnListDatasetsDeselected( wxListEvent& event )
{
////@begin wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LIST_DATASETS
 */

void wxDlgSimulateDatasets::OnListDatasetsItemActivated( wxListEvent& event )
{
////@begin wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_LEFT_DCLICK event handler for ID_LIST_DATASETS
 */

void wxDlgSimulateDatasets::OnLeftDClick( wxMouseEvent& event )
{
////@begin wxEVT_LEFT_DCLICK event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_LEFT_DCLICK event handler for ID_LIST_DATASETS in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC
 */

void wxDlgSimulateDatasets::OnCmdAddCcClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED
 */

void wxDlgSimulateDatasets::OnCmdAddPedClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_DEL_DATASET
 */

void wxDlgSimulateDatasets::OnCmdDelDatasetClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_DEL_DATASET in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_DEL_DATASET in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CHOICE_DISEASE
 */

void wxDlgSimulateDatasets::OnChoiceDiseaseSelected( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CHOICE_DISEASE in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CHOICE_DISEASE in wxDlgSimulateDatasets. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL
 */

void wxDlgSimulateDatasets::OnCmdConfigModelClick( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL in wxDlgSimulateDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL in wxDlgSimulateDatasets. 
}


/*!
 * Should we show tooltips?
 */

bool wxDlgSimulateDatasets::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgSimulateDatasets::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgSimulateDatasets bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgSimulateDatasets bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgSimulateDatasets::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgSimulateDatasets icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgSimulateDatasets icon retrieval
}
