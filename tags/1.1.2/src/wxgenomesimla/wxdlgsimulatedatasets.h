/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsimulatedatasets.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 10:18:34 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGSIMULATEDATASETS_H_
#define _WXDLGSIMULATEDATASETS_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/listctrl.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxListCtrl;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGSIMULATEDATASETS 10118
#define ID_TEXTCTRL 10000
#define ID_LIST_DATASETS 10006
#define ID_CHK_PED_HEADER 10002
#define ID_CHK_BINARY_DATASET 10003
#define ID_CMD_ADD_CC 10004
#define ID_CMD_ADD_PED 10005
#define ID_CMD_DEL_DATASET 10122
#define ID_CHOICE_DISEASE 10119
#define ID_CMD_CONFIG_MODEL 10120
#define SYMBOL_WXDLGSIMULATEDATASETS_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGSIMULATEDATASETS_TITLE _("Simulate Datasets")
#define SYMBOL_WXDLGSIMULATEDATASETS_IDNAME ID_WXDLGSIMULATEDATASETS
#define SYMBOL_WXDLGSIMULATEDATASETS_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGSIMULATEDATASETS_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgSimulateDatasets class declaration
 */

class wxDlgSimulateDatasets: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgSimulateDatasets )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgSimulateDatasets();
    wxDlgSimulateDatasets( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSIMULATEDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGSIMULATEDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGSIMULATEDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGSIMULATEDATASETS_SIZE, long style = SYMBOL_WXDLGSIMULATEDATASETS_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSIMULATEDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGSIMULATEDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGSIMULATEDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGSIMULATEDATASETS_SIZE, long style = SYMBOL_WXDLGSIMULATEDATASETS_STYLE );

    /// Destructor
    ~wxDlgSimulateDatasets();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgSimulateDatasets event handler declarations
    /// wxEVT_SIZE event handler for ID_WXDLGSIMULATEDATASETS
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LIST_DATASETS
    void OnListDatasetsSelected( wxListEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LIST_DATASETS
    void OnListDatasetsDeselected( wxListEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LIST_DATASETS
    void OnListDatasetsItemActivated( wxListEvent& event );

    /// wxEVT_LEFT_DCLICK event handler for ID_LIST_DATASETS
    void OnLeftDClick( wxMouseEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC
    void OnCmdAddCcClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED
    void OnCmdAddPedClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_DEL_DATASET
    void OnCmdDelDatasetClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CHOICE_DISEASE
    void OnChoiceDiseaseSelected( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL
    void OnCmdConfigModelClick( wxCommandEvent& event );

////@end wxDlgSimulateDatasets event handler declarations

////@begin wxDlgSimulateDatasets member function declarations
    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgSimulateDatasets member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgSimulateDatasets member variables
    wxTextCtrl* txtDatasetCount;
    wxListCtrl* lstDatasets;
    wxCheckBox* chkStandardPedigreeHeader;
    wxChoice* lstDiseaseModels;
////@end wxDlgSimulateDatasets member variables
};

#endif
    // _WXDLGSIMULATEDATASETS_H_
