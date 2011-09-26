/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigureccdatasets.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 21 Mar 2008 12:07:48 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGCONFIGURECCDATASETS_H_
#define _WXDLGCONFIGURECCDATASETS_H_

#include <string>
#include "wxuser.h"
/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGCONFIGURECCDATASETS 10110
#define ID_TEXTCTRL1 10001
#define ID_TEXTCTRL4 10111
#define ID_TEXTCTRL5 10112
#define ID_TEXTCTRL6 10114
#define ID_TEXTCTRL7 10115
#define ID_TEXTCTRL8 10117
#define SYMBOL_WXDLGCONFIGURECCDATASETS_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGCONFIGURECCDATASETS_TITLE _("Configure Case-Control Datasets")
#define SYMBOL_WXDLGCONFIGURECCDATASETS_IDNAME ID_WXDLGCONFIGURECCDATASETS
#define SYMBOL_WXDLGCONFIGURECCDATASETS_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGCONFIGURECCDATASETS_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgConfigureCCDatasets class declaration
 */

class wxDlgConfigureCCDatasets: public wxDialog, public WxUser
{    
    DECLARE_DYNAMIC_CLASS( wxDlgConfigureCCDatasets )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgConfigureCCDatasets();
    wxDlgConfigureCCDatasets( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGURECCDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGURECCDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGURECCDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGURECCDATASETS_SIZE, long style = SYMBOL_WXDLGCONFIGURECCDATASETS_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGURECCDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGURECCDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGURECCDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGURECCDATASETS_SIZE, long style = SYMBOL_WXDLGCONFIGURECCDATASETS_STYLE );

    /// Destructor
    ~wxDlgConfigureCCDatasets();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgConfigureCCDatasets event handler declarations

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
    void OnCancelClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
    void OnOkClick( wxCommandEvent& event );

////@end wxDlgConfigureCCDatasets event handler declarations

////@begin wxDlgConfigureCCDatasets member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgConfigureCCDatasets member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgConfigureCCDatasets member variables
    wxTextCtrl* txtLabel;
    wxTextCtrl* txtGenotypeErr;
    wxTextCtrl* txtPhenocopy;
    wxTextCtrl* txtMissingData;
    wxTextCtrl* txtCases;
    wxTextCtrl* txtControls;
////@end wxDlgConfigureCCDatasets member variables

	string GetModelDetails();
	string GetModelLabel();
	void ConfigureModelDetails(const char *details);
	string GetModelSummary();
};

}
}
#endif
    // _WXDLGCONFIGURECCDATASETS_H_
