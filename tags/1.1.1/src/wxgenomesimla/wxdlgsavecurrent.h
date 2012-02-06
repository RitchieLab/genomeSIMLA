/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsavecurrent.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Dec 2007 14:47:44 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGSAVECURRENT_H_
#define _WXDLGSAVECURRENT_H_


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

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGSAVECURRENT 10036
#define SYMBOL_WXDLGSAVECURRENT_STYLE wxCAPTION|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGSAVECURRENT_TITLE _("Save Current Configuration")
#define SYMBOL_WXDLGSAVECURRENT_IDNAME ID_WXDLGSAVECURRENT
#define SYMBOL_WXDLGSAVECURRENT_SIZE wxSize(200, 350)
#define SYMBOL_WXDLGSAVECURRENT_POSITION wxDefaultPosition
////@end control identifiers
namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgSaveCurrent class declaration
 */

class wxDlgSaveCurrent: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgSaveCurrent )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgSaveCurrent();
    wxDlgSaveCurrent( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSAVECURRENT_IDNAME, const wxString& caption = SYMBOL_WXDLGSAVECURRENT_TITLE, const wxPoint& pos = SYMBOL_WXDLGSAVECURRENT_POSITION, const wxSize& size = SYMBOL_WXDLGSAVECURRENT_SIZE, long style = SYMBOL_WXDLGSAVECURRENT_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSAVECURRENT_IDNAME, const wxString& caption = SYMBOL_WXDLGSAVECURRENT_TITLE, const wxPoint& pos = SYMBOL_WXDLGSAVECURRENT_POSITION, const wxSize& size = SYMBOL_WXDLGSAVECURRENT_SIZE, long style = SYMBOL_WXDLGSAVECURRENT_STYLE );

    /// Destructor
    ~wxDlgSaveCurrent();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgSaveCurrent event handler declarations

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CLOSE
    void OnCloseClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVEAS
    void OnSaveAsClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVE
    void OnSaveClick( wxCommandEvent& event );

////@end wxDlgSaveCurrent event handler declarations

////@begin wxDlgSaveCurrent member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgSaveCurrent member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgSaveCurrent member variables
////@end wxDlgSaveCurrent member variables
};

}

}

#endif
    // _WXDLGSAVECURRENT_H_
