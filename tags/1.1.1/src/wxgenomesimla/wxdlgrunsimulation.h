/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgrunsimulation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 14 Jan 2008 12:46:57 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGRUNSIMULATION_H_
#define _WXDLGRUNSIMULATION_H_


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
#define ID_WXDLGRUNSIMULATION 10059
#define ID_TEXTCTRL2 10060
#define ID_TEXTCTRL 10000
#define ID_TEXTCTRL1 10001
#define SYMBOL_WXDLGRUNSIMULATION_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGRUNSIMULATION_TITLE _("Run Simulation")
#define SYMBOL_WXDLGRUNSIMULATION_IDNAME ID_WXDLGRUNSIMULATION
#define SYMBOL_WXDLGRUNSIMULATION_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGRUNSIMULATION_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgRunSimulation class declaration
 */

class wxDlgRunSimulation: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgRunSimulation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgRunSimulation();
    wxDlgRunSimulation( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGRUNSIMULATION_IDNAME, const wxString& caption = SYMBOL_WXDLGRUNSIMULATION_TITLE, const wxPoint& pos = SYMBOL_WXDLGRUNSIMULATION_POSITION, const wxSize& size = SYMBOL_WXDLGRUNSIMULATION_SIZE, long style = SYMBOL_WXDLGRUNSIMULATION_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGRUNSIMULATION_IDNAME, const wxString& caption = SYMBOL_WXDLGRUNSIMULATION_TITLE, const wxPoint& pos = SYMBOL_WXDLGRUNSIMULATION_POSITION, const wxSize& size = SYMBOL_WXDLGRUNSIMULATION_SIZE, long style = SYMBOL_WXDLGRUNSIMULATION_STYLE );

    /// Destructor
    ~wxDlgRunSimulation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgRunSimulation event handler declarations

////@end wxDlgRunSimulation event handler declarations

////@begin wxDlgRunSimulation member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgRunSimulation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgRunSimulation member variables
////@end wxDlgRunSimulation member variables
};

#endif
    // _WXDLGRUNSIMULATION_H_
