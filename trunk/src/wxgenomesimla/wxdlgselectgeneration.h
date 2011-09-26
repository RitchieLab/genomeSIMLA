/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgselectgeneration.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 27 Feb 2008 15:41:39 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGSELECTGENERATION_H_
#define _WXDLGSELECTGENERATION_H_


/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxPageSelectGeneration;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGSELECTGENERATION 10084
#define ID_GENERATION_SELECTION 10085
#define SYMBOL_WXDLGSELECTGENERATION_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGSELECTGENERATION_TITLE _("wxDlgSelectGeneration")
#define SYMBOL_WXDLGSELECTGENERATION_IDNAME ID_WXDLGSELECTGENERATION
#define SYMBOL_WXDLGSELECTGENERATION_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGSELECTGENERATION_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgSelectGeneration class declaration
 */

class wxDlgSelectGeneration: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgSelectGeneration )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgSelectGeneration();
    wxDlgSelectGeneration( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSELECTGENERATION_IDNAME, const wxString& caption = SYMBOL_WXDLGSELECTGENERATION_TITLE, const wxPoint& pos = SYMBOL_WXDLGSELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXDLGSELECTGENERATION_SIZE, long style = SYMBOL_WXDLGSELECTGENERATION_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSELECTGENERATION_IDNAME, const wxString& caption = SYMBOL_WXDLGSELECTGENERATION_TITLE, const wxPoint& pos = SYMBOL_WXDLGSELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXDLGSELECTGENERATION_SIZE, long style = SYMBOL_WXDLGSELECTGENERATION_STYLE );

    /// Destructor
    ~wxDlgSelectGeneration();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgSelectGeneration event handler declarations
////@end wxDlgSelectGeneration event handler declarations

////@begin wxDlgSelectGeneration member function declarations
    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgSelectGeneration member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgSelectGeneration member variables
    wxPageSelectGeneration* treeGenerations;
////@end wxDlgSelectGeneration member variables
};

#endif
    // _WXDLGSELECTGENERATION_H_
