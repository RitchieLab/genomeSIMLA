/////////////////////////////////////////////////////////////////////////////
// Name:        wxbeginsimulation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 04 Feb 2008 15:56:51 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXBEGINSIMULATION_H_
#define _WXBEGINSIMULATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/statusbr.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxStatusBar;
////@end forward declarations
#include "appcontroller.h"

namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_BEGINSIMULATION 10072
#define ID_TXT_SIMLOG 10000
#define ID_GAUGE 10001
#define ID_CMD_RUN 10073
#define ID_STATUS_BAR 10074
#define SYMBOL_WXBEGINSIMULATION_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXBEGINSIMULATION_TITLE _("Begin Simulation")
#define SYMBOL_WXBEGINSIMULATION_IDNAME ID_BEGINSIMULATION
#define SYMBOL_WXBEGINSIMULATION_SIZE wxSize(400, 300)
#define SYMBOL_WXBEGINSIMULATION_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxBeginSimulation class declaration
 */

class wxBeginSimulation: public wxDialog 	{    
    DECLARE_DYNAMIC_CLASS( wxBeginSimulation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxBeginSimulation();
    wxBeginSimulation( wxWindow* parent, wxWindowID id = SYMBOL_WXBEGINSIMULATION_IDNAME, const wxString& caption = SYMBOL_WXBEGINSIMULATION_TITLE, const wxPoint& pos = SYMBOL_WXBEGINSIMULATION_POSITION, const wxSize& size = SYMBOL_WXBEGINSIMULATION_SIZE, long style = SYMBOL_WXBEGINSIMULATION_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXBEGINSIMULATION_IDNAME, const wxString& caption = SYMBOL_WXBEGINSIMULATION_TITLE, const wxPoint& pos = SYMBOL_WXBEGINSIMULATION_POSITION, const wxSize& size = SYMBOL_WXBEGINSIMULATION_SIZE, long style = SYMBOL_WXBEGINSIMULATION_STYLE );

    /// Destructor
    ~wxBeginSimulation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxBeginSimulation event handler declarations

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_RUN
    void OnCmdRunClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
    void OnCancelClick( wxCommandEvent& event );

////@end wxBeginSimulation event handler declarations

////@begin wxBeginSimulation member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxBeginSimulation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxBeginSimulation member variables
    wxTextCtrl* txtSimLog;
    wxGauge* guageSimCompletion;
    wxButton* cmdRun;
    wxButton* cmdCancel;
    wxStatusBar* statusBar;
////@end wxBeginSimulation member variables


	bool Initialize(AppController *appController);
	bool RunSimulation(bool doLoadFirst);
	bool isCompleted;
protected:
	AppController *appController;

	//The main simulation runs in it's own thread
	pthread_t simThread;
	bool continueRunning;
};

}

}

#endif
    // _WXBEGINSIMULATION_H_
