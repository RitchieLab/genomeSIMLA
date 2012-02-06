/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgprocessmonitor.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 18 Feb 2008 17:43:42 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGPROCESSMONITOR_H_
#define _WXDLGPROCESSMONITOR_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/statusbr.h"
////@end includes

#include "simwizard.h"


/*!
 * Forward declarations
 */

namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGPROCESSMONITOR 10086
#define ID_TXT_SIMLOG 10002
#define ID_GUAGE_PROCESS_COMPLETION 10003
#define ID_CMD_RUN 10004
#define ID_STATUSBAR 10005
#define SYMBOL_WXDLGPROCESSMONITOR_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGPROCESSMONITOR_TITLE _("Process Monitor")
#define SYMBOL_WXDLGPROCESSMONITOR_IDNAME ID_WXDLGPROCESSMONITOR
#define SYMBOL_WXDLGPROCESSMONITOR_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGPROCESSMONITOR_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgProcessMonitor class declaration
 */

class wxDlgProcessMonitor: public wxDialog, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxDlgProcessMonitor )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgProcessMonitor();
    wxDlgProcessMonitor( wxWindow* parent, AppController *ctrl, wxWindowID id = SYMBOL_WXDLGPROCESSMONITOR_IDNAME, const wxString& caption = SYMBOL_WXDLGPROCESSMONITOR_TITLE, const wxPoint& pos = SYMBOL_WXDLGPROCESSMONITOR_POSITION, const wxSize& size = SYMBOL_WXDLGPROCESSMONITOR_SIZE, long style = SYMBOL_WXDLGPROCESSMONITOR_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGPROCESSMONITOR_IDNAME, const wxString& caption = SYMBOL_WXDLGPROCESSMONITOR_TITLE, const wxPoint& pos = SYMBOL_WXDLGPROCESSMONITOR_POSITION, const wxSize& size = SYMBOL_WXDLGPROCESSMONITOR_SIZE, long style = SYMBOL_WXDLGPROCESSMONITOR_STYLE );

    /// Destructor
    ~wxDlgProcessMonitor();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgProcessMonitor event handler declarations

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_RUN
    void OnCmdRunClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
    void OnCancelClick( wxCommandEvent& event );

////@end wxDlgProcessMonitor event handler declarations

////@begin wxDlgProcessMonitor member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgProcessMonitor member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgProcessMonitor member variables
    wxTextCtrl* txtSimLog;
    wxGauge* guageSimCompletion;
    wxButton* cmdRun;
    wxButton* cmdCancel;
    wxStatusBar* statusBar;
////@end wxDlgProcessMonitor member variables
protected:
	virtual bool RunAnalysis();
		
	//The main simulation runs in it's own thread
	pthread_t simThread;
	bool continueRunning;
	bool isCompleted;
	
};

class wxDlgDatasetGenerationMonitor: public wxDlgProcessMonitor {
    DECLARE_DYNAMIC_CLASS( wxDlgProcessMonitor )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgDatasetGenerationMonitor() {}
    wxDlgDatasetGenerationMonitor( wxWindow* parent, AppController *ctrl, wxWindowID id = SYMBOL_WXDLGPROCESSMONITOR_IDNAME, const wxString& caption = SYMBOL_WXDLGPROCESSMONITOR_TITLE, const wxPoint& pos = SYMBOL_WXDLGPROCESSMONITOR_POSITION, const wxSize& size = SYMBOL_WXDLGPROCESSMONITOR_SIZE, long style = SYMBOL_WXDLGPROCESSMONITOR_STYLE );

	bool RunAnalysis();

};
}
}
#endif
    // _WXDLGPROCESSMONITOR_H_
