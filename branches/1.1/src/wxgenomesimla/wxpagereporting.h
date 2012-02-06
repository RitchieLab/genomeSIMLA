/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagereporting.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 21 Feb 2008 10:15:24 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGEREPORTING_H_
#define _WXPAGEREPORTING_H_

#include "appinterface.h"
#include "wxpanelreports.h"

/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */



namespace GenomeSIM {

namespace GUI {

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXPAGEREPORTING 10087
#define ID_PAGESELECTGENERATION 10000
#define SYMBOL_WXPAGEREPORTING_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGEREPORTING_TITLE _("Reporting")
#define SYMBOL_WXPAGEREPORTING_IDNAME ID_WXPAGEREPORTING
#define SYMBOL_WXPAGEREPORTING_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGEREPORTING_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxPageReporting class declaration
 */

class wxPageReporting: public wxPanel, public AppInterface {    
    DECLARE_DYNAMIC_CLASS( wxPageReporting )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageReporting();
    wxPageReporting( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEREPORTING_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEREPORTING_POSITION, const wxSize& size = SYMBOL_WXPAGEREPORTING_SIZE, long style = SYMBOL_WXPAGEREPORTING_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEREPORTING_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEREPORTING_POSITION, const wxSize& size = SYMBOL_WXPAGEREPORTING_SIZE, long style = SYMBOL_WXPAGEREPORTING_STYLE );

    /// Destructor
    ~wxPageReporting();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageReporting event handler declarations

    /// wxEVT_INIT_DIALOG event handler for ID_WXPAGEREPORTING
    void OnInitDialog( wxInitDialogEvent& event );

    /// wxEVT_SIZE event handler for ID_WXPAGEREPORTING
    void OnSize( wxSizeEvent& event );

////@end wxPageReporting event handler declarations

////@begin wxPageReporting member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageReporting member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageReporting member variables
    wxPanelReports* treeGenerations;
////@end wxPageReporting member variables

	void InitAppController(AppController *ctrl);

	/**
	 * @brief This is called prior to saving the configuration
	 */
	void Commit();

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings();

	bool HasChanged();
	void OnContinueSimulation(wxCommandEvent& event);
	void OnDetailedAnalysis(wxCommandEvent& event);
	void OpenReport(wxCommandEvent &event);
	void ExtractData(wxCommandEvent &event);
};

}

}

#endif
    // _WXPAGEREPORTING_H_
