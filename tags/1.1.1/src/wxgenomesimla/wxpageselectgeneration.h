/////////////////////////////////////////////////////////////////////////////
// Name:        wxpageselectgeneration.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 15 Feb 2008 11:08:28 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGESELECTGENERATION_H_
#define _WXPAGESELECTGENERATION_H_


/*!
 * Includes
 */

#include "simwizard.h"
#include "treelistctrl.h"
#include "utility/executionlog.h"
/*!
 * Forward declarations
 */

/*!
 * Control identifiers
 */

namespace GenomeSIM {

namespace GUI {

////@begin control identifiers
#define ID_WXPAGESELECTGENERATION 10076
#define ID_TREE_GEN_SELECTION 10077
#define SYMBOL_WXPAGESELECTGENERATION_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGESELECTGENERATION_TITLE _("Select Generation")
#define SYMBOL_WXPAGESELECTGENERATION_IDNAME ID_WXPAGESELECTGENERATION
#define SYMBOL_WXPAGESELECTGENERATION_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGESELECTGENERATION_POSITION wxDefaultPosition
////@end control identifiers

using namespace Utility;
/*!
 * wxPageSelectGeneration class declaration
 */

class LogEntryData : public wxTreeItemData {
public:
	LogEntryData(ExecutionLog::LogEntry &entry) : entry(entry) { }
	ExecutionLog::LogEntry entry;
};

class wxPageSelectGeneration: public wxPanel, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxPageSelectGeneration )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageSelectGeneration();
    wxPageSelectGeneration( wxWindow* parent, AppController *ctrl, wxWindowID id = SYMBOL_WXPAGESELECTGENERATION_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXPAGESELECTGENERATION_SIZE, long style = SYMBOL_WXPAGESELECTGENERATION_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGESELECTGENERATION_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXPAGESELECTGENERATION_SIZE, long style = SYMBOL_WXPAGESELECTGENERATION_STYLE );

    /// Destructor
    ~wxPageSelectGeneration();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageSelectGeneration event handler declarations

    /// wxEVT_SIZE event handler for ID_WXPAGESELECTGENERATION
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_TREE_SEL_CHANGED event handler for ID_TREE_GEN_SELECTION
    void OnTreeGenSelectionSelChanged( wxTreeEvent& event );

    /// wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREE_GEN_SELECTION
    void OnTreeGenSelectionItemActivated( wxTreeEvent& event );

    /// wxEVT_COMMAND_TREE_ITEM_EXPANDED event handler for ID_TREE_GEN_SELECTION
    void OnTreeGenSelectionItemExpanded( wxTreeEvent& event );

////@end wxPageSelectGeneration event handler declarations

////@begin wxPageSelectGeneration member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageSelectGeneration member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageSelectGeneration member variables
    wxTreeCtrl* treeGenerationSelection;
////@end wxPageSelectGeneration member variables


	void AddProjectEntries(const char *project, ExecutionLog::RunType *entries);
	void InitializeColumns();

	bool ValidSelection();

	/**
	 * @brief Get the generation selected by the user. Returns t/f indicating that the value was changed
	 */
	bool GetSelectedGeneration(ExecutionLog::LogEntry &entry);

	bool TransferDataFromWindow();

	ExecutionLog::LogEntry *entry;
protected:
	wxTreeItemId root;	

};

}

}
#endif
    // _WXPAGESELECTGENERATION_H_
