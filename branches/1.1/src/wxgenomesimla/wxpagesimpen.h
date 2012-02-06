/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagesimpen.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 26 Feb 2008 13:23:04 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGESIMPEN_H_
#define _WXPAGESIMPEN_H_


/*!
 * Includes
 */

////@begin includes
////@end includes
#include "wxdiseaseconfiguration.h"
#include "utility/lineparser.h"
#include "wx/grid.h"
/*!
 * Forward declarations
 */

////@begin forward declarations
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace Utility;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXPAGESIMPEN 10100
#define ID_TEXTCTRL 10006
#define ID_TXT_HERIT_TARGET 10102
#define ID_TXT_HERIT_WEIGHT 10000
#define ID_TXT_OR_TARGET 10001
#define ID_TXT_OR_WEIGHT 10002
#define ID_TXT_VARIANCE_TARGET 10003
#define ID_TXT_VARIANCE_WEIGHT 10004
#define ID_CMB_GA_SETTINGS 10103
#define ID_CMD_CONFIGURE_GA 10104
#define SYMBOL_WXPAGESIMPEN_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGESIMPEN_TITLE _("simPEN")
#define SYMBOL_WXPAGESIMPEN_IDNAME ID_WXPAGESIMPEN
#define SYMBOL_WXPAGESIMPEN_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGESIMPEN_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxPageSimPEN class declaration
 */

class wxPageSimPEN: public wxDiseaseConfiguration
{    
    DECLARE_DYNAMIC_CLASS( wxPageSimPEN )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageSimPEN();
    wxPageSimPEN( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGESIMPEN_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESIMPEN_POSITION, const wxSize& size = SYMBOL_WXPAGESIMPEN_SIZE, long style = SYMBOL_WXPAGESIMPEN_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGESIMPEN_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESIMPEN_POSITION, const wxSize& size = SYMBOL_WXPAGESIMPEN_SIZE, long style = SYMBOL_WXPAGESIMPEN_STYLE );

    /// Destructor
    ~wxPageSimPEN();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageSimPEN event handler declarations

    /// wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_GA_SETTINGS
    void OnCmbGaSettingsSelected( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIGURE_GA
    void OnCmdConfigureGaClick( wxCommandEvent& event );

////@end wxPageSimPEN event handler declarations

////@begin wxPageSimPEN member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageSimPEN member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageSimPEN member variables
    wxTextCtrl* txtPrevalence;
    wxTextCtrl* txtHeritTarget;
    wxTextCtrl* txtHeritWeight;
    wxTextCtrl* txtORTarget;
    wxTextCtrl* txtORWeight;
    wxTextCtrl* txtVarTarget;
    wxTextCtrl* txtVarianceWeight;
    wxChoice* cmbGASettings;
////@end wxPageSimPEN member variables

	std::string GetPreferredExtension() { return "simpen"; }
	virtual bool CanSave() { return true; }
	virtual bool CanSaveAs() { return true; }
	virtual bool CanImport() { return true; }
	virtual bool CanExport() { return false; }
	virtual bool CanEvaluate() { return false; }
	string GetImportFileTypes() { return "SIMpen Cfg(*.simpen)|*.simpen"; }


	virtual bool Import(const char *filename, string &desc);

	virtual bool Save(const char *filename, const char *description);
	virtual bool Evaluate() { assert(0); }

	int GetLocusCount() { return -1; }

	FileToMap fileContents;
	

};

}

}




#endif
    // _WXPAGESIMPEN_H_
