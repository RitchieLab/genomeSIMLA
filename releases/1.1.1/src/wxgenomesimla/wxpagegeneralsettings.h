/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagegeneralsettings.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 20 Dec 2007 15:26:39 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGEGENERALSETTINGS_H_
#define _WXPAGEGENERALSETTINGS_H_


/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxImgGrowthChart;
////@end forward declarations
#include <string>
#include "appinterface.h"
#include "wximggrowthchart.h"
#include "wx/wizard.h"

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXPAGEGENERALSETTINGS 10032
#define ID_TXT_RANDOM_SEED 10000
#define ID_TXT_SIMULTANEOUS_CHROMS 10047
#define IS_TXT_THR_PER_CHROM 10050
#define ID_IMGGROWTHCHART 10004
#define ID_TXT_GROWTH_CFG 10005
#define ID_CMD_CONFIG_GROWTH 10006
#define SYMBOL_WXPAGEGENERALSETTINGS_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGEGENERALSETTINGS_TITLE _("General Settings")
#define SYMBOL_WXPAGEGENERALSETTINGS_IDNAME ID_WXPAGEGENERALSETTINGS
#define SYMBOL_WXPAGEGENERALSETTINGS_SIZE wxSize(550, 550)
#define SYMBOL_WXPAGEGENERALSETTINGS_POSITION wxDefaultPosition
////@end control identifiers

namespace GenomeSIM {

namespace GUI {
/*!
 * wxPageGeneralSettings class declaration
 */

class wxPageGeneralSettings: public wxPanel, public AppInterface 	{    
    DECLARE_DYNAMIC_CLASS( wxPageGeneralSettings )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageGeneralSettings();
    wxPageGeneralSettings( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEGENERALSETTINGS_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEGENERALSETTINGS_POSITION, const wxSize& size = SYMBOL_WXPAGEGENERALSETTINGS_SIZE, long style = SYMBOL_WXPAGEGENERALSETTINGS_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEGENERALSETTINGS_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEGENERALSETTINGS_POSITION, const wxSize& size = SYMBOL_WXPAGEGENERALSETTINGS_SIZE, long style = SYMBOL_WXPAGEGENERALSETTINGS_STYLE );

    /// Destructor
    ~wxPageGeneralSettings();

    /// Initialises member variables
    void Init();

	void InitAppController(AppController *ctrl) { appController = ctrl; }

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageGeneralSettings event handler declarations

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_RANDOM_SEED
    void OnTxtRandomSeedUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_RANDOM_SEED
    void OnTxtRandomSeedEnter( wxCommandEvent& event );

    /// wxEVT_SIZE event handler for ID_IMGGROWTHCHART
    void OnSize( wxSizeEvent& event );

    /// wxEVT_LEFT_DOWN event handler for ID_IMGGROWTHCHART
    void OnLeftDown( wxMouseEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_GROWTH
    void OnCmdConfigGrowthClick( wxCommandEvent& event );

////@end wxPageGeneralSettings event handler declarations

////@begin wxPageGeneralSettings member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageGeneralSettings member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageGeneralSettings member variables
    wxTextCtrl* txtRandomSeed;
    wxTextCtrl* txtSimultaneousChrom;
    wxTextCtrl* txtThreadsPerChrom;
    wxImgGrowthChart* imgGrowthChart;
    wxTextCtrl* txtGrowthParams;
////@end wxPageGeneralSettings member variables

	void Init(AppController *appController);

	
	/**
	 * @brief This is called prior to saving the configuration
	 */
	void Commit();

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings();

	bool HasChanged();

	void OnWizConfirmationPageChanged( wxWizardEvent& event );
	bool VerifyForRun() { return true; }
private:
		bool hasChanged;
};

}

}

#endif
    // _WXPAGEGENERALSETTINGS_H_
