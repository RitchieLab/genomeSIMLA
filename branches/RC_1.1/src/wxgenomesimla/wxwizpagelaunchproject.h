/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchproject.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:12:50 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGELAUNCHPROJECT_H_
#define _WXWIZPAGELAUNCHPROJECT_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes
#include "simwizard.h"
#include "utility/executionlog.h"
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageLaunchProject;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_PROJECT_SETTINGS 10060
#define ID_CMB_PROJECT_NAME 10065
#define ID_TXT_PROJECT_PATH 10066
#define ID_CMD_SELECT_PROJECT_PATH 10067
////@end control identifiers


/*!
 * wxWizPageLaunchProject class declaration
 */

class wxWizPageLaunchProject: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageLaunchProject )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageLaunchProject();

    wxWizPageLaunchProject( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageLaunchProject();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageLaunchProject event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_PROJECT_SETTINGS
    void OnWizProjectSettingsPageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_PROJECT_SETTINGS
    void OnWizProjectSettingsPageChanging( wxWizardEvent& event );

    /// wxEVT_SIZE event handler for ID_WIZ_PROJECT_SETTINGS
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_COMBOBOX_SELECTED event handler for ID_CMB_PROJECT_NAME
    void OnCmbProjectNameSelected( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_SELECT_PROJECT_PATH
    void OnCmdSelectProjectPathClick( wxCommandEvent& event );

////@end wxWizPageLaunchProject event handler declarations

////@begin wxWizPageLaunchProject member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageLaunchProject member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageLaunchProject member variables
    wxComboBox* cmbProjectName;
    wxTextCtrl* txtProjectPath;
////@end wxWizPageLaunchProject member variables
	wxStaticText *lblDescription;
	wxWizardPage *GetNext() const;

	
	virtual void SetAppController(AppController *appController) {
		this->appController=appController;
		projectEntries = appController->GetProjectEntries();
	}

protected:
	ExecutionLog::RunType *projectEntries;
	vector<string> projectFilenames;
	void SetProjectDir(const char *filename);
};

}
}
#endif
    // _WXWIZPAGELAUNCHPROJECT_H_
