/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchselectgeneration.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:02 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGELAUNCHSELECTGENERATION_H_
#define _WXWIZPAGELAUNCHSELECTGENERATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/listctrl.h"
////@end includes
#include "simwizard.h"
#include "wxwizpagelaunchselectgeneration.h"
/*!
 * Forward declarations
 */


namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_SELECT_PREVIOUS_RUN 10062
#define ID_LST_SELECT_GENERATION 10063
////@end control identifiers


/*!
 * wxWizPageLaunchSelectGeneration class declaration
 */

class wxWizPageLaunchSelectGeneration: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageLaunchSelectGeneration )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageLaunchSelectGeneration();

    wxWizPageLaunchSelectGeneration( wxWizard* parent, AppController *ctrl );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageLaunchSelectGeneration();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageLaunchSelectGeneration event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_SELECT_PREVIOUS_RUN
    void OnWizSelectPreviousRunPageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_SELECT_PREVIOUS_RUN
    void OnWizSelectPreviousRunPageChanging( wxWizardEvent& event );

    /// wxEVT_WIZARD_FINISHED event handler for ID_WIZ_SELECT_PREVIOUS_RUN
    void OnWizSelectPreviousRunFinished( wxWizardEvent& event );

    /// wxEVT_SIZE event handler for ID_WIZ_SELECT_PREVIOUS_RUN
    void OnSize( wxSizeEvent& event );

////@end wxWizPageLaunchSelectGeneration event handler declarations

////@begin wxWizPageLaunchSelectGeneration member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageLaunchSelectGeneration member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageLaunchSelectGeneration member variables
    wxStaticText* lblDescription;
    wxListCtrl* listGenerations;
////@end wxWizPageLaunchSelectGeneration member variables
	void InitializeColumns();
	void AddEntry(ExecutionLog::LogEntry *entry = NULL);
};
}
}
#endif
    // _WXWIZPAGELAUNCHSELECTGENERATION_H_
