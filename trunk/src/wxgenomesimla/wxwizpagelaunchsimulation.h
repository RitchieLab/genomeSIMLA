/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchsimulation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:26 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGELAUNCHSIMULATION_H_
#define _WXWIZPAGELAUNCHSIMULATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes
#include "simwizard.h"
#include "wxwizpagelaunchsimulation.h"
/*!
 * Forward declarations
 */


namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_MONITOR_RUN 10068
#define ID_TXT_SIM_LOG 10069
#define ID_SIM_COMPLETION 10070
////@end control identifiers


/*!
 * wxWizPageLaunchSimulation class declaration
 */

class wxWizPageLaunchSimulation: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageLaunchSimulation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageLaunchSimulation();

    wxWizPageLaunchSimulation( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageLaunchSimulation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageLaunchSimulation event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_MONITOR_RUN
    void OnWizMonitorRunPageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_MONITOR_RUN
    void OnWizMonitorRunPageChanging( wxWizardEvent& event );

    /// wxEVT_WIZARD_CANCEL event handler for ID_WIZ_MONITOR_RUN
    void OnWizMonitorRunCancel( wxWizardEvent& event );

////@end wxWizPageLaunchSimulation event handler declarations

////@begin wxWizPageLaunchSimulation member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageLaunchSimulation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageLaunchSimulation member variables
    wxTextCtrl* txtSimOutput;
    wxGauge* guageSimCompletion;
////@end wxWizPageLaunchSimulation member variables

	virtual bool TransferDataFromWindow() { return isComplete; }

	bool isComplete;
};

}

}

#endif
    // _WXWIZPAGELAUNCHSIMULATION_H_
