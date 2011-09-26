/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchconfirmation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:19 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGELAUNCHCONFIRMATION_H_
#define _WXWIZPAGELAUNCHCONFIRMATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/richtext/richtextctrl.h"
////@end includes
#include "wxwizpagelaunchconfirmation.h"
#include "simwizard.h"
/*!
 * Forward declarations
 */

class wxRichTextCtrl;

namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_CONFIRMATION 10059
#define ID_TXT_SIM_DETAILS 10071
////@end control identifiers


/*!
 * wxWizPageLaunchConfirmation class declaration
 */

class wxWizPageLaunchConfirmation: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageLaunchConfirmation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageLaunchConfirmation();

    wxWizPageLaunchConfirmation( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageLaunchConfirmation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageLaunchConfirmation event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_CONFIRMATION
    void OnWizConfirmationPageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_CONFIRMATION
    void OnWizConfirmationPageChanging( wxWizardEvent& event );

////@end wxWizPageLaunchConfirmation event handler declarations

////@begin wxWizPageLaunchConfirmation member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageLaunchConfirmation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageLaunchConfirmation member variables
    wxRichTextCtrl* txtSimDetails;
////@end wxWizPageLaunchConfirmation member variables


	

protected:
	bool isComplete;
	
};
}
}
#endif
    // _WXWIZPAGELAUNCHCONFIRMATION_H_
