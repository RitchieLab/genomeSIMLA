/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchdroppoints.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:13:11 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGELAUNCHDROPPOINTS_H_
#define _WXWIZPAGELAUNCHDROPPOINTS_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/statline.h"
////@end includes
#include "simwizard.h"
namespace GenomeSIM {

namespace GUI {

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageLaunchDropPoints;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_DROP_POINTS 10064
#define ID_TXT_DROP_FIRST 10000
#define ID_TXT_DROP_FREQ 10004
#define ID_TXT_DROP_COUNT 10005
////@end control identifiers


/*!
 * wxWizPageLaunchDropPoints class declaration
 */

class wxWizPageLaunchDropPoints: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageLaunchDropPoints )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageLaunchDropPoints();

    wxWizPageLaunchDropPoints( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageLaunchDropPoints();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageLaunchDropPoints event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_DROP_POINTS
    void OnWizDropPointsPageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_DROP_POINTS
    void OnWizDropPointsPageChanging( wxWizardEvent& event );

    /// wxEVT_WIZARD_CANCEL event handler for ID_WIZ_DROP_POINTS
    void OnWizDropPointsCancel( wxWizardEvent& event );

    /// wxEVT_WIZARD_FINISHED event handler for ID_WIZ_DROP_POINTS
    void OnWizDropPointsFinished( wxWizardEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FIRST
    void OnTxtDropFirstUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FIRST
    void OnTxtDropFirstEnter( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_FREQ
    void OnTxtDropFreqUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_FREQ
    void OnTxtDropFreqEnter( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DROP_COUNT
    void OnTxtDropCountUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_DROP_COUNT
    void OnTxtDropCountEnter( wxCommandEvent& event );

////@end wxWizPageLaunchDropPoints event handler declarations

////@begin wxWizPageLaunchDropPoints member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageLaunchDropPoints member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageLaunchDropPoints member variables
    wxTextCtrl* txtDropInitial;
    wxTextCtrl* txtDropFrequency;
    wxTextCtrl* txtDropCount;
////@end wxWizPageLaunchDropPoints member variables
	
	
	wxWizardPage *GetPrev() const;
};

}
}
#endif
    // _WXWIZPAGELAUNCHDROPPOINTS_H_
