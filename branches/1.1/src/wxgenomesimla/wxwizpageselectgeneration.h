/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpageselectgeneration.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 18 Feb 2008 13:55:24 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGESELECTGENERATION_H_
#define _WXWIZPAGESELECTGENERATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes

#include "wxwizpageselectgeneration.h"
#include "wxpageselectgeneration.h"
#include "simwizard.h"

namespace GenomeSIM {

namespace GUI {
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZARDPAGE1 10078
#define ID_PANEL1 10079
#define ID_SKIP_ANALYSIS 10081
////@end control identifiers


/*!
 * wxWizPageSelectGeneration class declaration
 */

class wxWizPageSelectGeneration: public wxWizardPageSimple, public SimWizard 
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageSelectGeneration )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageSelectGeneration();

    wxWizPageSelectGeneration( wxWizard* parent, AppController *ctrl );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageSelectGeneration();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageSelectGeneration event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZARDPAGE1
    void OnWizardpage1PageChanged( wxWizardEvent& event );

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZARDPAGE1
    void OnWizardpage1PageChanging( wxWizardEvent& event );

////@end wxWizPageSelectGeneration event handler declarations

////@begin wxWizPageSelectGeneration member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageSelectGeneration member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageSelectGeneration member variables
    wxPageSelectGeneration* treeGenerationSelection;
////@end wxWizPageSelectGeneration member variables


};

}

}
#endif
    // _WXWIZPAGESELECTGENERATION_H_
