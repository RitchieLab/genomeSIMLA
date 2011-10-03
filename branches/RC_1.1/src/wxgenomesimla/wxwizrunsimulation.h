/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizrunsimulation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 29 Jan 2008 11:51:24 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZRUNSIMULATION_H_
#define _WXWIZRUNSIMULATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes
#include "appcontroller.h"
#include "simwizard.h"

#include "wxwizpagelaunchproject.h"
#include "wxwizpagelaunchdroppoints.h"
#include "wxwizpagelaunchselectgeneration.h"
#include "wxwizpagelaunchconfirmation.h"
#include "wxwizpagelaunchsimulation.h"
#include "wxwizpagereview.h"



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
#define ID_WXWIZRUNSIMULATION 10061
#define SYMBOL_WXWIZRUNSIMULATION_IDNAME ID_WXWIZRUNSIMULATION
////@end control identifiers





class wxWizRunSimulation: public wxWizard, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizRunSimulation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizRunSimulation();
    wxWizRunSimulation( wxWindow* parent, bool skipGenSelection, AppController *appCtrl, wxWindowID id = SYMBOL_WXWIZRUNSIMULATION_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Creation
    bool Create( wxWindow* parent, bool skipGenSelection, wxWindowID id = SYMBOL_WXWIZRUNSIMULATION_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Destructor
    ~wxWizRunSimulation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls(bool skipGenSelection);

////@begin wxWizRunSimulation event handler declarations

////@end wxWizRunSimulation event handler declarations

////@begin wxWizRunSimulation member function declarations

    /// Runs the wizard
    bool Run();

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizRunSimulation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();


////@begin wxWizRunSimulation member variables
    wxWizPageLaunchProject* pageProjectSettings;
//    wxWizPageLaunchSelectGeneration* pageSelectPreviousRun;
    wxWizPageLaunchDropPoints* pageDropPoints;
    wxWizPageLaunchConfirmation* pageConfirmation;
    wxWizPageReview* pageReview;
////@end wxWizRunSimulation member variables

//	void Initialize();

};


}
}


#endif
    // _WXWIZRUNSIMULATION_H_
