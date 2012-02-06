/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizperformanalysis.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 27 Feb 2008 15:42:32 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPERFORMANALYSIS_H_
#define _WXWIZPERFORMANALYSIS_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageSelectGeneration;
class wxPageSelectGeneration;
class wxWizPageMonitorAnalysis;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXWIZPERFORMANALYSIS 10052
#define SYMBOL_WXWIZPERFORMANALYSIS_IDNAME ID_WXWIZPERFORMANALYSIS
////@end control identifiers


/*!
 * wxWizPerformAnalysis class declaration
 */

class wxWizPerformAnalysis: public wxWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPerformAnalysis )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPerformAnalysis();
    wxWizPerformAnalysis( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZPERFORMANALYSIS_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZPERFORMANALYSIS_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Destructor
    ~wxWizPerformAnalysis();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPerformAnalysis event handler declarations
////@end wxWizPerformAnalysis event handler declarations

////@begin wxWizPerformAnalysis member function declarations
    /// Runs the wizard
    bool Run();

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPerformAnalysis member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPerformAnalysis member variables
    wxWizPageSelectGeneration* pageSelectGeneration;
    wxWizPageMonitorAnalysis* pageMonitor;
////@end wxWizPerformAnalysis member variables
};

#endif
    // _WXWIZPERFORMANALYSIS_H_
