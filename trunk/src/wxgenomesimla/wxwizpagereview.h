/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagereview.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 04 Feb 2008 16:58:16 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGEREVIEW_H_
#define _WXWIZPAGEREVIEW_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/treectrl.h"
////@end includes
#include "utility/executionlog.h"
#include "simwizard.h"

#ifdef WIN32
#define FILE_PREDICATE "file:///"
#else
#define FILE_PREDICATE "file://"
#endif

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageReview;
class wxTreeListCtrl;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace Utility;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WIZ_REVIEW_RESULTS 10068
#define ID_TREECTRL1 10069
////@end control identifiers



/*!
 * wxWizPageReview class declaration
 */

class wxWizPageReview: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageReview )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageReview();

    wxWizPageReview( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageReview();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageReview event handler declarations

    /// wxEVT_COMMAND_TREE_ITEM_ACTIVATED event handler for ID_TREECTRL1
    void OnTreectrl1ItemActivated( wxTreeEvent& event );

////@end wxWizPageReview event handler declarations

////@begin wxWizPageReview member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageReview member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageReview member variables
    wxTreeListCtrl* reportTree;
////@end wxWizPageReview member variables

	void LoadResults();

	void InitSimResults();

	void AddNewEntry(ExecutionLog::LogEntry &entry);

	void ClearReports();
	
protected:
	bool doLoadReports;
	wxTreeItemId root;	
};

}

}


#endif
    // _WXWIZPAGEREVIEW_H_
