/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagereviewsimulation.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 28 Mar 2008 11:10:59 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGEREVIEWSIMULATION_H_
#define _WXWIZPAGEREVIEWSIMULATION_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/html/htmlwin.h"
////@end includes
#include "wxwizardpagedatasim.h"
#include <sstream>
#include <iostream>
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageReviewSimulation;
class wxHtmlWindow;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PAGE_SIMULATION_REVIEW 10127
#define ID_HTMLWINDOW2 10130
#define ID_CMD_GENERATE_DATA 10129
////@end control identifiers
#define ID_SAVE_PENETRANCE 30129

/*!
 * wxWizPageReviewSimulation class declaration
 */

class wxWizPageReviewSimulation: public wxWizardPageDataSummary
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageReviewSimulation )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageReviewSimulation();

    wxWizPageReviewSimulation( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageReviewSimulation();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageReviewSimulation event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGED event handler for ID_PAGE_SIMULATION_REVIEW
    void OnPageSimulationReviewPageChanged( wxWizardEvent& event );

    /// wxEVT_SIZE event handler for ID_PAGE_SIMULATION_REVIEW
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_GENERATE_DATA
    void OnCmdGenerateDataClick( wxCommandEvent& event );

////@end wxWizPageReviewSimulation event handler declarations

////@begin wxWizPageReviewSimulation member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageReviewSimulation member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();


////@begin wxWizPageReviewSimulation member variables
    wxHtmlWindow* txtReview;
////@end wxWizPageReviewSimulation member variables
	wxButton *wxSavePenetrance;
	void Clear(const char *txt);
	void AddHeader(const char *txt);
	void AddNote(const char *txt);
	void AddConfigLine(const char *txt);
	
	void OnSavePenetrance( wxCommandEvent& event);

	/**
	 * @brief This is called prior to saving the configuration
	 */
	void Commit() { }

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings() { }

	/**
	 * @brief Compares current state with that of the saved version at least against what is in memory
	 * @return True indicates that something has changed 
	 */
	bool HasChanged() { return false;} 
protected:
	std::stringstream configuration;

};

}

}

#endif
    // _WXWIZPAGEREVIEWSIMULATION_H_
