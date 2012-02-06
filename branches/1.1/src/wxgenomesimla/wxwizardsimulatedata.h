/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizardsimulatedata.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 10:39:25 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZARDSIMULATEDATA_H_
#define _WXWIZARDSIMULATEDATA_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes

#include <vector>
#include "modellist.h"

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageDefineDatasets;
class wxListCtrl;
class wxWizPageSelectLoci;
class wxGrid;
class wxSpinCtrl;
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
#define ID_WXWIZARDSIMULATEDATA 10119
#define SYMBOL_WXWIZARDSIMULATEDATA_IDNAME ID_WXWIZARDSIMULATEDATA
////@end control identifiers


/*!
 * wxWizardSimulateData class declaration
 */

class wxWizardSimulateData: public wxWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizardSimulateData )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizardSimulateData();
    wxWizardSimulateData( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZARDSIMULATEDATA_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZARDSIMULATEDATA_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Destructor
    ~wxWizardSimulateData();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizardSimulateData event handler declarations

////@end wxWizardSimulateData event handler declarations

////@begin wxWizardSimulateData member function declarations

    /// Runs the wizard
    bool Run();

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizardSimulateData member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizardSimulateData member variables
    wxWizPageDefineDatasets* pageDefineDatasets;
    wxWizPageSelectLoci* pageSelectLoci;
    wxWizPageReviewSimulation* pageSimulationReview;
////@end wxWizardSimulateData member variables
	bool SetChromosomes(vector<FileBasedChromosome *> *ch);
	void InitAppController(AppController *appController);
protected:
	vector<FileBasedChromosome *> *chromosomes;
};

}

}


#endif
    // _WXWIZARDSIMULATEDATA_H_
