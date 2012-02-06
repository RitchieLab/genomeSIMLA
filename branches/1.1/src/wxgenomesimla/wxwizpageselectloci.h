/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpageselectloci.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 11:10:40 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGESELECTLOCI_H_
#define _WXWIZPAGESELECTLOCI_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/grid.h"
#include "wx/spinctrl.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageSelectLoci;
class wxGrid;
class wxSpinCtrl;
////@end forward declarations

#include <vector>
#include "simulation/locus.h"
#include "locusmanager.h"
#include <map>
#include "wxwizardpagedatasim.h"
#include "diseasemodeldetails.h"

namespace GenomeSIM {

namespace GUI {

using namespace Simulation;
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PAGE_SELECT_LOCI 10125
#define ID_GRD_LOCI 10126
#define ID_SPINCTRL_LOCUSCOUNT 10007
////@end control identifiers


/*!
 * wxWizPageSelectLoci class declaration
 */

class wxWizPageSelectLoci: public wxWizardPageDataSim
{
    DECLARE_DYNAMIC_CLASS( wxWizPageSelectLoci )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageSelectLoci();

    wxWizPageSelectLoci( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageSelectLoci();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageSelectLoci event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_PAGE_SELECT_LOCI
    void OnPageSelectLociPageChanging( wxWizardEvent& event );

    /// wxEVT_SIZE event handler for ID_PAGE_SELECT_LOCI
    void OnSize( wxSizeEvent& event );

    /// wxEVT_GRID_CELL_LEFT_DCLICK event handler for ID_GRD_LOCI
    void OnLeftDClick( wxGridEvent& event );

    /// wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_SPINCTRL
    void OnSpinctrlUpdated( wxSpinEvent& event );

////@end wxWizPageSelectLoci event handler declarations

////@begin wxWizPageSelectLoci member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageSelectLoci member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageSelectLoci member variables
    wxGrid* grdLocusSelection;
    wxSpinCtrl* spnLocusCount;
////@end wxWizPageSelectLoci member variables

	void RefreshSize();
	//void InitListBox();
	void SetLocusCount(int locusCount);
	void SetLocus(int idx, Locus *locus);
	void SetModel(ModelList::ModelFileItem &item);
	bool Initialize(vector<FileBasedChromosome *> *chrom);
	void PrepSummary(wxWizardPageDataSummary *summary);
	void PrepConfig(ostream& config);
	bool TransferDataFromWindow();

	/**
	 * @brief This is called prior to saving the configuration
	 */
	void Commit();

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings();

	/**
	 * @brief Compares current state with that of the saved version at least against what is in memory
	 * @return True indicates that something has changed 
	 */
	bool HasChanged();

	static void *SummarizeDiseaseModel(void *argv);

public:
	map<string, FileBasedChromosome *> chromosomes;
	wxString *labels;							///<Seems to be required for the combo box
	DiseaseModelDetails *modelDetails;
	bool directionOfTransition;					///<State information to avoid running TransferDataFromWindow when going back

};

}

}
#endif
    // _WXWIZPAGESELECTLOCI_H_
