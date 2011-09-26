/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagedefinedatasets.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 11:10:58 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGEDEFINEDATASETS_H_
#define _WXWIZPAGEDEFINEDATASETS_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
#include "wx/listctrl.h"
////@end includes
#include <string>
#include <vector>
#include "modellist.h"
#include "wxwizardpagedatasim.h"
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageDefineDatasets;
class wxListCtrl;
////@end forward declarations
namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PAGE_DEFINE_DATASETS 10124
#define ID_TXT_DATASET_COUNT 10000
#define ID_LST_DATASETS 10002
#define ID_CHK_STD_HEADER 10003
#define ID_CHK_BINARY_DATASETS 10004
#define ID_CMD_ADD_CC 10005
#define ID_CMD_ADD_PED 10006
#define ID_DEL_DATASET 10015
#define ID_CMB_MODEL 10016
#define ID_CMD_CONFIG_MODEL 10017
////@end control identifiers


/*!
 * wxWizPageDefineDatasets class declaration
 */

class wxWizPageDefineDatasets: public wxWizardPageDataSim
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageDefineDatasets )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageDefineDatasets();

    wxWizPageDefineDatasets( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~wxWizPageDefineDatasets();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageDefineDatasets event handler declarations

    /// wxEVT_WIZARD_PAGE_CHANGING event handler for ID_PAGE_DEFINE_DATASETS
    void OnPageDefineDatasetsPageChanging( wxWizardEvent& event );

    /// wxEVT_SIZE event handler for ID_PAGE_DEFINE_DATASETS
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_DATASETS
    void OnLstDatasetsSelected( wxListEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_DATASETS
    void OnLstDatasetsDeselected( wxListEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LST_DATASETS
    void OnLstDatasetsItemActivated( wxListEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC
    void OnCmdAddCcClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED
    void OnCmdAddPedClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_DEL_DATASET
    void OnDelDatasetClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_MODEL
    void OnCmbModelSelected( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL
    void OnCmdConfigModelClick( wxCommandEvent& event );

////@end wxWizPageDefineDatasets event handler declarations

////@begin wxWizPageDefineDatasets member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageDefineDatasets member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageDefineDatasets member variables
    wxTextCtrl* txtDatasetCount;
    wxListCtrl* lstDatasets;
    wxCheckBox* chkStandardPedigreeHeader;
    wxCheckBox* chkWriteBinary;
    wxChoice* lstDiseaseModels;
////@end wxWizPageDefineDatasets member variables
	virtual bool TransferDataFromWindow();
	
	void RefreshSize();
	void InitListBox();
	/**
	 * @brief creates an array of the filenames contained within the file, modelList
	 */
	void InitializeDiseaseModels();	

	void PrepSummary(wxWizardPageDataSummary *summary);
	void PrepConfig(ostream& config);


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

protected:
	long selectedModel;

	struct DatasetDetails {
		string label;
		string details;

		DatasetDetails() : label(""), details("") { }
		DatasetDetails(const char *label, const char *details) : label(label), details(details) { }
	};
	vector<DatasetDetails> datasets;
	ModelList *modelLibrary;
};

}

}

#endif
    // _WXWIZPAGEDEFINEDATASETS_H_
