/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigurediseasemodel.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 12:56:57 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGCONFIGUREDISEASEMODEL_H_
#define _WXDLGCONFIGUREDISEASEMODEL_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/notebook.h"
////@end includes
#include "wxpagepenetrancemodel.h"
#include "wxpagesimpen.h"
#include "wxdlgsimla.h"
#include "simulation/locus.h"

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxNotebook;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace Simulation;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGCONFIGUREDISEASEMODEL 10094
#define ID_TEXTCTRL_DESCRIPTION 10001
#define ID_NB_DiseaseCfg 10095
#define ID_CMD_EVALUATE 10121
#define SYMBOL_WXDLGCONFIGUREDISEASEMODEL_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGCONFIGUREDISEASEMODEL_TITLE _("Configure Disease Model")
#define SYMBOL_WXDLGCONFIGUREDISEASEMODEL_IDNAME ID_WXDLGCONFIGUREDISEASEMODEL
#define SYMBOL_WXDLGCONFIGUREDISEASEMODEL_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGCONFIGUREDISEASEMODEL_POSITION wxDefaultPosition
////@end control identifiers

class wxGridDiseaseAssignment : public wxGridTableBase {
public:
	wxGridDiseaseAssignment() {}
	~wxGridDiseaseAssignment() {}
	virtual int GetNumberRows() { return locusCount; }
	/**
	 * @brief Label, freq(1), freq(2), Details
	 */
	virtual int GetNumberCols() { return 4; } 

	virtual bool IsEmptyCell(int row, int col);
	virtual wxString GetValue(int row, int col);
	virtual void SetValue(int row, int col, const wxString& value);

	virtual wxString GetColLabelValue( int col );
 	virtual wxString GetRowLabelValue( int row );
	virtual wxString GetTypeName( int row, int col );

	void SetLocusCount(int locCount);
	/**
	 * @brief Set the locus at index, idx. 
	 */
	void SetLocus(int idx, Locus& loc, const char *details);
protected:
	int locusCount;

	struct LocusEntry {
		Locus l;
		string details;
		LocusEntry() {}
		LocusEntry(Locus& l, const char *details) : l(l), details(details) { }
	};
	vector<LocusEntry> loci;

};

/*!
 * wxDlgConfigureDiseaseModel class declaration
 */

class wxDlgConfigureDiseaseModel: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgConfigureDiseaseModel )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgConfigureDiseaseModel();
    wxDlgConfigureDiseaseModel( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_SIZE, long style = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_SIZE, long style = SYMBOL_WXDLGCONFIGUREDISEASEMODEL_STYLE );

    /// Destructor
    ~wxDlgConfigureDiseaseModel();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgConfigureDiseaseModel event handler declarations

    /// wxEVT_COMMAND_NOTEBOOK_PAGE_CHANGED event handler for ID_NB_DiseaseCfg
    void OnNBDiseaseCfgPageChanged( wxNotebookEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_EVALUATE
    void OnCmdEvaluateClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVE
    void OnSaveClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVEAS
    void OnSaveasClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OPEN
    void OnOpenClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
    void OnCancelClick( wxCommandEvent& event );

////@end wxDlgConfigureDiseaseModel event handler declarations

////@begin wxDlgConfigureDiseaseModel member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgConfigureDiseaseModel member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgConfigureDiseaseModel member variables
    wxTextCtrl* txtDescription;
    wxNotebook* nbDiseaseConfiguration;
    wxButton* cmdEvaluate;
    wxButton* cmdSave;
    wxButton* cmdSaveAs;
    wxButton* cmdOpen;
    wxButton* cmdCancel;
////@end wxDlgConfigureDiseaseModel member variables


	bool LoadModel(const char *filetype, const char *filename, bool editable);

	int GetLocusCount();

	wxPagePenetranceModel *pagePenetrance;
	wxPageSimPEN *pageSimpen;
	wxDlgSIMLA *pageSimla;
	bool newFile;
	string filename;
	string filetype;
	string GetLabel();
	wxDiseaseConfiguration *cfg;
};

}

}

#endif
    // _WXDLGCONFIGUREDISEASEMODEL_H_
