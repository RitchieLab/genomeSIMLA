/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgconfigurepeddatasets.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 21 Mar 2008 12:39:38 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGCONFIGUREPEDDATASETS_H_
#define _WXDLGCONFIGUREPEDDATASETS_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/grid.h"
#include "wx/toolbar.h"
////@end includes
#include "wxuser.h"
#include <string>
#include <vector>
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxGrid;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace std;


/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGCONFIGUREPEDDATASETS 10116
#define ID_DISEASE_LABEL 10004
#define ID_DISEASE_GT_ERROR 10005
#define ID_DISEASE_PHENO_ERR 10006
#define ID_DISEASE_MISSING 10007
#define ID_FAMTYPE_GRID 10003
#define ID_BTN_ADD_FAMILY_TYPE 10000
#define ID_BTN_DEL_FAMILY_TYPE 10001
#define SYMBOL_WXDLGCONFIGUREPEDDATASETS_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGCONFIGUREPEDDATASETS_TITLE _("Configure Pedigree Datasets")
#define SYMBOL_WXDLGCONFIGUREPEDDATASETS_IDNAME ID_WXDLGCONFIGUREPEDDATASETS
#define SYMBOL_WXDLGCONFIGUREPEDDATASETS_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGCONFIGUREPEDDATASETS_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDlgConfigurePedDatasets class declaration
 */
class wxDlgConfigurePedDatasets: public wxDialog, public WxUser
{    
    DECLARE_DYNAMIC_CLASS( wxDlgConfigurePedDatasets )
    DECLARE_EVENT_TABLE()

public:

struct FamType {
	int affected;
	int unaffected;
	int extrasibs;
	int famcount;

	FamType(int affected, int unaffected, int extrasibs, int famcount) : 
			affected(affected), unaffected(unaffected), extrasibs(extrasibs), famcount(famcount) { }
	FamType() : affected(1), unaffected(0), extrasibs(0), famcount(100) { }
	FamType(const FamType& other) : affected(other.affected), unaffected(other.unaffected), 
		extrasibs(other.extrasibs), famcount(other.famcount) { }
	string GetDetails();
	bool Parse(const char *line);
};

class wxGridFamilyTypes : public wxGridTableBase {
public:
	wxGridFamilyTypes() {}
	~wxGridFamilyTypes() {}
	virtual int GetNumberRows() { return familyTypes.size(); }
	/**
	 * @brief aff, unaff, extra, fam. count
	 */
	virtual int GetNumberCols() { return 4; } 

	virtual bool IsEmptyCell(int row, int col);
	virtual void SetValue(int row, int col, const wxString& value);
	virtual wxString GetValue(int row, int col);

	virtual long GetValueAsLong(int row, int col);
	virtual void SetValueAsLong(int row, int col, long value);

    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );

	virtual wxString GetColLabelValue( int col );
 	virtual wxString GetRowLabelValue( int row );
	virtual wxString GetTypeName( int row, int col );

	void AddFamilyType();
	void AddFamilyType(const char *details);
	void DelFamilyType(int idx);
	void Refresh();

	string GetDetails();
	void ConfigureModelDetails(istream &modelDetails);
	string GetSummary();

	vector<FamType> familyTypes;
};



    /// Constructors
    wxDlgConfigurePedDatasets();
    wxDlgConfigurePedDatasets( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGUREPEDDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGUREPEDDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGUREPEDDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGUREPEDDATASETS_SIZE, long style = SYMBOL_WXDLGCONFIGUREPEDDATASETS_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGCONFIGUREPEDDATASETS_IDNAME, const wxString& caption = SYMBOL_WXDLGCONFIGUREPEDDATASETS_TITLE, const wxPoint& pos = SYMBOL_WXDLGCONFIGUREPEDDATASETS_POSITION, const wxSize& size = SYMBOL_WXDLGCONFIGUREPEDDATASETS_SIZE, long style = SYMBOL_WXDLGCONFIGUREPEDDATASETS_STYLE );

    /// Destructor
    ~wxDlgConfigurePedDatasets();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgConfigurePedDatasets event handler declarations

    /// wxEVT_SIZE event handler for ID_WXDLGCONFIGUREPEDDATASETS
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_BTN_ADD_FAMILY_TYPE
    void OnBtnAddFamilyTypeClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_BTN_DEL_FAMILY_TYPE
    void OnBtnDelFamilyTypeClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CANCEL
    void OnCancelClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
    void OnOkClick( wxCommandEvent& event );

////@end wxDlgConfigurePedDatasets event handler declarations

////@begin wxDlgConfigurePedDatasets member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgConfigurePedDatasets member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

	string GetModelDetails();

////@begin wxDlgConfigurePedDatasets member variables
    wxTextCtrl* txtLabel;
    wxTextCtrl* txtGenotypeErr;
    wxTextCtrl* txtPhenocopy;
    wxTextCtrl* txtMissingData;
    wxGrid* grdFamTypes;
////@end wxDlgConfigurePedDatasets member variables
	string GetModelLabel();
	int GetFamilyTypes() { return gridMaster->GetNumberRows(); }
	FamType GetFamilyType(int idx) { assert(idx < gridMaster->GetNumberRows()); return gridMaster->familyTypes[idx]; }
	string GetModelSummary();
	void ConfigureModelDetails(const char *details);
protected:
	wxGridFamilyTypes *gridMaster;
};



}

}

#endif
    // _WXDLGCONFIGUREPEDDATASETS_H_
