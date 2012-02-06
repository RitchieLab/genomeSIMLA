/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagepenetrancemodel.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 13:00:53 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGEPENETRANCEMODEL_H_
#define _WXPAGEPENETRANCEMODEL_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/spinctrl.h"
#include "wx/grid.h"
////@end includes
#include "wxdiseaseconfiguration.h"
#include <vector>
#include "utility/lineparser.h"

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxSpinCtrl;
class wxGrid;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXPAGEPENETRANCEMODEL 10096
#define ID_TXT_FREQ_THRESH 10113
#define ID_LOCUS_COUNT 10099
#define ID_GRD_ALLELE_FREQ 10098
#define ID_GRD_PENETRANCE 10101
#define SYMBOL_WXPAGEPENETRANCEMODEL_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGEPENETRANCEMODEL_TITLE _("Penetrance Model")
#define SYMBOL_WXPAGEPENETRANCEMODEL_IDNAME ID_WXPAGEPENETRANCEMODEL
#define SYMBOL_WXPAGEPENETRANCEMODEL_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGEPENETRANCEMODEL_POSITION wxDefaultPosition
////@end control identifiers

class wxPenGrid : public wxGridTableBase {
public:
	wxPenGrid(int modelSize);
	virtual int GetNumberRows();
	virtual int GetNumberCols();
	virtual bool IsEmptyCell(int row, int col);
	virtual wxString GetValue(int row, int col);

	virtual void SetValue(int row, int col, const wxString& value);
	virtual wxString GetRowLabelValue( int col );

    virtual wxString GetTypeName( int row, int col );
    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );

    virtual double GetValueAsDouble( int row, int col );
	virtual void SetValueAsDouble( int row, int col, double value );

	void SetLocusCount(int locusCount);

	void Refresh();
	void Reset(bool clearPenetrances = false);

	void BuildGenotypeLabels(uint modelSize);
	bool InvertAllele(int locusToInvert);
	void ConvertFromGenotypes(vector<int>& genotypes, int genotypeCount, int& mlGenotype);
	void ConvertToGenotypes(vector<int>& genotypes, int genotypeCount, int mlGenotype);
	
	
	bool Load(int lociCount, Utility::FileToMap& file);
	bool Save(const char *filename);
protected:
	vector<double> penetrances;
	int lociCount;
	vector<string> penLabels;


	string BuildGenotypeLabel(uint genotype, uint position);
	string BuildGenotypeLabel(uint *genotypes, uint modelSize);
};
class wxFreqGrid : public wxGridTableBase {
public:
	struct Freqs {
		float al1;
		float al2;
		Freqs() : al1(0.5), al2(0.5) { }
		Freqs(float al1, float al2) : al1(al1), al2(al2) { }
	};

	wxFreqGrid(int modelSize);
	virtual int GetNumberRows();
	virtual int GetNumberCols();
	virtual bool IsEmptyCell(int row, int col);
	virtual wxString GetValue(int row, int col);

	virtual void SetValue(int row, int col, const wxString& value);
	virtual wxString GetColLabelValue( int col );
	virtual wxString GetRowLabelValue( int row );
    virtual wxString GetTypeName( int row, int col );
    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );

    virtual double GetValueAsDouble( int row, int col );
	virtual void SetValueAsDouble( int row, int col, double value );
	
	void InvertFrequencies(int row);

	void SetLocusCount(size_t locusCount);
	size_t GetLocusCount();

	void Refresh();
	void Reset();
	
	void Load(vector<string>& freqs);


protected:
	vector<Freqs> frequencies;
	int lociCount;


};





/*!
 * wxPagePenetranceModel class declaration
 */

class wxPagePenetranceModel: public wxDiseaseConfiguration
{    
    DECLARE_DYNAMIC_CLASS( wxPagePenetranceModel )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPagePenetranceModel();
    wxPagePenetranceModel( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEPENETRANCEMODEL_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEPENETRANCEMODEL_POSITION, const wxSize& size = SYMBOL_WXPAGEPENETRANCEMODEL_SIZE, long style = SYMBOL_WXPAGEPENETRANCEMODEL_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEPENETRANCEMODEL_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEPENETRANCEMODEL_POSITION, const wxSize& size = SYMBOL_WXPAGEPENETRANCEMODEL_SIZE, long style = SYMBOL_WXPAGEPENETRANCEMODEL_STYLE );

    /// Destructor
    ~wxPagePenetranceModel();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPagePenetranceModel event handler declarations

    /// wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_LOCUS_COUNT
    void OnLocusCountUpdated( wxSpinEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_LOCUS_COUNT
    void OnLocusCountTextUpdated( wxCommandEvent& event );

    /// wxEVT_GRID_LABEL_RIGHT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnLabelRightClick( wxGridEvent& event );

    /// wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_ALLELE_FREQ
    void OnEditorShown( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnGrdAlleleFreqLabelLeftClick( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnGrdAlleleFreqLabelRightClick( wxGridEvent& event );

    /// wxEVT_SIZE event handler for ID_GRD_ALLELE_FREQ
    void OnSize( wxSizeEvent& event );

    /// wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_PENETRANCE
    void OnGrdPenetranceLabelLeftClick( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE
    void OnGrdPenetranceLabelRightClick( wxGridEvent& event );

////@end wxPagePenetranceModel event handler declarations

////@begin wxPagePenetranceModel member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPagePenetranceModel member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPagePenetranceModel member variables
    wxTextCtrl* txtFreqThreshold;
    wxSpinCtrl* dialLocusCount;
    wxGrid* grdAlleleFrequencies;
    wxGrid* grdPenetrances;
////@end wxPagePenetranceModel member variables

	void InitGrids();

	void ShowContextMenu(const wxPoint& pos, int snp);
	void InvertAlleles(wxCommandEvent& event);

	std::string GetPreferredExtension(){ return "pen";}
	virtual bool CanSave() { return !isLocked; }
	virtual bool CanSaveAs() { return true; }
	virtual bool CanImport() { return true; }
	virtual bool CanExport() { return true; }
	virtual bool CanEvaluate() { return true; }

	virtual bool Save(const char *filename, const char *description);
	virtual bool Import(const char *filename, string &desc);
	virtual bool Export(const char *filename, const char *desc) { assert(0); } 
	virtual bool Evaluate();
	string GetImportFileTypes() { return "Pen. Table(*.pen)|*.pen"; }
	int GetLocusCount() { return freqGridMaster->GetNumberRows(); }
	

protected:
	bool isLocked;					///<Indicates that the file can not be saved
	wxFreqGrid *freqGridMaster;
	wxPenGrid  *penGridMaster;
	int selectedSnp;
};

}

}
/*!
 * wxPagePenetranceModel class declaration
 */

class wxPagePenetranceModel: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( wxPagePenetranceModel )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPagePenetranceModel();
    wxPagePenetranceModel( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEPENETRANCEMODEL_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEPENETRANCEMODEL_POSITION, const wxSize& size = SYMBOL_WXPAGEPENETRANCEMODEL_SIZE, long style = SYMBOL_WXPAGEPENETRANCEMODEL_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGEPENETRANCEMODEL_IDNAME, const wxPoint& pos = SYMBOL_WXPAGEPENETRANCEMODEL_POSITION, const wxSize& size = SYMBOL_WXPAGEPENETRANCEMODEL_SIZE, long style = SYMBOL_WXPAGEPENETRANCEMODEL_STYLE );

    /// Destructor
    ~wxPagePenetranceModel();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPagePenetranceModel event handler declarations
    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL3
    void OnTextctrl3TextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_LOCUS_COUNT
    void OnLocusCountUpdated( wxSpinEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_LOCUS_COUNT
    void OnLocusCountTextUpdated( wxCommandEvent& event );

    /// wxEVT_GRID_LABEL_RIGHT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnLabelRightClick( wxGridEvent& event );

    /// wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_ALLELE_FREQ
    void OnEditorShown( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnGrdAlleleFreqLabelLeftClick( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_ALLELE_FREQ
    void OnGrdAlleleFreqLabelRightClick( wxGridEvent& event );

    /// wxEVT_SIZE event handler for ID_GRD_ALLELE_FREQ
    void OnSize( wxSizeEvent& event );

    /// wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_PENETRANCE
    void OnGrdPenetranceLabelLeftClick( wxGridEvent& event );

    /// wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE
    void OnGrdPenetranceLabelRightClick( wxGridEvent& event );

////@end wxPagePenetranceModel event handler declarations

////@begin wxPagePenetranceModel member function declarations
    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPagePenetranceModel member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPagePenetranceModel member variables
    wxTextCtrl* txtFreqThreshold;
    wxSpinCtrl* dialLocusCount;
    wxGrid* grdAlleleFrequencies;
    wxGrid* grdPenetrances;
////@end wxPagePenetranceModel member variables
};

#endif
    // _WXPAGEPENETRANCEMODEL_H_
