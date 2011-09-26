/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsimla.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 15:55:49 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGSIMLA_H_
#define _WXDLGSIMLA_H_

#include <string>
#include <vector>
#include "wxdiseaseconfiguration.h"
#include "utility/lineparser.h"

/*!
 * Includes
 */

////@begin includes
#include "wx/spinctrl.h"
#include "wx/grid.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxSpinCtrl;
class wxGrid;
////@end forward declarations
#include <map>
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGSIMLA 10105
#define ID_SPINCTRL 10001
#define ID_TARGETPREV 10000
#define ID_GRID 10002
#define ID_GRID1 10003
#define SYMBOL_WXDLGSIMLA_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGSIMLA_TITLE _("Configure SIMLA")
#define SYMBOL_WXDLGSIMLA_IDNAME ID_WXDLGSIMLA
#define SYMBOL_WXDLGSIMLA_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGSIMLA_POSITION wxDefaultPosition
////@end control identifiers

#define SIMLA_LOCUS "SIMLA_LOCUS"
#define SIMLA_INTERACTION "SIMLA_INTERACTION"
#define SIMLA_PREVALENCE "SIMLA_PREVALENCE"

namespace GenomeSIM {

namespace GUI {

using namespace std;
using namespace Utility;


class wxGridBetaWise : public wxGridTableBase {
public:
	double ToBeta(double value) { 		assert(value > 0.0); return log(value);	}
	double FromBeta(double value) { 	return exp(value); }
	wxGridBetaWise(int locCount = 0) : locusCount(locCount) { }
	void SetLocusCount(int locusCount);
	int GetLocusCount();
	virtual void ResizeContents(int locusCount) = 0;

protected:
	int locusCount;
};

class wxGridSimlaLoci : public wxGridBetaWise {
public:
struct LocusDetails {
	string ID;
	bool diseaseAtMajor;
	double beta;
	double modelType;
	LocusDetails() :ID(""), diseaseAtMajor(false), beta(1.0), modelType(0.0) { }
	LocusDetails(string id, bool atMajor, double& beta, double& modelType) : 
		ID(id), diseaseAtMajor(atMajor), beta(beta), modelType(modelType) { }
};
	wxGridSimlaLoci(int locCount) : wxGridBetaWise(locCount) {  ResizeContents(locCount);}
	virtual int GetNumberRows();
	virtual int GetNumberCols();
	virtual bool IsEmptyCell(int row, int col);
	virtual wxString GetValue(int row, int col);

	virtual void SetValue(int row, int col, const wxString& value);
	virtual wxString GetRowLabelValue( int col );
	virtual wxString GetColLabelValue( int col );

    virtual wxString GetTypeName( int row, int col );
    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );
	virtual bool GetValueAsBool( int row, int col);
	virtual void SetValueAsBool (int row, int col, bool isTrue);
    virtual double GetValueAsDouble( int row, int col );
	virtual void SetValueAsDouble( int row, int col, double value );

	void ResizeContents(int locusCount);

	void Refresh();
	void Reset();
	void InitCheckBoxes();

	void Load(FileToMap& file);
	void Save(ostream& os);
	bool Verify();
protected:
	vector<LocusDetails> loci;
};


class wxGridSimlaInteractions : public wxGridBetaWise {
public:
struct Interactions {
	string ID;
	double beta;
	Interactions(const string &ID, const double &beta) : ID(ID), beta(beta) { }
	Interactions() : ID(""), beta(1.0){ }
	bool operator<(const Interactions& other) const { return ID.length() < other.ID.length() || ID < other.ID; }
};

	wxGridSimlaInteractions(int locCount) : wxGridBetaWise(locCount) { ResizeContents(locCount);}
	virtual int GetNumberRows();
	virtual int GetNumberCols();
	virtual bool IsEmptyCell(int row, int col);
	virtual wxString GetValue(int row, int col);

	virtual void SetValue(int row, int col, const wxString& value);
	virtual wxString GetRowLabelValue( int col );
	virtual wxString GetColLabelValue( int col );

    virtual wxString GetTypeName( int row, int col );
    virtual bool CanGetValueAs( int row, int col, const wxString& typeName );
    virtual bool CanSetValueAs( int row, int col, const wxString& typeName );

    virtual double GetValueAsDouble( int row, int col );
	virtual void SetValueAsDouble( int row, int col, double value );

	void ResizeContents(int locusCount);
	void Append(string prev, char letter, int curIdx);

	void Refresh();
	void Reset();
	
	void Load(FileToMap& file, int lociCount);	
	void Save(ostream& file);

	void WriteToStream(ostream &of);
	void LoadFromStream(istream &file);

	bool Verify();
protected:
	vector<string> interactionNames;
	map<string, Interactions> interactions;
};

/*!
 * wxDlgSIMLA class declaration
 */

class wxDlgSIMLA: public wxDiseaseConfiguration {    
    DECLARE_DYNAMIC_CLASS( wxDlgSIMLA )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgSIMLA();
    wxDlgSIMLA( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSIMLA_IDNAME, const wxPoint& pos = SYMBOL_WXDLGSIMLA_POSITION, const wxSize& size = SYMBOL_WXDLGSIMLA_SIZE, long style = SYMBOL_WXDLGSIMLA_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSIMLA_IDNAME, const wxPoint& pos = SYMBOL_WXDLGSIMLA_POSITION, const wxSize& size = SYMBOL_WXDLGSIMLA_SIZE, long style = SYMBOL_WXDLGSIMLA_STYLE );

    /// Destructor
    ~wxDlgSIMLA();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgSIMLA event handler declarations

    /// wxEVT_SIZE event handler for ID_WXDLGSIMLA
    void OnSize( wxSizeEvent& event );

////@end wxDlgSIMLA event handler declarations

////@begin wxDlgSIMLA member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgSIMLA member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgSIMLA member variables
    wxSpinCtrl* dialLocusCount;
    wxTextCtrl* txtTargetPrevalence;
    wxGrid* grdLocusDefinition;
    wxGrid* grdInteractions;
////@end wxDlgSIMLA member variables
	void OnDialLocusSpin(wxSpinEvent& event);
	wxGridSimlaLoci *locGrdMaster;
	wxGridSimlaInteractions *intGrdMaster;
	void InitGrids();
	string GetImportFileTypes() { return "SIMLA Cfg(*.simla)|*.simla"; }

	int GetLocusCount() { return locGrdMaster->GetNumberRows(); }
	std::string GetPreferredExtension(){ return "simla"; }
	bool CanImport() { return true; }
	bool Save(const char *filename, const char *description);
	virtual bool Evaluate();
	virtual bool CanEvaluate() { return false; }

	bool Import(const char *filename, std::string& desc);
};

}

}


#endif
    // _WXDLGSIMLA_H_
