/////////////////////////////////////////////////////////////////////////////
// Name:        wxpageloadchromosome.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 06 Dec 2007 09:20:09 AM CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// Generated by DialogBlocks (Personal Edition), Thu 06 Dec 2007 09:20:09 AM CST

#ifndef _WXPAGELOADCHROMOSOME_H_
#define _WXPAGELOADCHROMOSOME_H_

#if defined(__GNUG__) && !defined(__APPLE__)
#pragma interface "wxpageloadchromosome.cpp"
#endif

/*!
 * Includes
 */

////@begin includes
#include "wx/grid.h"
////@end includes
#include "wxchromcfgdialog.h"
#include "locusmanager.h"
#include "wxgridloci.h"
#include "wxpanelchromosomerep.h"

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxGrid;
class wxPanelChromosomeRep;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PAGE_LOAD_CHROMOSOME 10040
#define ID_TXT_LOCUS_REPORT_FILENAME 10001
#define ID_CMD_SELECT_LOCUS_REPORT 10002
#define ID_TXT_CHROM_LABEL 10000
#define ID_GRID1 10041
#define SYMBOL_WXPAGELOADCHROMOSOME_STYLE wxRESIZE_BORDER
#define SYMBOL_WXPAGELOADCHROMOSOME_TITLE _("Predefined Chromosome")
#define SYMBOL_WXPAGELOADCHROMOSOME_IDNAME ID_PAGE_LOAD_CHROMOSOME
#define SYMBOL_WXPAGELOADCHROMOSOME_SIZE wxSize(600, 400)
#define SYMBOL_WXPAGELOADCHROMOSOME_POSITION wxDefaultPosition
////@end control identifiers

/*!
 * Compatibility
 */

#ifndef wxCLOSE_BOX
#define wxCLOSE_BOX 0x1000
#endif
#ifndef wxFIXED_MINSIZE
#define wxFIXED_MINSIZE 0
#endif

namespace GenomeSIM {

namespace GUI {

/*!
 * wxPageLoadChromosome class declaration
 */

class wxPageLoadChromosome: public wxPanel, public wxChromCfgDialog {    
    DECLARE_DYNAMIC_CLASS( wxPageLoadChromosome )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageLoadChromosome( );
    wxPageLoadChromosome( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGELOADCHROMOSOME_IDNAME, const wxPoint& pos = SYMBOL_WXPAGELOADCHROMOSOME_POSITION, const wxSize& size = SYMBOL_WXPAGELOADCHROMOSOME_SIZE, long style = SYMBOL_WXPAGELOADCHROMOSOME_STYLE );


	~wxPageLoadChromosome();

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGELOADCHROMOSOME_IDNAME, const wxPoint& pos = SYMBOL_WXPAGELOADCHROMOSOME_POSITION, const wxSize& size = SYMBOL_WXPAGELOADCHROMOSOME_SIZE, long style = SYMBOL_WXPAGELOADCHROMOSOME_STYLE );

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageLoadChromosome event handler declarations

    /// wxEVT_SIZE event handler for ID_PAGE_LOAD_CHROMOSOME
    void OnSize( wxSizeEvent& event );

    /// wxEVT_UPDATE_UI event handler for wxID_STATIC
    void OnStaticUpdate( wxUpdateUIEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_LOCUS_REPORT_FILENAME
    void OnTxtLocusReportFilenameUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_ENTER event handler for ID_TXT_LOCUS_REPORT_FILENAME
    void OnTxtLocusReportFilenameEnter( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_SELECT_LOCUS_REPORT
    void OnCmdSelectLocusReportClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_CHROM_LABEL
    void OnTxtChromLabelTextUpdated( wxCommandEvent& event );

    /// wxEVT_GRID_EDITOR_HIDDEN event handler for ID_GRID1
    void OnEditorHidden( wxGridEvent& event );

////@end wxPageLoadChromosome event handler declarations

////@begin wxPageLoadChromosome member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageLoadChromosome member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageLoadChromosome member variables
    wxTextCtrl* txtLocusReportFilename;
    wxTextCtrl* txtChromLabel;
    wxGrid* gridLoci;
    wxPanelChromosomeRep* panelChromosome;
////@end wxPageLoadChromosome member variables
	void InitLabel(const char *label, uint chromID);
	
	bool SetChromosome(FileBasedChromosome *chrom);
	void RefreshSize();
	bool SelectLocusReport();

	void Commit();
	void RefreshSettings() {}
	bool HasChanged();
protected:
	bool InitGrid();

	void SetChromFilename(const char *filename);

	wxGridLoci 			*gridMaster;	///<The structure driving the grid contents
	FileBasedChromosome *chrom;			///<Manages the loci

};

}

}

/*!
 * wxPanelChromosomeRep class declaration
 */

class wxPanelChromosomeRep: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( wxPanelChromosomeRep )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPanelChromosomeRep();
    wxPanelChromosomeRep(wxWindow* parent, wxWindowID id = ID_CHROMOSOME_VISUAL, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize(-1, 32), long style = wxRAISED_BORDER|wxTAB_TRAVERSAL);

    /// Creation
    bool Create(wxWindow* parent, wxWindowID id = ID_CHROMOSOME_VISUAL, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize(-1, 32), long style = wxRAISED_BORDER|wxTAB_TRAVERSAL);

    /// Destructor
    ~wxPanelChromosomeRep();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPanelChromosomeRep event handler declarations

    /// wxEVT_SIZE event handler for ID_CHROMOSOMEREP
    void OnSize( wxSizeEvent& event );

    /// wxEVT_PAINT event handler for ID_CHROMOSOMEREP
    void OnPaint( wxPaintEvent& event );

    /// wxEVT_MOTION event handler for ID_CHROMOSOMEREP
    void OnMotion( wxMouseEvent& event );

////@end wxPanelChromosomeRep event handler declarations

////@begin wxPanelChromosomeRep member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPanelChromosomeRep member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPanelChromosomeRep member variables
////@end wxPanelChromosomeRep member variables
};

#endif
    // _WXPAGELOADCHROMOSOME_H_
