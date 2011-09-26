/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgselectlocusrange.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 12 Dec 2007 15:13:17 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGSELECTLOCUSRANGE_H_
#define _WXDLGSELECTLOCUSRANGE_H_

#include "wxlocuslistcontrol.h"

typedef GenomeSIM::GUI::wxLocusListControl wxLocusListCtrl;

/*!
 * Includes
 */

////@begin includes
////@end includes


#include "locusmanager.h"
#include "wxgridloci.h"

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxLocusListControl;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGSELECTLOCUSRANGE 10035
#define ID_LOCUSLISTCONTROL 10044
#define SYMBOL_WXDLGSELECTLOCUSRANGE_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGSELECTLOCUSRANGE_TITLE _("Select Locus Range")
#define SYMBOL_WXDLGSELECTLOCUSRANGE_IDNAME ID_WXDLGSELECTLOCUSRANGE
#define SYMBOL_WXDLGSELECTLOCUSRANGE_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGSELECTLOCUSRANGE_POSITION wxDefaultPosition
////@end control identifiers

#define ID_MENU_BEGIN 11100
#define ID_MENU_END   11101
#define ID_LST_LOCUS_LIST 11102

namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgSelectLocusRange class declaration
 */

class wxDlgSelectLocusRange: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgSelectLocusRange )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgSelectLocusRange();
    wxDlgSelectLocusRange( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSELECTLOCUSRANGE_IDNAME, const wxString& caption = SYMBOL_WXDLGSELECTLOCUSRANGE_TITLE, const wxPoint& pos = SYMBOL_WXDLGSELECTLOCUSRANGE_POSITION, const wxSize& size = SYMBOL_WXDLGSELECTLOCUSRANGE_SIZE, long style = SYMBOL_WXDLGSELECTLOCUSRANGE_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGSELECTLOCUSRANGE_IDNAME, const wxString& caption = SYMBOL_WXDLGSELECTLOCUSRANGE_TITLE, const wxPoint& pos = SYMBOL_WXDLGSELECTLOCUSRANGE_POSITION, const wxSize& size = SYMBOL_WXDLGSELECTLOCUSRANGE_SIZE, long style = SYMBOL_WXDLGSELECTLOCUSRANGE_STYLE );

    /// Destructor
    ~wxDlgSelectLocusRange();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgSelectLocusRange event handler declarations
	/// wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LOCUSLISTCONTROL
	void OnLocuslistcontrolSelected( wxListEvent& event );

	/// wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LOCUSLISTCONTROL
	void OnLocuslistcontrolItemActivated( wxListEvent& event );

    /// wxEVT_SIZE event handler for ID_LOCUSLISTCONTROL
    void OnSize( wxSizeEvent& event );

////@end wxDlgSelectLocusRange event handler declarations

////@begin wxDlgSelectLocusRange member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgSelectLocusRange member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgSelectLocusRange member variables
    wxStaticText* txtChromosomeLabel;
    wxLocusListControl* lstLoci;
////@end wxDlgSelectLocusRange member variables

	void Initialize(FileBasedChromosome *chrom, Locus *loc);
	void RefreshSize();

	void SetEnd(wxCommandEvent &event);
	void SetBegin(wxCommandEvent &event);

	void EndOK();
	
	Locus *GetSelection();
protected:

	FileBasedChromosome *chrom;			///<Manages the loci
	Locus 				*curSelection;	///<The locus of the current selection
};

}

}

#endif
    // _WXDLGSELECTLOCUSRANGE_H_
