/////////////////////////////////////////////////////////////////////////////
// Name:        wxpanelchromosomerep.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 13:27:12 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPANELCHROMOSOMEREP_H_
#define _WXPANELCHROMOSOMEREP_H_


/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */



////@begin forward declarations
class wxPanelChromosomeRep;
////@end forward declarations
#include <vector>
#include "simulation/locus.h"
namespace GenomeSIM {

namespace GUI {

using namespace Simulation;

using namespace std;

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_CHROMOSOME_VISUAL 10003
#define SYMBOL_WXPANELCHROMOSOMEREP_STYLE wxRAISED_BORDER|wxTAB_TRAVERSAL
#define SYMBOL_WXPANELCHROMOSOMEREP_IDNAME ID_CHROMOSOME_VISUAL
#define SYMBOL_WXPANELCHROMOSOMEREP_SIZE wxSize(-1, 32)
#define SYMBOL_WXPANELCHROMOSOMEREP_POSITION wxDefaultPosition
////@end control identifiers


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
    wxPanelChromosomeRep(wxWindow* parent, wxWindowID id = ID_CHROMOSOME_VISUAL, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxSUNKEN_BORDER|wxTAB_TRAVERSAL);

    /// Creation
    bool Create(wxWindow* parent, wxWindowID id = ID_CHROMOSOME_VISUAL, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxSUNKEN_BORDER|wxTAB_TRAVERSAL);

    /// Destructor
    ~wxPanelChromosomeRep();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPanelChromosomeRep event handler declarations

    /// wxEVT_SIZE event handler for ID_CHROMOSOME_VISUAL
    void OnSize( wxSizeEvent& event );

    /// wxEVT_PAINT event handler for ID_CHROMOSOME_VISUAL
    void OnPaint( wxPaintEvent& event );

    /// wxEVT_MOTION event handler for ID_CHROMOSOME_VISUAL
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

	int GetLocusCount() { return loci.size(); }
	void SetLoci(vector<Locus>& loci);
	void Refresh_Image();
protected:
	vector<Locus> loci;
};

}

}

#endif
    // _WXPANELCHROMOSOMEREP_H_
