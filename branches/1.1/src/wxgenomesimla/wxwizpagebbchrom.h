/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagebbchrom.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 11:41:48 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZPAGEBBCHROM_H_
#define _WXWIZPAGEBBCHROM_H_
#include "wxpagebbchrom.h"

/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxWizPageBBChrom;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PANEL_CHROM_CFG 10132
#define SYMBOL_WXWIZPAGEBBCHROM_STYLE wxSUNKEN_BORDER|wxTAB_TRAVERSAL
#define SYMBOL_WXWIZPAGEBBCHROM_IDNAME ID_PANEL_CHROM_CFG
#define SYMBOL_WXWIZPAGEBBCHROM_SIZE wxDefaultSize
#define SYMBOL_WXWIZPAGEBBCHROM_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxWizPageBBChrom class declaration
 */

class wxWizPageBBChrom: public wxPageBBChrom
{    
    DECLARE_DYNAMIC_CLASS( wxWizPageBBChrom )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizPageBBChrom();
    wxWizPageBBChrom(wxWindow* parent, wxWindowID id = ID_PANEL_CHROM_CFG, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxSUNKEN_BORDER|wxTAB_TRAVERSAL);

    /// Creation
    bool Create(wxWindow* parent, wxWindowID id = ID_PANEL_CHROM_CFG, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxSUNKEN_BORDER|wxTAB_TRAVERSAL);

    /// Destructor
    ~wxWizPageBBChrom();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizPageBBChrom event handler declarations

////@end wxWizPageBBChrom event handler declarations

////@begin wxWizPageBBChrom member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizPageBBChrom member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizPageBBChrom member variables
////@end wxWizPageBBChrom member variables
};

}

}
#endif
    // _WXWIZPAGEBBCHROM_H_
