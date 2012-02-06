/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagebbchrom.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 11:41:48 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "wxwizpagebbchrom.h"
#include "simulation/chrompool.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageBBChrom type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageBBChrom, wxPageBBChrom )


/*!
 * wxWizPageBBChrom event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageBBChrom, wxPageBBChrom )

////@begin wxWizPageBBChrom event table entries
////@end wxWizPageBBChrom event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageBBChrom constructors
 */

wxWizPageBBChrom::wxWizPageBBChrom()
{
    Init();
}

wxWizPageBBChrom::wxWizPageBBChrom(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style)
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxWizPageBBChrom creator
 */

bool wxWizPageBBChrom::Create(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style)
{
////@begin wxWizPageBBChrom creation
    wxPageBBChrom::Create(parent, id, pos, size, style);
    CreateControls();
////@end wxWizPageBBChrom creation
    return true;
}


/*!
 * wxWizPageBBChrom destructor
 */

wxWizPageBBChrom::~wxWizPageBBChrom()
{
////@begin wxWizPageBBChrom destruction
////@end wxWizPageBBChrom destruction
}


/*!
 * Member initialisation
 */

void wxWizPageBBChrom::Init()
{
////@begin wxWizPageBBChrom member initialisation
////@end wxWizPageBBChrom member initialisation
}


/*!
 * Control creation for wxWizPageBBChrom
 */

void wxWizPageBBChrom::CreateControls()
{    
////@begin wxWizPageBBChrom content construction
////@end wxWizPageBBChrom content construction
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageBBChrom::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageBBChrom::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageBBChrom bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageBBChrom bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageBBChrom::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageBBChrom icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageBBChrom icon retrieval
}

}

}
