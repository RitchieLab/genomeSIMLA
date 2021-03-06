/////////////////////////////////////////////////////////////////////////////
// Name:        wxpageinfo.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 06 Dec 2007 09:38:49 AM CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// Generated by DialogBlocks (Personal Edition), Thu 06 Dec 2007 09:38:49 AM CST

#if defined(__GNUG__) && !defined(__APPLE__)
#pragma implementation "wxpageinfo.h"
#endif

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

#include "wxpageinfo.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxPageInfo type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPageInfo, wxPanel )

/*!
 * wxPageInfo event table definition
 */

BEGIN_EVENT_TABLE( wxPageInfo, wxPanel )

////@begin wxPageInfo event table entries
////@end wxPageInfo event table entries

END_EVENT_TABLE()

/*!
 * wxPageInfo constructors
 */

wxPageInfo::wxPageInfo( )
{
}

wxPageInfo::wxPageInfo( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Create(parent, id, pos, size, style);
}

/*!
 * wxPageInfo creator
 */

bool wxPageInfo::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPageInfo member initialisation
    txtMessage = NULL;
////@end wxPageInfo member initialisation

////@begin wxPageInfo creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPageInfo creation
    return TRUE;
}

/*!
 * Control creation for wxPageInfo
 */

void wxPageInfo::CreateControls()
{    
////@begin wxPageInfo content construction
    // Generated by DialogBlocks, Mon 10 Dec 2007 14:13:31 CST (Personal Edition)

    wxPageInfo* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer4, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMessage = new wxStaticText( itemPanel1, wxID_STATIC, _("Set Message Using SetMessage(str)"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(txtMessage, 0, wxALIGN_CENTER_HORIZONTAL|wxALL|wxADJUST_MINSIZE, 5);

////@end wxPageInfo content construction
}

/*!
 * Should we show tooltips?
 */

bool wxPageInfo::ShowToolTips()
{
    return TRUE;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPageInfo::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPageInfo bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPageInfo bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPageInfo::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPageInfo icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPageInfo icon retrieval
}

}

}

