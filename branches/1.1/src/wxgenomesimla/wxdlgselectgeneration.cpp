/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgselectgeneration.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 27 Feb 2008 15:41:39 CST
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

#include "wxdlgselectgeneration.h"

////@begin XPM images

////@end XPM images


/*!
 * wxDlgSelectGeneration type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgSelectGeneration, wxDialog )


/*!
 * wxDlgSelectGeneration event table definition
 */

BEGIN_EVENT_TABLE( wxDlgSelectGeneration, wxDialog )

////@begin wxDlgSelectGeneration event table entries
////@end wxDlgSelectGeneration event table entries

END_EVENT_TABLE()


/*!
 * wxDlgSelectGeneration constructors
 */

wxDlgSelectGeneration::wxDlgSelectGeneration()
{
    Init();
}

wxDlgSelectGeneration::wxDlgSelectGeneration( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgSelectGeneration creator
 */

bool wxDlgSelectGeneration::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgSelectGeneration creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgSelectGeneration creation
    return true;
}


/*!
 * wxDlgSelectGeneration destructor
 */

wxDlgSelectGeneration::~wxDlgSelectGeneration()
{
////@begin wxDlgSelectGeneration destruction
////@end wxDlgSelectGeneration destruction
}


/*!
 * Member initialisation
 */

void wxDlgSelectGeneration::Init()
{
////@begin wxDlgSelectGeneration member initialisation
    treeGenerations = NULL;
////@end wxDlgSelectGeneration member initialisation
}


/*!
 * Control creation for wxDlgSelectGeneration
 */

void wxDlgSelectGeneration::CreateControls()
{    
////@begin wxDlgSelectGeneration content construction
    wxDlgSelectGeneration* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemDialog1, wxID_STATIC, _("Select the generation of interest from the project(s) below. "), wxDefaultPosition, wxDefaultSize, wxST_NO_AUTORESIZE );
    itemBoxSizer2->Add(itemStaticText3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    treeGenerations = new wxPageSelectGeneration( itemDialog1, ID_GENERATION_SELECTION, wxDefaultPosition, wxSize(300, 200), wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
    itemBoxSizer2->Add(treeGenerations, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer5, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton6 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer5->Add(itemButton6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton7 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer5->Add(itemButton7, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgSelectGeneration content construction
}


/*!
 * Should we show tooltips?
 */

bool wxDlgSelectGeneration::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgSelectGeneration::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgSelectGeneration bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgSelectGeneration bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgSelectGeneration::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgSelectGeneration icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgSelectGeneration icon retrieval
}
