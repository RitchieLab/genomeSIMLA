/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsavecurrent.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Dec 2007 14:47:44 CST
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

#include "wxdlgsavecurrent.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgSaveCurrent type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgSaveCurrent, wxDialog )


/*!
 * wxDlgSaveCurrent event table definition
 */

BEGIN_EVENT_TABLE( wxDlgSaveCurrent, wxDialog )

////@begin wxDlgSaveCurrent event table entries
    EVT_BUTTON( wxID_CLOSE, wxDlgSaveCurrent::OnCloseClick )

    EVT_BUTTON( wxID_SAVEAS, wxDlgSaveCurrent::OnSaveAsClick )

    EVT_BUTTON( wxID_SAVE, wxDlgSaveCurrent::OnSaveClick )

////@end wxDlgSaveCurrent event table entries

END_EVENT_TABLE()


/*!
 * wxDlgSaveCurrent constructors
 */

wxDlgSaveCurrent::wxDlgSaveCurrent()
{
    Init();
}

wxDlgSaveCurrent::wxDlgSaveCurrent( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgSaveCurrent creator
 */

bool wxDlgSaveCurrent::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgSaveCurrent creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgSaveCurrent creation
    return true;
}


/*!
 * wxDlgSaveCurrent destructor
 */

wxDlgSaveCurrent::~wxDlgSaveCurrent()
{
////@begin wxDlgSaveCurrent destruction
////@end wxDlgSaveCurrent destruction
}


/*!
 * Member initialisation
 */

void wxDlgSaveCurrent::Init()
{
////@begin wxDlgSaveCurrent member initialisation
////@end wxDlgSaveCurrent member initialisation
}


/*!
 * Control creation for wxDlgSaveCurrent
 */

void wxDlgSaveCurrent::CreateControls()
{    
////@begin wxDlgSaveCurrent content construction
    wxDlgSaveCurrent* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemDialog1, wxID_STATIC, _("Changes have been made to the current file. Do you want to save those changes?"), wxDefaultPosition, wxSize(300, -1), 0 );
    itemBoxSizer2->Add(itemStaticText3, 0, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxButton* itemButton5 = new wxButton( itemDialog1, wxID_CLOSE, _("Don't Save"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton6 = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer4->Add(25, 5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton8 = new wxButton( itemDialog1, wxID_SAVEAS, _("Save &As..."), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton8, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton9 = new wxButton( itemDialog1, wxID_SAVE, _("&Save"), wxDefaultPosition, wxDefaultSize, 0 );
    itemButton9->SetName(_T("cmdSave"));
    itemBoxSizer4->Add(itemButton9, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer2->Add(25, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

////@end wxDlgSaveCurrent content construction
}


/*!
 * Should we show tooltips?
 */

bool wxDlgSaveCurrent::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgSaveCurrent::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgSaveCurrent bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgSaveCurrent bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgSaveCurrent::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgSaveCurrent icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgSaveCurrent icon retrieval
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CLOSE
 */

void wxDlgSaveCurrent::OnCloseClick( wxCommandEvent& event )
{
	EndModal(wxID_CLOSE);
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CLOSE in wxDlgSaveCurrent.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_CLOSE in wxDlgSaveCurrent. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVEAS
 */

void wxDlgSaveCurrent::OnSaveAsClick( wxCommandEvent& event )
{
	EndModal(wxID_SAVEAS);
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVEAS in wxDlgSaveCurrent.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVEAS in wxDlgSaveCurrent. 
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVE
 */

void wxDlgSaveCurrent::OnSaveClick( wxCommandEvent& event )
{
	EndModal(wxID_SAVE);
////@begin wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVE in wxDlgSaveCurrent.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_SAVE in wxDlgSaveCurrent. 
}
}

}


