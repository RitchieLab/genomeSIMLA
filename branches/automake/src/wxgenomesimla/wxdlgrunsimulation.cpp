/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgrunsimulation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 14 Jan 2008 12:46:57 CST
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

#include "wxdlgrunsimulation.h"

////@begin XPM images
////@end XPM images


/*!
 * wxDlgRunSimulation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgRunSimulation, wxDialog )


/*!
 * wxDlgRunSimulation event table definition
 */

BEGIN_EVENT_TABLE( wxDlgRunSimulation, wxDialog )

////@begin wxDlgRunSimulation event table entries
////@end wxDlgRunSimulation event table entries

END_EVENT_TABLE()


/*!
 * wxDlgRunSimulation constructors
 */

wxDlgRunSimulation::wxDlgRunSimulation()
{
    Init();
}

wxDlgRunSimulation::wxDlgRunSimulation( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgRunSimulation creator
 */

bool wxDlgRunSimulation::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgRunSimulation creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgRunSimulation creation
    return true;
}


/*!
 * wxDlgRunSimulation destructor
 */

wxDlgRunSimulation::~wxDlgRunSimulation()
{
////@begin wxDlgRunSimulation destruction
////@end wxDlgRunSimulation destruction
}


/*!
 * Member initialisation
 */

void wxDlgRunSimulation::Init()
{
////@begin wxDlgRunSimulation member initialisation
////@end wxDlgRunSimulation member initialisation
}


/*!
 * Control creation for wxDlgRunSimulation
 */

void wxDlgRunSimulation::CreateControls()
{    
////@begin wxDlgRunSimulation content construction
    wxDlgRunSimulation* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxStaticBox* itemStaticBoxSizer3Static = new wxStaticBox(itemDialog1, wxID_ANY, _("Drop Points"));
    wxStaticBoxSizer* itemStaticBoxSizer3 = new wxStaticBoxSizer(itemStaticBoxSizer3Static, wxVERTICAL);
    itemBoxSizer2->Add(itemStaticBoxSizer3, 1, wxALIGN_CENTER_HORIZONTAL|wxALL, 10);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer3->Add(itemBoxSizer4, 1, wxALIGN_RIGHT|wxALL, 5);

    wxStaticText* itemStaticText5 = new wxStaticText( itemDialog1, wxID_STATIC, _("First Drop Point:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemStaticText5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxTextCtrl* itemTextCtrl6 = new wxTextCtrl( itemDialog1, ID_TEXTCTRL2, _T(""), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemTextCtrl6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer3->Add(itemBoxSizer7, 1, wxALIGN_RIGHT|wxALL, 5);

    wxStaticText* itemStaticText8 = new wxStaticText( itemDialog1, wxID_STATIC, _("Drop Frequency:"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer7->Add(itemStaticText8, 1, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxTextCtrl* itemTextCtrl9 = new wxTextCtrl( itemDialog1, ID_TEXTCTRL, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer7->Add(itemTextCtrl9, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer10 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer3->Add(itemBoxSizer10, 1, wxALIGN_RIGHT|wxALL, 5);

    wxStaticText* itemStaticText11 = new wxStaticText( itemDialog1, wxID_STATIC, _("Drop Count:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer10->Add(itemStaticText11, 1, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxTextCtrl* itemTextCtrl12 = new wxTextCtrl( itemDialog1, ID_TEXTCTRL1, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer10->Add(itemTextCtrl12, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

////@end wxDlgRunSimulation content construction
}


/*!
 * Should we show tooltips?
 */

bool wxDlgRunSimulation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgRunSimulation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgRunSimulation bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgRunSimulation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgRunSimulation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgRunSimulation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgRunSimulation icon retrieval
}
