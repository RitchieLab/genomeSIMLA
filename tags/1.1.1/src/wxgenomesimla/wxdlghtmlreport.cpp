/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlghtmlreport.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 14 Mar 2008 14:52:06 CDT
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

#include "wxdlghtmlreport.h"

////@begin XPM images

////@end XPM images
namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgHtmlReport type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgHtmlReport, wxDialog )


/*!
 * wxDlgHtmlReport event table definition
 */

BEGIN_EVENT_TABLE( wxDlgHtmlReport, wxDialog )

////@begin wxDlgHtmlReport event table entries
    EVT_BUTTON( wxID_OK, wxDlgHtmlReport::OnOkClick )

////@end wxDlgHtmlReport event table entries

END_EVENT_TABLE()


/*!
 * wxDlgHtmlReport constructors
 */

wxDlgHtmlReport::wxDlgHtmlReport()
{
    Init();
}

wxDlgHtmlReport::wxDlgHtmlReport( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgHtmlReport creator
 */

bool wxDlgHtmlReport::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgHtmlReport creation
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgHtmlReport creation
    return true;
}


/*!
 * wxDlgHtmlReport destructor
 */

wxDlgHtmlReport::~wxDlgHtmlReport()
{
////@begin wxDlgHtmlReport destruction
////@end wxDlgHtmlReport destruction
}


/*!
 * Member initialisation
 */

void wxDlgHtmlReport::Init()
{
////@begin wxDlgHtmlReport member initialisation
    htmlReportWindow = NULL;
////@end wxDlgHtmlReport member initialisation
}


/*!
 * Control creation for wxDlgHtmlReport
 */

void wxDlgHtmlReport::CreateControls()
{    
////@begin wxDlgHtmlReport content construction
    wxDlgHtmlReport* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    htmlReportWindow = new wxHtmlWindow( itemDialog1, ID_HTMLWINDOW1, wxDefaultPosition, wxSize(500, 300), wxHW_SCROLLBAR_AUTO|wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    itemBoxSizer2->Add(htmlReportWindow, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton5 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgHtmlReport content construction
	Refresh();
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
 */

void wxDlgHtmlReport::OnOkClick( wxCommandEvent& event )
{
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}


/*!
 * Should we show tooltips?
 */

bool wxDlgHtmlReport::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgHtmlReport::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgHtmlReport bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgHtmlReport bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgHtmlReport::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgHtmlReport icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgHtmlReport icon retrieval
}



void wxDlgHtmlReport::Append(const char *reportData) {
	report<<reportData<<"\n";
	Refresh();
}


void wxDlgHtmlReport::Clear() {
	report.str("");
}


void wxDlgHtmlReport::Refresh() {
	if 	(htmlReportWindow) 
		htmlReportWindow->SetPage(wxT(report.str().c_str()));

}

}

}
