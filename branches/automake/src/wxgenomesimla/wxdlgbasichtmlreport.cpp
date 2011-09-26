/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgbasichtmlreport.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 14 Mar 2008 15:21:17 CDT
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

#include "wxdlgbasichtmlreport.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxDlgBasicHtmlReport type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgBasicHtmlReport, wxDialog )


/*!
 * wxDlgBasicHtmlReport event table definition
 */

BEGIN_EVENT_TABLE( wxDlgBasicHtmlReport, wxDialog )

////@begin wxDlgBasicHtmlReport event table entries
////@end wxDlgBasicHtmlReport event table entries

END_EVENT_TABLE()


/*!
 * wxDlgBasicHtmlReport constructors
 */

wxDlgBasicHtmlReport::wxDlgBasicHtmlReport()
{
    Init();
}

wxDlgBasicHtmlReport::wxDlgBasicHtmlReport( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgBasicHtmlReport creator
 */

bool wxDlgBasicHtmlReport::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgBasicHtmlReport creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgBasicHtmlReport creation
    return true;
}


/*!
 * wxDlgBasicHtmlReport destructor
 */

wxDlgBasicHtmlReport::~wxDlgBasicHtmlReport()
{
////@begin wxDlgBasicHtmlReport destruction
////@end wxDlgBasicHtmlReport destruction
}


/*!
 * Member initialisation
 */

void wxDlgBasicHtmlReport::Init()
{
////@begin wxDlgBasicHtmlReport member initialisation
    htmlReportWindow = NULL;
////@end wxDlgBasicHtmlReport member initialisation
}


/*!
 * Control creation for wxDlgBasicHtmlReport
 */

void wxDlgBasicHtmlReport::CreateControls()
{    
////@begin wxDlgBasicHtmlReport content construction
    wxDlgBasicHtmlReport* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    htmlReportWindow = new wxHtmlWindow( itemDialog1, ID_HTMLWINDOW, wxDefaultPosition, wxSize(500, 300), wxHW_SCROLLBAR_AUTO|wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    itemBoxSizer2->Add(htmlReportWindow, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxALIGN_RIGHT|wxALL, 5);

    wxButton* itemButton5 = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxDlgBasicHtmlReport content construction
}


/*!
 * Should we show tooltips?
 */

bool wxDlgBasicHtmlReport::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgBasicHtmlReport::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgBasicHtmlReport bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgBasicHtmlReport bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgBasicHtmlReport::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgBasicHtmlReport icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgBasicHtmlReport icon retrieval
}



void wxDlgBasicHtmlReport::Append(const char *reportData) {
	report<<reportData<<"\n";
	Refresh();
}


void wxDlgBasicHtmlReport::Clear() {
	report.str("");
}


void wxDlgBasicHtmlReport::Refresh() {
	if 	(htmlReportWindow) 
		htmlReportWindow->SetPage(wxT(report.str().c_str()));

}

}

}
