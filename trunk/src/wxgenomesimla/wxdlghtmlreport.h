/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlghtmlreport.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 14 Mar 2008 14:52:06 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGHTMLREPORT_H_
#define _WXDLGHTMLREPORT_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/html/htmlwin.h"
////@end includes
#include <sstream>
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxHtmlWindow;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDLGHTMLREPORT 10097
#define ID_HTMLWINDOW1 10106
#define SYMBOL_WXDLGHTMLREPORT_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGHTMLREPORT_TITLE _("HTML Report")
#define SYMBOL_WXDLGHTMLREPORT_IDNAME ID_WXDLGHTMLREPORT
#define SYMBOL_WXDLGHTMLREPORT_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGHTMLREPORT_POSITION wxDefaultPosition
////@end control identifiers
namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxDlgHtmlReport class declaration
 */

class wxDlgHtmlReport: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgHtmlReport )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgHtmlReport();
    wxDlgHtmlReport( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGHTMLREPORT_IDNAME, const wxString& caption = SYMBOL_WXDLGHTMLREPORT_TITLE, const wxPoint& pos = SYMBOL_WXDLGHTMLREPORT_POSITION, const wxSize& size = SYMBOL_WXDLGHTMLREPORT_SIZE, long style = SYMBOL_WXDLGHTMLREPORT_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGHTMLREPORT_IDNAME, const wxString& caption = SYMBOL_WXDLGHTMLREPORT_TITLE, const wxPoint& pos = SYMBOL_WXDLGHTMLREPORT_POSITION, const wxSize& size = SYMBOL_WXDLGHTMLREPORT_SIZE, long style = SYMBOL_WXDLGHTMLREPORT_STYLE );

    /// Destructor
    ~wxDlgHtmlReport();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgHtmlReport event handler declarations
    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
    void OnOkClick( wxCommandEvent& event );

////@end wxDlgHtmlReport event handler declarations

////@begin wxDlgHtmlReport member function declarations
    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgHtmlReport member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgHtmlReport member variables
    wxHtmlWindow* htmlReportWindow;
////@end wxDlgHtmlReport member variables

	void Append(const char *reportData);
	void Clear();
	void Refresh();
protected:
	stringstream report;
};
}
}
#endif
    // _WXDLGHTMLREPORT_H_
