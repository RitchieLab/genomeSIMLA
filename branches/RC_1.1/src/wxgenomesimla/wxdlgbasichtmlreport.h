/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgbasichtmlreport.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 14 Mar 2008 15:21:17 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGBASICHTMLREPORT_H_
#define _WXDLGBASICHTMLREPORT_H_


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
#define ID_WXDLGBASICHTMLREPORT 10109
#define ID_HTMLWINDOW 10000
#define SYMBOL_WXDLGBASICHTMLREPORT_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXDLGBASICHTMLREPORT_TITLE _("Html Report")
#define SYMBOL_WXDLGBASICHTMLREPORT_IDNAME ID_WXDLGBASICHTMLREPORT
#define SYMBOL_WXDLGBASICHTMLREPORT_SIZE wxSize(400, 300)
#define SYMBOL_WXDLGBASICHTMLREPORT_POSITION wxDefaultPosition
////@end control identifiers

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxDlgBasicHtmlReport class declaration
 */

class wxDlgBasicHtmlReport: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDlgBasicHtmlReport )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDlgBasicHtmlReport();
    wxDlgBasicHtmlReport( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGBASICHTMLREPORT_IDNAME, const wxString& caption = SYMBOL_WXDLGBASICHTMLREPORT_TITLE, const wxPoint& pos = SYMBOL_WXDLGBASICHTMLREPORT_POSITION, const wxSize& size = SYMBOL_WXDLGBASICHTMLREPORT_SIZE, long style = SYMBOL_WXDLGBASICHTMLREPORT_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDLGBASICHTMLREPORT_IDNAME, const wxString& caption = SYMBOL_WXDLGBASICHTMLREPORT_TITLE, const wxPoint& pos = SYMBOL_WXDLGBASICHTMLREPORT_POSITION, const wxSize& size = SYMBOL_WXDLGBASICHTMLREPORT_SIZE, long style = SYMBOL_WXDLGBASICHTMLREPORT_STYLE );

    /// Destructor
    ~wxDlgBasicHtmlReport();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDlgBasicHtmlReport event handler declarations

////@end wxDlgBasicHtmlReport event handler declarations

////@begin wxDlgBasicHtmlReport member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDlgBasicHtmlReport member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxDlgBasicHtmlReport member variables
    wxHtmlWindow* htmlReportWindow;
////@end wxDlgBasicHtmlReport member variables

	void Append(const char *reportData);
	void Clear();
	void Refresh();
protected:
	stringstream report;
};

}

}
#endif
    // _WXDLGBASICHTMLREPORT_H_
