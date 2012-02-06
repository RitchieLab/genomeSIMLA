/////////////////////////////////////////////////////////////////////////////
// Name:        wxdialogtaskalert.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed  9 Jul 16:35:46 2008
// RCS-ID:      
// Copyright:   Maylyn Ritchie
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDIALOGTASKALERT_H_
#define _WXDIALOGTASKALERT_H_
#include <pthread.h>

#include <iostream>
namespace GenomeSIM {

namespace GUI {

/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXDIALOGTASKALERT 10002
#define ID_TEXTCTRL 10003
#define SYMBOL_WXDIALOGTASKALERT_STYLE wxCAPTION|wxRESIZE_BORDER
#define SYMBOL_WXDIALOGTASKALERT_TITLE _("DialogTaskAlert")
#define SYMBOL_WXDIALOGTASKALERT_IDNAME ID_WXDIALOGTASKALERT
#define SYMBOL_WXDIALOGTASKALERT_SIZE wxSize(600, 400)
#define SYMBOL_WXDIALOGTASKALERT_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * wxDialogTaskAlert class declaration
 */



class wxDialogTaskAlert: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( wxDialogTaskAlert )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxDialogTaskAlert();
    wxDialogTaskAlert( wxWindow* parent, bool showLog=false, wxWindowID id = SYMBOL_WXDIALOGTASKALERT_IDNAME, const wxString& caption = SYMBOL_WXDIALOGTASKALERT_TITLE, const wxPoint& pos = SYMBOL_WXDIALOGTASKALERT_POSITION, const wxSize& size = SYMBOL_WXDIALOGTASKALERT_SIZE, long style = SYMBOL_WXDIALOGTASKALERT_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXDIALOGTASKALERT_IDNAME, const wxString& caption = SYMBOL_WXDIALOGTASKALERT_TITLE, const wxPoint& pos = SYMBOL_WXDIALOGTASKALERT_POSITION, const wxSize& size = SYMBOL_WXDIALOGTASKALERT_SIZE, long style = SYMBOL_WXDIALOGTASKALERT_STYLE );

    /// Destructor
    ~wxDialogTaskAlert();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxDialogTaskAlert event handler declarations

    /// wxEVT_PAINT event handler for ID_WXDIALOGTASKALERT
    void OnPaint( wxPaintEvent& event );

////@end wxDialogTaskAlert event handler declarations

////@begin wxDialogTaskAlert member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxDialogTaskAlert member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

	void SetMessage(const char* message);
	void SetWindowTitle(const char *title);
	void WriteLog(const char *message);

	wxStaticText *messageText;
	wxTextCtrl *log;

	bool wasCancelled;

	void ShowCompleted();

	wxButton *cmdOK;
	wxButton *cmdSaveLog;
	wxButton *cmdCancel;


	void OnCancelClick(wxCommandEvent& event);
	void OnSaveLog(wxCommandEvent& event);
	bool hasFinished;

protected:
	pthread_t dialogThread;
	bool continueRunning;
	bool showLog;
};


}

}

#endif
