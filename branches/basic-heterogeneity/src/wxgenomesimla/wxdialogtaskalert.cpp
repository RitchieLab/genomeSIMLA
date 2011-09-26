//
// C++ Implementation: wxdialogtaskalert
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <wx/wx.h>
#include "wxdialogtaskalert.h"
#include <iostream>
#include <fstream>
#include <wx/animate.h>
#include <wx/filename.h>
#include <iostream>
namespace GenomeSIM {

namespace GUI {

using namespace std;

//wxDialogTaskAlert *wxDialogTaskAlert::alertDialog = NULL;


/*!
 * wxDialogTaskAlert type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDialogTaskAlert, wxDialog )


/*!
 * wxDialogTaskAlert event table definition
 */

BEGIN_EVENT_TABLE( wxDialogTaskAlert, wxDialog )

	EVT_BUTTON( wxID_CANCEL, wxDialogTaskAlert::OnCancelClick )
	EVT_BUTTON( 1032, wxDialogTaskAlert::OnSaveLog )

END_EVENT_TABLE()


/*!
 * wxDialogTaskAlert constructors
 */

wxDialogTaskAlert::wxDialogTaskAlert() 
{
    Init();
}

wxDialogTaskAlert::wxDialogTaskAlert( wxWindow* parent, bool showLog, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : showLog(showLog)
{
    Init();
    Create(parent, id, caption, pos, size, style);
}

/*!
 * wxDialogTaskAlert destructor
 */

wxDialogTaskAlert::~wxDialogTaskAlert()
{
}


/*!
 * wxDialogTaskAlert creator
 */

bool wxDialogTaskAlert::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDialogTaskAlert creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDialogTaskAlert creation
    return true;
}




/*!
 * Member initialisation
 */

void wxDialogTaskAlert::Init()
{
	messageText = NULL;
	log 		= NULL;
	wasCancelled = false;
	hasFinished  = false;
	cmdOK		= NULL;
	cmdSaveLog	= NULL;
	cmdCancel 	= NULL;
}


/*!
 * Control creation for wxDialogTaskAlert
 */

void wxDialogTaskAlert::CreateControls()
{    

    wxDialogTaskAlert* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 0, wxGROW|wxALL, 5);
	
	messageText = new wxStaticText( itemDialog1, wxID_STATIC, _("Long mesage that will change size of dialog"), wxDefaultPosition, wxDefaultSize, 0 );
	itemBoxSizer3->Add(messageText, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

	if (showLog) {
		log = new wxTextCtrl( itemDialog1, ID_TEXTCTRL, _(""), wxDefaultPosition, wxSize(300, 200), wxTE_MULTILINE|wxTE_READONLY|wxTE_DONTWRAP);
		itemBoxSizer2->Add(log, 1, wxGROW|wxALL, 5);

		wxBoxSizer *itemBoxSizer6 = new wxBoxSizer(wxHORIZONTAL);
		itemBoxSizer2->Add(itemBoxSizer6, 0, wxALIGN_RIGHT|wxALL, 5);

		cmdSaveLog = new wxButton(itemDialog1, 1032, _("Save Log"));
		cmdSaveLog->Enable(false);
		itemBoxSizer6->Add(cmdSaveLog, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

		cmdCancel = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
		itemBoxSizer6->Add(cmdCancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
	
		cmdOK = new wxButton( itemDialog1, wxID_OK, _("OK"), wxDefaultPosition, wxDefaultSize, 0 );
		cmdOK->Enable(false);
		itemBoxSizer6->Add(cmdOK, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
	}
}


/*!
 * Should we show tooltips?
 */

bool wxDialogTaskAlert::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDialogTaskAlert::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
}

void wxDialogTaskAlert::ShowCompleted() {
	if (cmdCancel)
		cmdCancel->Enable(false);
	
	if (cmdOK)
		cmdOK->Enable(true);
	
	if (cmdSaveLog)
		cmdSaveLog->Enable(true);
}

/*!
 * Get icon resources
 */

wxIcon wxDialogTaskAlert::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDialogTaskAlert icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDialogTaskAlert icon retrieval
}


/*!
 * wxEVT_PAINT event handler for ID_WXDIALOGTASKALERT
 */

void wxDialogTaskAlert::OnPaint( wxPaintEvent& event )
{
	std::cout<<"OnPaint()\n";
    // Before editing this code, remove the block markers.
    wxPaintDC dc(this);
}


void wxDialogTaskAlert::WriteLog(const char *message) {
	if (log)
		*log<<message;
}
void wxDialogTaskAlert::SetMessage(const char* message) {
	std::cout<<"Setting message to: "<<message<<"\n";
	messageText->SetLabel(_(message));
	messageText->Refresh();
	Update();
	wxYield();
}
void wxDialogTaskAlert::SetWindowTitle(const char *title) {
	SetTitle(_(title));
	hasFinished = true;
}


void wxDialogTaskAlert::OnCancelClick(wxCommandEvent& event) {
	SetReturnCode(-1);
}

void wxDialogTaskAlert::OnSaveLog(wxCommandEvent &event) {
	cout<<"Saving to some special place!\n";
	string buffer = log->GetValue().c_str();
	wxFileDialog getLogFilename(this, _("Log filename"), _(""), _(""), _("Log File (*.log)|*.log|All Files (*.*)|*.*"), wxSAVE|wxOVERWRITE_PROMPT);
	if (getLogFilename.ShowModal() == wxID_OK) {
		string filename(getLogFilename.GetPath());
		ofstream file(filename.c_str());
		file<<buffer;
	}
}

}

}
