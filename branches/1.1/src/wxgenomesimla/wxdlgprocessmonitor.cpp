/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgprocessmonitor.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 18 Feb 2008 17:43:42 CST
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

#include "wxdlgprocessmonitor.h"
#include <sstream>
////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxDlgProcessMonitor type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgProcessMonitor, wxDialog )


/*!
 * wxDlgProcessMonitor event table definition
 */

BEGIN_EVENT_TABLE( wxDlgProcessMonitor, wxDialog )

////@begin wxDlgProcessMonitor event table entries
    EVT_BUTTON( ID_CMD_RUN, wxDlgProcessMonitor::OnCmdRunClick )

    EVT_BUTTON( wxID_CANCEL, wxDlgProcessMonitor::OnCancelClick )

////@end wxDlgProcessMonitor event table entries

END_EVENT_TABLE()


/*!
 * wxDlgProcessMonitor constructors
 */

wxDlgProcessMonitor::wxDlgProcessMonitor()
{
    Init();
}

wxDlgProcessMonitor::wxDlgProcessMonitor( wxWindow* parent, AppController *ctrl, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
	appController = ctrl;
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgProcessMonitor creator
 */

bool wxDlgProcessMonitor::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgProcessMonitor creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgProcessMonitor creation
    return true;
}


/*!
 * wxDlgProcessMonitor destructor
 */

wxDlgProcessMonitor::~wxDlgProcessMonitor()
{
////@begin wxDlgProcessMonitor destruction
////@end wxDlgProcessMonitor destruction
}


/*!
 * Member initialisation
 */

void wxDlgProcessMonitor::Init()
{
	isCompleted = false;
////@begin wxDlgProcessMonitor member initialisation
    txtSimLog = NULL;
    guageSimCompletion = NULL;
    cmdRun = NULL;
    cmdCancel = NULL;
    statusBar = NULL;
////@end wxDlgProcessMonitor member initialisation
}


/*!
 * Control creation for wxDlgProcessMonitor
 */

void wxDlgProcessMonitor::CreateControls()
{
////@begin wxDlgProcessMonitor content construction
    wxDlgProcessMonitor* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    txtSimLog = new wxTextCtrl( itemDialog1, ID_TXT_SIMLOG, _T(""), wxDefaultPosition, wxSize(450, 250), wxTE_MULTILINE );
    txtSimLog->SetFont(wxFont(10, wxTELETYPE, wxNORMAL, wxNORMAL, false, wxT("Courier 10 Pitch")));
    itemBoxSizer2->Add(txtSimLog, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxGROW|wxALL, 5);

    guageSimCompletion = new wxGauge( itemDialog1, ID_GUAGE_PROCESS_COMPLETION, 100, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
    guageSimCompletion->SetValue(1);
    itemBoxSizer4->Add(guageSimCompletion, 1, wxGROW|wxALL, 5);

    cmdRun = new wxButton( itemDialog1, ID_CMD_RUN, _("Run"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(cmdRun, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdCancel = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(cmdCancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    statusBar = new wxStatusBar( itemDialog1, ID_STATUSBAR, wxST_SIZEGRIP|wxNO_BORDER );
    statusBar->SetFieldsCount(1);
    itemBoxSizer2->Add(statusBar, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

////@end wxDlgProcessMonitor content construction
	
	stringstream ss;
	appController->SummarizeSimulation(ss);
	txtSimLog->WriteText(_T(ss.str().c_str()));

}


/*!
 * Should we show tooltips?
 */

bool wxDlgProcessMonitor::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgProcessMonitor::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgProcessMonitor bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgProcessMonitor bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgProcessMonitor::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgProcessMonitor icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgProcessMonitor icon retrieval
}


/*!
* wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_RUN
 */

void wxDlgProcessMonitor::OnCmdRunClick( wxCommandEvent& event )
{
	if (isCompleted) {
		EndModal(1);
		return;
	}
	
	if (appController)  {
		cmdRun->Enable(false);
		statusBar->PushStatusText(_T("genomeSIMLA is running. Please be patient while the simulation completes."));
		try {
			RunAnalysis();
		/**
		 * @todo Provide a file selection dialog to allow users to select the appropriate pool. 
	 	 */
		} catch (Exception::FileNotFound& e) {
			stringstream ss;
			ss<<"The pool you selected to draw data from, "<<e.filename<<", was not found. Please make sure that the filenames haven't been changed.";
			wxMessageBox(_(e.GetErrorMessage().c_str()), _T("Exception encountered"), wxICON_WARNING | wxOK, this);
			return;
		} catch (...) {
			cout<<"Unknown error when opening a previous pool. We are going away\n";
			assert(0);
		}


	}
	cmdRun->Enable(true);
	cmdRun->SetLabel(_T("OK"));
}


bool wxDlgProcessMonitor::RunAnalysis() {
	stringstream initstream;
	continueRunning = true;

	bool currSampledLDSetting = appController->parameters.FastLD();
	appController->parameters.FastLD(false);
	
	SetCursor(wxCursor(wxCURSOR_WAIT));
	appController->InitializeExecution(initstream);
	txtSimLog->WriteText(_T(initstream.str().c_str()));

	Update();
		
	appController->PerformAnalysis();
	do {
		txtSimLog->WriteText(_T(appController->GetSimOutput().c_str()));
		Update();
		wxYield();
		wxMilliSleep(50);
	} while (appController->SimIsRunning());
		
	SetCursor(wxCursor(wxCURSOR_ARROW));
	
	if (!continueRunning)
		EndModal(0);
	else {
		isCompleted = true;
		txtSimLog->WriteText(_T("Analysis Completed. \n"));
		cmdRun->Enable(true);
		statusBar->PopStatusText();
	}
	isCompleted = true;

	appController->parameters.FastLD(currSampledLDSetting);

}




/*!
* wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON1
 */

void wxDlgProcessMonitor::OnCancelClick( wxCommandEvent& event )
{
	if (appController->SimIsRunning()) {
		statusBar->PopStatusText();
		statusBar->PushStatusText(_T("cancelling run. Please wait while genomeSIMLA can cleanly stop processing"));
		continueRunning = false;
		appController->HaltExecution();
	}
	else
		event.Skip();
}



IMPLEMENT_DYNAMIC_CLASS( wxDlgDatasetGenerationMonitor, wxDialog )
/*!
 * wxDlgProcessMonitor event table definition
 */

BEGIN_EVENT_TABLE( wxDlgDatasetGenerationMonitor, wxDialog )

////@begin wxDlgProcessMonitor event table entries
	EVT_BUTTON( ID_CMD_RUN, wxDlgDatasetGenerationMonitor::OnCmdRunClick )
	
    EVT_BUTTON( wxID_CANCEL, wxDlgDatasetGenerationMonitor::OnCancelClick )	
////@end wxDlgProcessMonitor event table entries

END_EVENT_TABLE()

wxDlgDatasetGenerationMonitor::wxDlgDatasetGenerationMonitor( wxWindow* parent, AppController *ctrl, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : wxDlgProcessMonitor(parent, ctrl, id, caption, pos, size, style) 
{
}


bool wxDlgDatasetGenerationMonitor::RunAnalysis() {
	stringstream initstream;
	continueRunning = true;

	SetCursor(wxCursor(wxCURSOR_WAIT));
//	size_t startingGeneration = appController->InitializeExecution(initstream);

//	cout<<"Starting Generation: "<<startingGeneration<<"\n";
	
	txtSimLog->WriteText(_T(initstream.str().c_str()));
	guageSimCompletion->SetValue(0);
	Update();
		
	appController->parameters.WriteConfiguration();
	appController->WriteDatasets();
	do {
		txtSimLog->WriteText(_T(appController->GetSimOutput().c_str()));
		guageSimCompletion->SetValue((int)(appController->GetPercentCompleted() * 100.0));
		Update();
		wxYield();
		wxMilliSleep(50);
	} while (appController->SimIsRunning());
	guageSimCompletion->SetValue(100);
	SetCursor(wxCursor(wxCURSOR_ARROW));
	
	if (!continueRunning)
		EndModal(0);
	else {
		isCompleted = true;
		txtSimLog->WriteText(_T("Datasets Completed. \n"));
		cmdRun->Enable(true);
	}
	isCompleted = true;


}

}

}
