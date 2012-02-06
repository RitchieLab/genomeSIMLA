/////////////////////////////////////////////////////////////////////////////
// Name:        wxbeginsimulation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Mon 04 Feb 2008 15:56:51 CST
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

#include "wxbeginsimulation.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxBeginSimulation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxBeginSimulation, wxDialog )


/*!
 * wxBeginSimulation event table definition
 */

BEGIN_EVENT_TABLE( wxBeginSimulation, wxDialog )

////@begin wxBeginSimulation event table entries
    EVT_BUTTON( ID_CMD_RUN, wxBeginSimulation::OnCmdRunClick )

    EVT_BUTTON( wxID_CANCEL, wxBeginSimulation::OnCancelClick )

////@end wxBeginSimulation event table entries

END_EVENT_TABLE()


/*!
 * wxBeginSimulation constructors
 */

wxBeginSimulation::wxBeginSimulation() : isCompleted(false)
{
    Init();
}

wxBeginSimulation::wxBeginSimulation( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : isCompleted(false), continueRunning(true)
{
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxBeginSimulation creator
 */

bool wxBeginSimulation::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxBeginSimulation creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxBeginSimulation creation
    return true;
}


/*!
 * wxBeginSimulation destructor
 */

wxBeginSimulation::~wxBeginSimulation()
{
////@begin wxBeginSimulation destruction
////@end wxBeginSimulation destruction
}


/*!
 * Member initialisation
 */

void wxBeginSimulation::Init()
{
////@begin wxBeginSimulation member initialisation
    txtSimLog = NULL;
    guageSimCompletion = NULL;
    cmdRun = NULL;
    cmdCancel = NULL;
    statusBar = NULL;
////@end wxBeginSimulation member initialisation
}


/*!
 * Control creation for wxBeginSimulation
 */

void wxBeginSimulation::CreateControls()
{    
////@begin wxBeginSimulation content construction
    wxBeginSimulation* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    txtSimLog = new wxTextCtrl( itemDialog1, ID_TXT_SIMLOG, _T(""), wxDefaultPosition, wxSize(450, 250), wxTE_MULTILINE );
    itemBoxSizer2->Add(txtSimLog, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxGROW|wxALL, 5);

    guageSimCompletion = new wxGauge( itemDialog1, ID_GAUGE, 100, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
    guageSimCompletion->SetValue(1);
    itemBoxSizer4->Add(guageSimCompletion, 1, wxGROW|wxALL, 5);

    cmdRun = new wxButton( itemDialog1, ID_CMD_RUN, _("Run"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(cmdRun, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    cmdCancel = new wxButton( itemDialog1, wxID_CANCEL, _("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(cmdCancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    statusBar = new wxStatusBar( itemDialog1, ID_STATUS_BAR, wxST_SIZEGRIP|wxNO_BORDER );
    statusBar->SetFieldsCount(1);
    itemBoxSizer2->Add(statusBar, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

////@end wxBeginSimulation content construction
}


/*!
 * Should we show tooltips?
 */

bool wxBeginSimulation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxBeginSimulation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxBeginSimulation bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxBeginSimulation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxBeginSimulation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxBeginSimulation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxBeginSimulation icon retrieval
}


bool wxBeginSimulation::Initialize(AppController *ctrl) {
	appController=ctrl;

	statusBar->SetStatusText(_("genomeSIMLA is ready to begin. Please click \"Run\" to start the simulation."));
 
	txtSimLog->Clear();
	stringstream ss;
	appController->SummarizeSimulation(ss);
	txtSimLog->WriteText(_T(ss.str().c_str()));

	//Set the range to include twice the number of drop points, 
	//since we will be ticking away for reports and drops
	guageSimCompletion->SetRange(ctrl->parameters.GetDropCount() * 2);
	guageSimCompletion->SetValue(0);

}

bool wxBeginSimulation::RunSimulation(bool doLoadFirst) {
	stringstream initstream;
	continueRunning = true;

	SetCursor(wxCursor(wxCURSOR_WAIT));
	appController->ClearCurrentLogEntries();
	size_t startingGeneration = appController->InitializeExecution(initstream);
	size_t dropPoints = appController->parameters.GetDropCount();
	txtSimLog->WriteText(_T(initstream.str().c_str()));
	int dropsMade = 0;
	bool simWasRun = false;
//	int reps = 0;
	while (dropPoints-- > 0 && continueRunning) {
		bool doAnalysis = appController->RunSimulation(doLoadFirst);
		do {
			txtSimLog->WriteText(_T(appController->GetSimOutput().c_str()));
			Update();
			wxYield();
			wxMilliSleep(50);
		} while (appController->SimIsRunning());

		guageSimCompletion->SetValue(++dropsMade);
		Update();
			
		if (doAnalysis && continueRunning)
			appController->PerformAnalysis();
		do {
			txtSimLog->WriteText(_T(appController->GetSimOutput().c_str()));
			Update();
			wxYield();
			wxMilliSleep(50);
		} while (appController->SimIsRunning());


		guageSimCompletion->SetValue(++dropsMade);
		Update();
		simWasRun = true;
	}	
	SetCursor(wxCursor(wxCURSOR_ARROW));

	if (!continueRunning)
		EndModal(0);
	else {
		isCompleted = true;
		txtSimLog->WriteText(_T("\nSimulation Completed. \n"));
	}
	return simWasRun;
}



/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_RUN
 */

void wxBeginSimulation::OnCmdRunClick( wxCommandEvent& event )
{
	if (isCompleted) 
		EndModal(1);
	else {
		if (appController)  {
			cmdRun->Enable(false);
			statusBar->PushStatusText(_T("genomeSIMLA is running. Please be patient while the simulation completes."));
			RunSimulation(true);
			statusBar->PopStatusText();
			
		}
		cmdRun->Enable(true);
		cmdRun->SetLabel(_T("OK"));
	}
}





/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON1
 */

void wxBeginSimulation::OnCancelClick( wxCommandEvent& event )
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
}

}



