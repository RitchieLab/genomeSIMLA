/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagereviewsimulation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 28 Mar 2008 11:10:59 CDT
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

#include "wxwizpagereviewsimulation.h"
#include "wxdlgprocessmonitor.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageReviewSimulation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageReviewSimulation, wxWizardPageDataSummary )


/*!
 * wxWizPageReviewSimulation event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageReviewSimulation, wxWizardPageDataSummary )

////@begin wxWizPageReviewSimulation event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageReviewSimulation::OnPageSimulationReviewPageChanged )
    EVT_SIZE( wxWizPageReviewSimulation::OnSize )

    EVT_BUTTON( ID_CMD_GENERATE_DATA, wxWizPageReviewSimulation::OnCmdGenerateDataClick )

////@end wxWizPageReviewSimulation event table entries
	EVT_BUTTON(ID_SAVE_PENETRANCE, wxWizPageReviewSimulation::OnSavePenetrance )

END_EVENT_TABLE()


/*!
 * wxWizPageReviewSimulation constructors
 */

wxWizPageReviewSimulation::wxWizPageReviewSimulation()
{
    Init();
}

wxWizPageReviewSimulation::wxWizPageReviewSimulation( wxWizard* parent )
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageReviewSimulation creator
 */

bool wxWizPageReviewSimulation::Create( wxWizard* parent )
{
////@begin wxWizPageReviewSimulation creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageDataSummary::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageReviewSimulation creation
    return true;
}


/*!
 * wxWizPageReviewSimulation destructor
 */

wxWizPageReviewSimulation::~wxWizPageReviewSimulation()
{
////@begin wxWizPageReviewSimulation destruction
////@end wxWizPageReviewSimulation destruction
}


/*!
 * Member initialisation
 */

void wxWizPageReviewSimulation::Init()
{
	wxSavePenetrance = NULL;
////@begin wxWizPageReviewSimulation member initialisation
    txtReview = NULL;
////@end wxWizPageReviewSimulation member initialisation
}


/*!
 * Control creation for wxWizPageReviewSimulation
 */

void wxWizPageReviewSimulation::CreateControls()
{    
////@begin wxWizPageReviewSimulation content construction
    wxWizPageReviewSimulation* itemWizardPageDataSummary1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageDataSummary1->SetSizer(itemBoxSizer2);

    txtReview = new wxHtmlWindow( itemWizardPageDataSummary1, ID_HTMLWINDOW2, wxDefaultPosition, wxSize(200, 150), wxHW_SCROLLBAR_AUTO|wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    itemBoxSizer2->Add(txtReview, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer4, 0, wxALIGN_RIGHT|wxALL, 5);
////@end wxWizPageReviewSimulation content construction

	wxSavePenetrance = new wxButton( itemWizardPageDataSummary1, ID_SAVE_PENETRANCE , _("Save Penetrance"), 
		wxDefaultPosition, wxDefaultSize, 0);
	itemBoxSizer4->Add(wxSavePenetrance, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton5 = new wxButton( itemWizardPageDataSummary1, ID_CMD_GENERATE_DATA, _("Generate Data"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemButton5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);


}


/*!
 * Should we show tooltips?
 */

bool wxWizPageReviewSimulation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageReviewSimulation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageReviewSimulation bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageReviewSimulation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageReviewSimulation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageReviewSimulation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageReviewSimulation icon retrieval
}

void wxWizPageReviewSimulation::Clear(const char *title) {
	txtReview->SetPage("<HTML><HEAD><TITLE>");
	txtReview->AppendToPage(title);
	txtReview->AppendToPage("</TITLE></HEAD>");
}

void wxWizPageReviewSimulation::AddHeader(const char *txt) {
	txtReview->AppendToPage("<H3>");
	txtReview->AppendToPage(txt);
	txtReview->AppendToPage("</H3>");
}

void wxWizPageReviewSimulation::AddNote(const char *txt) {
	txtReview->AppendToPage(txt);
	txtReview->AppendToPage("<P>");
}

void wxWizPageReviewSimulation::AddConfigLine(const char *txt) {
	configuration<<txt<<"\n";
}

/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_PAGE_SIMULATION_REVIEW
 */

void wxWizPageReviewSimulation::OnPageSimulationReviewPageChanged( wxWizardEvent& event )
{
	((wxWizardPageDataSim*)GetPrev())->PrepSummary(this);
	wxSavePenetrance->Enable(diseaseModel.modelType == "SIMLA" || diseaseModel.modelType == "SIMPEN");
}


/*!
 * wxEVT_SIZE event handler for ID_PAGE_SIMULATION_REVIEW
 */

void wxWizPageReviewSimulation::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_PAGE_SIMULATION_REVIEW in wxWizPageReviewSimulation.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_PAGE_SIMULATION_REVIEW in wxWizPageReviewSimulation. 
}
/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_GENERATE_DATA
 */

void wxWizPageReviewSimulation::OnCmdGenerateDataClick( wxCommandEvent& event )
{
	stringstream ss;
	((wxWizardPageDataSim*)GetPrev())->PrepConfig(ss);
	cout<<ss.str()<<"\n";

	wxDlgDatasetGenerationMonitor dlg(this, appController);
	cout<<"Showing dataset generation dialog\n";
	dlg.ShowModal();

}

void wxWizPageReviewSimulation::OnSavePenetrance( wxCommandEvent& event) {
	string penFilename;
	wxFileDialog newPenFile(this, _("Save Penetrance File"), wxT(""), _(""), wxT("Penetrance Table (*.pen)|*.pen"), wxSAVE|wxOVERWRITE_PROMPT);

	if (newPenFile.ShowModal() == wxID_OK) {
		wxFileName filename(newPenFile.GetPath());
		AppConfig *config = AppConfig::GetPrimaryInstance();
		cout<<"Saving penetrance file: "<<filename.GetFullPath().c_str()<<"\n";
	
		PenetranceModel *model = config->statusSettings.model;
		model->WritePenetranceFile(filename.GetFullPath().c_str(), loci, 0.001);
	}
}



}

}




