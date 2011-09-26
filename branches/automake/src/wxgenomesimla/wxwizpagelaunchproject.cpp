/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagelaunchproject.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 01 Feb 2008 17:12:50 CST
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

#include "wxwizpagelaunchproject.h"
#include "wx/filename.h"
////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxWizPageLaunchProject type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageLaunchProject, wxWizardPageSimple )


/*!
 * wxWizPageLaunchProject event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageLaunchProject, wxWizardPageSimple )

////@begin wxWizPageLaunchProject event table entries
    EVT_WIZARD_PAGE_CHANGED( -1, wxWizPageLaunchProject::OnWizProjectSettingsPageChanged )
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageLaunchProject::OnWizProjectSettingsPageChanging )
    EVT_SIZE( wxWizPageLaunchProject::OnSize )

    EVT_COMBOBOX( ID_CMB_PROJECT_NAME, wxWizPageLaunchProject::OnCmbProjectNameSelected )

    EVT_BUTTON( ID_CMD_SELECT_PROJECT_PATH, wxWizPageLaunchProject::OnCmdSelectProjectPathClick )

////@end wxWizPageLaunchProject event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageLaunchProject constructors
 */

wxWizPageLaunchProject::wxWizPageLaunchProject() {
    Init();
}

wxWizPageLaunchProject::wxWizPageLaunchProject( wxWizard* parent ) {
    Init();
    Create( parent );
}


/*!
 * wxWizPageLaunchProject creator
 */

bool wxWizPageLaunchProject::Create( wxWizard* parent ) {
////@begin wxWizPageLaunchProject creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageLaunchProject creation
    return true;
}


/*!
 * wxWizPageLaunchProject destructor
 */

wxWizPageLaunchProject::~wxWizPageLaunchProject() {
////@begin wxWizPageLaunchProject destruction
////@end wxWizPageLaunchProject destruction
}


/*!
 * Member initialisation
 */

void wxWizPageLaunchProject::Init() {
////@begin wxWizPageLaunchProject member initialisation
    cmbProjectName = NULL;
    txtProjectPath = NULL;
////@end wxWizPageLaunchProject member initialisation
}


/*!
 * Control creation for wxWizPageLaunchProject
 */

void wxWizPageLaunchProject::CreateControls() {    
////@begin wxWizPageLaunchProject content construction
    wxWizPageLaunchProject* itemWizardPageSimple1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageSimple1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer4, 2, wxGROW|wxALL, 0);

    wxStaticText* itemStaticText5 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Project Settings:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemStaticText5, 0, wxALIGN_LEFT|wxALL, 5);

    wxStaticText* itemStaticText6 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Select the name of the project name for this execution. If the project already exists, you will be able to continue running from one of the previous stopping points or start over once again. \n\nThe project name will be used as part of the filenames generated during the execution. "), wxDefaultPosition, wxSize(450, -1), wxST_NO_AUTORESIZE );
    itemBoxSizer4->Add(itemStaticText6, 0, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer7 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer4->Add(itemBoxSizer7, 0, wxGROW|wxALL, 10);

    wxStaticText* itemStaticText8 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Project Name:"), wxDefaultPosition, wxSize(150, -1), wxALIGN_LEFT );
    itemBoxSizer7->Add(itemStaticText8, 0, wxALIGN_LEFT|wxALL, 5);

    wxBoxSizer* itemBoxSizer9 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer7->Add(itemBoxSizer9, 0, wxGROW|wxALL, 0);

    wxArrayString cmbProjectNameStrings;
    cmbProjectName = new wxComboBox( itemWizardPageSimple1, ID_CMB_PROJECT_NAME, _T(""), wxDefaultPosition, wxSize(200, -1), cmbProjectNameStrings, wxCB_DROPDOWN );
    itemBoxSizer9->Add(cmbProjectName, 1, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    itemBoxSizer4->Add(5, 25, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticText* itemStaticText12 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("The project path is where the files will be created. "), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer4->Add(itemStaticText12, 0, wxALIGN_LEFT|wxALL, 5);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer4->Add(itemBoxSizer13, 0, wxGROW|wxALL, 10);

    wxStaticText* itemStaticText14 = new wxStaticText( itemWizardPageSimple1, wxID_STATIC, _("Project Directory:"), wxDefaultPosition, wxSize(150, -1), wxALIGN_LEFT );
    itemBoxSizer13->Add(itemStaticText14, 0, wxALIGN_LEFT|wxALL, 5);

    wxBoxSizer* itemBoxSizer15 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer13->Add(itemBoxSizer15, 1, wxGROW|wxALL, 0);

    txtProjectPath = new wxTextCtrl( itemWizardPageSimple1, ID_TXT_PROJECT_PATH, _T(""), wxDefaultPosition, wxSize(200, -1), 0 );
    itemBoxSizer15->Add(txtProjectPath, 1, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxButton* itemButton17 = new wxButton( itemWizardPageSimple1, ID_CMD_SELECT_PROJECT_PATH, _("..."), wxDefaultPosition, wxSize(25, -1), 0 );
    itemBoxSizer15->Add(itemButton17, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    itemBoxSizer4->Add(5, 25, 1, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    itemBoxSizer2->Add(0, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

////@end wxWizPageLaunchProject content construction

	wxString curPath = wxGetCwd();
	SetProjectDir(curPath.c_str());
}

void wxWizPageLaunchProject::SetProjectDir(const char *filename) {
	wxFileName relFilename (wxString(_(filename)));	
	relFilename.MakeRelativeTo(wxGetCwd());

	txtProjectPath->ChangeValue(relFilename.GetPath().c_str());

}

/*!
 * wxEVT_WIZARD_PAGE_CHANGED event handler for ID_WIZ_PROJECT_SETTINGS
 */

void wxWizPageLaunchProject::OnWizProjectSettingsPageChanged( wxWizardEvent& event ) {
	if (event.GetDirection()) {
		string filename;
		
		if (appController->parameters.GetProjectName() == "")
			filename = appController->parameters.GetConfigFilename();
		else
			filename = appController->parameters.GetProjectName();
		string file, path, ext;
		vector<string> projects;

		appController->GetProjectList(projects);

		vector<string>::iterator itr = projects.begin();
		vector<string>::iterator end = projects.end();
		while (itr != end) {
			wxFileName filename (_(itr->c_str()));
			cmbProjectName->Insert(filename.GetFullName(),projectFilenames.size());
			projectFilenames.push_back(itr->c_str());
			
			itr++;
		}
		if (projects.size() > 0) {
			cmbProjectName->SetSelection(0);
			SetProjectDir(cmbProjectName->GetValue().c_str());
		}
		else {
			SetProjectDir(filename.c_str());
			cmbProjectName->SetValue(ExtractFilename(filename.c_str()).c_str());
		}
	}
}


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_WIZ_PROJECT_SETTINGS
 */
void wxWizPageLaunchProject::OnWizProjectSettingsPageChanging( wxWizardEvent& event ) {
	if (event.GetPage() == this && event.GetDirection()) {
		wxFileName projectFilename(txtProjectPath->GetValue(), cmbProjectName->GetValue());
		projectFilename.MakeRelativeTo(wxGetCwd());

		appController->parameters.SetProjectName(projectFilename.GetFullPath().c_str());
		if (cmbProjectName->GetSelection() < 0)
			appController->parameters.SetStartingGeneration(0);
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_SELECT_PROJECT_PATH
 */
void wxWizPageLaunchProject::OnCmdSelectProjectPathClick( wxCommandEvent& event ) {
	wxString curPath = txtProjectPath->GetValue();
	wxFileName projectFilename;
	
	if (wxDirExists(curPath)) 
		projectFilename.AssignDir(curPath);
	else
		projectFilename.Assign(curPath);
	
	wxDirDialog dialog(this, wxT("Select the project directory"), projectFilename.GetPath(), wxDD_NEW_DIR_BUTTON);
	
	if (dialog.ShowModal() == wxID_OK) {
		txtProjectPath->SetValue(dialog.GetPath());
	}
}


/*!
 * Should we show tooltips?
 */

bool wxWizPageLaunchProject::ShowToolTips() {
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageLaunchProject::GetBitmapResource( const wxString& name ) {
    // Bitmap retrieval
////@begin wxWizPageLaunchProject bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageLaunchProject bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageLaunchProject::GetIconResource( const wxString& name ){
    // Icon retrieval
////@begin wxWizPageLaunchProject icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageLaunchProject icon retrieval
}

wxWizardPage *wxWizPageLaunchProject::GetNext() const{
	wxWizardPage *next = wxWizardPageSimple::GetNext();
/*	if (next && cmbProjectName->GetSelection() < 0)
		return next->GetNext();
	else
*/		return next;	
}



/*!
 * wxEVT_COMMAND_COMBOBOX_SELECTED event handler for ID_CMB_PROJECT_NAME
 */

void wxWizPageLaunchProject::OnCmbProjectNameSelected( wxCommandEvent& event )
{
	SetProjectDir(projectFilenames[cmbProjectName->GetSelection()].c_str());
}



/*!
 * wxEVT_SIZE event handler for ID_WIZ_PROJECT_SETTINGS
 */

void wxWizPageLaunchProject::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_WIZ_PROJECT_SETTINGS in wxWizPageLaunchProject.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_WIZ_PROJECT_SETTINGS in wxWizPageLaunchProject. 
}

}

}


