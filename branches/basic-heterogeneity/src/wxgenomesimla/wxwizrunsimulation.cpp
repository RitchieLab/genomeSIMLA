/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizrunsimulation.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 29 Jan 2008 11:51:23 CST
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
#include "wxwizpagelaunchproject.h"
#include "wxwizpagelaunchselectgeneration.h"
#include "wxwizpagelaunchdroppoints.h"
#include "wxwizpagelaunchconfirmation.h"
#include "wxwizpagereview.h"
////@end includes
#include <wx/filename.h>
#include "wxwizrunsimulation.h"
#include "utility/strings.h"

////@begin XPM images
#include "img/DnaSketch.xpm"
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizRunSimulation type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizRunSimulation, wxWizard )


/*!
 * wxWizRunSimulation event table definition
 */

BEGIN_EVENT_TABLE( wxWizRunSimulation, wxWizard )

////@begin wxWizRunSimulation event table entries
////@end wxWizRunSimulation event table entries

END_EVENT_TABLE()


/*!
 * wxWizRunSimulation constructors
 */

wxWizRunSimulation::wxWizRunSimulation() 
{
    Init();
}

wxWizRunSimulation::wxWizRunSimulation( wxWindow* parent, bool skipGenSelection, AppController *appCtrl, wxWindowID id, const wxPoint& pos ) 
{
	appController = appCtrl;
    Init();
    Create(parent, skipGenSelection, id, pos);
}


/*!
 * wxWizRunSimulation creator
 */

bool wxWizRunSimulation::Create( wxWindow* parent, bool skipGenSelection, wxWindowID id, const wxPoint& pos )
{
////@begin wxWizRunSimulation creation
    SetExtraStyle(wxWIZARD_EX_HELPBUTTON);
    wxBitmap wizardBitmap(GetBitmapResource(wxT("img/DnaSketch.xpm")));
    wxWizard::Create( parent, id, _("Run Simulation"), wizardBitmap, pos, wxDEFAULT_DIALOG_STYLE|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX );

    CreateControls(skipGenSelection);
////@end wxWizRunSimulation creation
    return true;
}


/*!
 * wxWizRunSimulation destructor
 */

wxWizRunSimulation::~wxWizRunSimulation()
{
////@begin wxWizRunSimulation destruction
////@end wxWizRunSimulation destruction
}


/*!
 * Member initialisation
 */

void wxWizRunSimulation::Init()
{
////@begin wxWizRunSimulation member initialisation
////@end wxWizRunSimulation member initialisation
	pageProjectSettings = NULL;
//	pageSelectPreviousRun = NULL;
	pageDropPoints = NULL;
	pageConfirmation = NULL;
	pageReview = NULL;
	
}


/*!
 * Control creation for wxWizRunSimulation
 */

void wxWizRunSimulation::CreateControls(bool skipGenSelection)
{    
////@ - begin wxWiz - RunSimulation content construction
    wxWizard* itemWizard1 = this;

	if (!skipGenSelection) {
		pageProjectSettings = new wxWizPageLaunchProject( itemWizard1 );
		itemWizard1->GetPageAreaSizer()->Add(pageProjectSettings);
	
//		This page is redundant, since we select the starting generation from the report page
//		 = new wxWizPageLaunchSelectGeneration( itemWizard1, appController );
//		itemWizard1->GetPageAreaSizer()->Add(pageSelectPreviousRun);
	}
    pageDropPoints = new wxWizPageLaunchDropPoints( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageDropPoints);

    pageConfirmation = new wxWizPageLaunchConfirmation( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageConfirmation);

    pageReview = new wxWizPageReview( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageReview);

    wxWizardPageSimple* lastPage = NULL;

	if (!skipGenSelection) {
		if (lastPage)
			wxWizardPageSimple::Chain(lastPage, pageProjectSettings);
		lastPage = pageProjectSettings;
//		if (lastPage)
//			wxWizardPageSimple::Chain(lastPage, pageSelectPreviousRun);
//		lastPage = pageSelectPreviousRun;
	}
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageDropPoints);
    lastPage = pageDropPoints;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageConfirmation);
    lastPage = pageConfirmation;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageReview);
    lastPage = pageReview;
//// @ end wx - WizRunSimulation content construction

	if (pageConfirmation)
		pageConfirmation->SetAppController(appController);

	if (pageDropPoints)
		pageDropPoints->SetAppController(appController);

	if (pageProjectSettings)
		pageProjectSettings->SetAppController(appController);

//	if (pageSelectPreviousRun)
//		pageSelectPreviousRun->SetAppController(appController);

	if (pageReview)
		pageReview->SetAppController(appController);
}



/*!
 * Runs the wizard.
 */

bool wxWizRunSimulation::Run()
{
    wxWindowList::compatibility_iterator node = GetChildren().GetFirst();
    while (node)
    {
        wxWizardPage* startPage = wxDynamicCast(node->GetData(), wxWizardPage);
        if (startPage) return RunWizard(startPage);
        node = node->GetNext();
    }
    return false;
}


/*!
 * Should we show tooltips?
 */

bool wxWizRunSimulation::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizRunSimulation::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizRunSimulation bitmap retrieval
    wxUnusedVar(name);
    if (name == _("img/DnaSketch.xpm"))
    {
        wxBitmap bitmap( DnaSketch_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxWizRunSimulation bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizRunSimulation::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizRunSimulation icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizRunSimulation icon retrieval
}



}

}




