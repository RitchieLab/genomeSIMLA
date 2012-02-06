/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizardsimulatedata.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 10:39:25 CDT
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
#include "wxwizpagedefinedatasets.h"
#include "wxwizpageselectloci.h"
#include "wxwizpagereviewsimulation.h"
////@end includes

#include "wxwizardsimulatedata.h"

////@begin XPM images
#include "img/DnaSketch.xpm"
////@end XPM images


namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizardSimulateData type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizardSimulateData, wxWizard )


/*!
 * wxWizardSimulateData event table definition
 */

BEGIN_EVENT_TABLE( wxWizardSimulateData, wxWizard )

////@begin wxWizardSimulateData event table entries
////@end wxWizardSimulateData event table entries

END_EVENT_TABLE()


/*!
 * wxWizardSimulateData constructors
 */

wxWizardSimulateData::wxWizardSimulateData()
{
    Init();
}

wxWizardSimulateData::wxWizardSimulateData( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
    Init();
    Create(parent, id, pos);
}


/*!
 * wxWizardSimulateData creator
 */

bool wxWizardSimulateData::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
////@begin wxWizardSimulateData creation
    SetExtraStyle(wxWIZARD_EX_HELPBUTTON);
    wxBitmap wizardBitmap(GetBitmapResource(wxT("img/DnaSketch.xpm")));
    wxWizard::Create( parent, id, _("Simulate Data"), wizardBitmap, pos, wxDEFAULT_DIALOG_STYLE|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX );

    CreateControls();
////@end wxWizardSimulateData creation
    return true;
}


/*!
 * wxWizardSimulateData destructor
 */

wxWizardSimulateData::~wxWizardSimulateData()
{
////@begin wxWizardSimulateData destruction
////@end wxWizardSimulateData destruction
}


/*!
 * Member initialisation
 */

void wxWizardSimulateData::Init()
{
////@begin wxWizardSimulateData member initialisation
    pageDefineDatasets = NULL;
    pageSelectLoci = NULL;
    pageSimulationReview = NULL;
////@end wxWizardSimulateData member initialisation
}


/*!
 * Control creation for wxWizardSimulateData
 */

void wxWizardSimulateData::CreateControls()
{    
////@begin wxWizardSimulateData content construction
    wxWizardSimulateData* itemWizard1 = this;

    pageDefineDatasets = new wxWizPageDefineDatasets( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageDefineDatasets);

    pageSelectLoci = new wxWizPageSelectLoci( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageSelectLoci);

    pageSimulationReview = new wxWizPageReviewSimulation( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageSimulationReview);

    wxWizardPageSimple* lastPage = NULL;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageDefineDatasets);
    lastPage = pageDefineDatasets;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageSelectLoci);
    lastPage = pageSelectLoci;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageSimulationReview);
    lastPage = pageSimulationReview;
////@end wxWizardSimulateData content construction
}

void wxWizardSimulateData::InitAppController(AppController *appController) {
	pageDefineDatasets->InitAppController(appController);
	pageSelectLoci->InitAppController(appController);
	pageSimulationReview->InitAppController(appController);
}

bool  wxWizardSimulateData::SetChromosomes(vector<FileBasedChromosome *> *ch){  
	bool success=false;
	chromosomes=ch;

	if (pageSelectLoci)
		success = pageSelectLoci->Initialize(chromosomes);
	return success;
}

/*!
 * Runs the wizard.
 */

bool wxWizardSimulateData::Run()
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

bool wxWizardSimulateData::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizardSimulateData::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizardSimulateData bitmap retrieval
    wxUnusedVar(name);
    if (name == _T("img/DnaSketch.xpm"))
    {
        wxBitmap bitmap( DnaSketch_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxWizardSimulateData bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizardSimulateData::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizardSimulateData icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizardSimulateData icon retrieval
}

}

}
