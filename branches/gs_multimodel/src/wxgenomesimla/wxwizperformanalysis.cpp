/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizperformanalysis.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 27 Feb 2008 15:42:32 CST
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
#include "wxwizpageselectgeneration.h"
#include "wxWizPageMonitorAnalysis.h"
////@end includes

#include "wxwizperformanalysis.h"

////@begin XPM images

#include "img/GlowingStrand4.xpm"
////@end XPM images


/*!
 * wxWizPerformAnalysis type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPerformAnalysis, wxWizard )


/*!
 * wxWizPerformAnalysis event table definition
 */

BEGIN_EVENT_TABLE( wxWizPerformAnalysis, wxWizard )

////@begin wxWizPerformAnalysis event table entries
////@end wxWizPerformAnalysis event table entries

END_EVENT_TABLE()


/*!
 * wxWizPerformAnalysis constructors
 */

wxWizPerformAnalysis::wxWizPerformAnalysis()
{
    Init();
}

wxWizPerformAnalysis::wxWizPerformAnalysis( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
    Init();
    Create(parent, id, pos);
}


/*!
 * wxWizPerformAnalysis creator
 */

bool wxWizPerformAnalysis::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
////@begin wxWizPerformAnalysis creation
    SetExtraStyle(wxWIZARD_EX_HELPBUTTON);
    wxBitmap wizardBitmap(GetBitmapResource(wxT("img/GlowingStrand4.xpm")));
    wxWizard::Create( parent, id, _("Perform Analysis"), wizardBitmap, pos, wxDEFAULT_DIALOG_STYLE|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX );

    CreateControls();
////@end wxWizPerformAnalysis creation
    return true;
}


/*!
 * wxWizPerformAnalysis destructor
 */

wxWizPerformAnalysis::~wxWizPerformAnalysis()
{
////@begin wxWizPerformAnalysis destruction
////@end wxWizPerformAnalysis destruction
}


/*!
 * Member initialisation
 */

void wxWizPerformAnalysis::Init()
{
////@begin wxWizPerformAnalysis member initialisation
    pageSelectGeneration = NULL;
    pageMonitor = NULL;
////@end wxWizPerformAnalysis member initialisation
}


/*!
 * Control creation for wxWizPerformAnalysis
 */

void wxWizPerformAnalysis::CreateControls()
{    
////@begin wxWizPerformAnalysis content construction
    wxWizPerformAnalysis* itemWizard1 = this;

    pageSelectGeneration = new wxWizPageSelectGeneration( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageSelectGeneration);

    pageMonitor = new wxWizPageMonitorAnalysis( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageMonitor);

    wxWizardPageSimple* lastPage = NULL;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageSelectGeneration);
    lastPage = pageSelectGeneration;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageMonitor);
    lastPage = pageMonitor;
////@end wxWizPerformAnalysis content construction
}


/*!
 * Runs the wizard.
 */

bool wxWizPerformAnalysis::Run()
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

bool wxWizPerformAnalysis::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPerformAnalysis::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPerformAnalysis bitmap retrieval
    wxUnusedVar(name);
    if (name == _("img/GlowingStrand4.xpm"))
    {
        wxBitmap bitmap( GlowingStrand4_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxWizPerformAnalysis bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPerformAnalysis::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPerformAnalysis icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPerformAnalysis icon retrieval
}
