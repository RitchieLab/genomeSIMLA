/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgselectlocusrange.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Wed 12 Dec 2007 15:13:17 CST
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

#include "wxdlgselectlocusrange.h"
#include <iostream>

////@begin XPM images
////@end XPM images
namespace GenomeSIM {

namespace GUI {

using namespace std;


/*!
 * wxDlgSelectLocusRange type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgSelectLocusRange, wxDialog )


/*!
 * wxDlgSelectLocusRange event table definition
 */

BEGIN_EVENT_TABLE( wxDlgSelectLocusRange, wxDialog )

////@begin wxDlgSelectLocusRange event table entries
////@end wxDlgSelectLocusRange event table entries
	EVT_LIST_ITEM_SELECTED( ID_LOCUSLISTCONTROL, wxDlgSelectLocusRange::OnLocuslistcontrolSelected )
	EVT_LIST_ITEM_ACTIVATED( ID_LOCUSLISTCONTROL, wxDlgSelectLocusRange::OnLocuslistcontrolItemActivated )

END_EVENT_TABLE()


/*!
 * wxDlgSelectLocusRange constructors
 */

wxDlgSelectLocusRange::wxDlgSelectLocusRange(): chrom(NULL), curSelection(NULL) {
    Init();
}

wxDlgSelectLocusRange::wxDlgSelectLocusRange( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) : chrom(NULL), curSelection(NULL) {
    Init();
    Create(parent, id, caption, pos, size, style);
}


/*!
 * wxDlgSelectLocusRange creator
 */

bool wxDlgSelectLocusRange::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style ) 
{
////@begin wxDlgSelectLocusRange creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
    wxDialog::Create( parent, id, caption, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgSelectLocusRange creation
    return true;
}


/*!
 * wxDlgSelectLocusRange destructor
 */

wxDlgSelectLocusRange::~wxDlgSelectLocusRange()
{
////@begin wxDlgSelectLocusRange destruction
////@end wxDlgSelectLocusRange destruction
}


void wxDlgSelectLocusRange::Initialize(FileBasedChromosome *chrom, Locus *loc) {
	this->chrom = chrom;
	curSelection = loc;
	int idx = chrom->GetLocusIdx(loc->GetLabel().c_str());

	txtChromosomeLabel->SetLabel(chrom->label.c_str());


	lstLoci->Initialize(&(chrom->GetLoci()), idx);

	RefreshSize();
	
	lstLoci->SetItemState(idx, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED);
	lstLoci->EnsureVisible(idx);
}

void wxDlgSelectLocusRange::RefreshSize() {
	int width, height;

	lstLoci->GetSize(&width, &height);
	width-=65;
	width/=4;
	
	lstLoci->SetColumnWidth(0, 50);
	for (int i=1; i<5; i++) 
		lstLoci->SetColumnWidth(i, width);

	//lstLoci->FindItem(1, wxT(curSelection->GetLabel().c_str()));
	//


}

void wxDlgSelectLocusRange::EndOK() {
	if (Validate() && TransferDataFromWindow() ) {
		if ( IsModal() )
			EndModal(wxID_OK);
		else {
			SetReturnCode(wxID_OK);
			this->Show(false);
		}
	}
}


/*!
 * Member initialisation
 */

void wxDlgSelectLocusRange::Init()
{
////@begin wxDlgSelectLocusRange member initialisation
    txtChromosomeLabel = NULL;
    lstLoci = NULL;
////@end wxDlgSelectLocusRange member initialisation
}


/*!
 * Control creation for wxDlgSelectLocusRange
 */

void wxDlgSelectLocusRange::CreateControls()
{    
////@begin wxDlgSelectLocusRange content construction
    wxDlgSelectLocusRange* itemDialog1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);

    txtChromosomeLabel = new wxStaticText( itemDialog1, wxID_STATIC, _("Chromosome-Label"), wxDefaultPosition, wxDefaultSize, 0 );
    txtChromosomeLabel->SetFont(wxFont(12, wxSWISS, wxNORMAL, wxBOLD, false, wxT("Sans")));
    itemBoxSizer2->Add(txtChromosomeLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    lstLoci = new wxLocusListControl( itemDialog1, ID_LOCUSLISTCONTROL, wxDefaultPosition, wxSize(500, 250), wxLC_REPORT|wxLC_VIRTUAL  );
    itemBoxSizer2->Add(lstLoci, 1, wxGROW|wxALL, 5);

    // Connect events and objects
    lstLoci->Connect(ID_LOCUSLISTCONTROL, wxEVT_SIZE, wxSizeEventHandler(wxDlgSelectLocusRange::OnSize), NULL, this);
////@end wxDlgSelectLocusRange content construction
 
}




Locus *wxDlgSelectLocusRange::GetSelection() {
	return curSelection;
}

/*!
 * Should we show tooltips?
 */

bool wxDlgSelectLocusRange::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgSelectLocusRange::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgSelectLocusRange bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgSelectLocusRange bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgSelectLocusRange::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgSelectLocusRange icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgSelectLocusRange icon retrieval
}









/*!
 * wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LOCUSLISTCONTROL
 */

void wxDlgSelectLocusRange::OnLocuslistcontrolSelected( wxListEvent& event )
{
	curSelection = chrom->GetLocus(event.GetIndex());
//	EndOK();
////@begin wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange. 
}





/*!
 * wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LOCUSLISTCONTROL
 */

void wxDlgSelectLocusRange::OnLocuslistcontrolItemActivated( wxListEvent& event )
{
	curSelection = chrom->GetLocus(event.GetIndex());

	EndOK();

////@begin wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange. 
}



/*!
 * wxEVT_SIZE event handler for ID_LOCUSLISTCONTROL
 */

void wxDlgSelectLocusRange::OnSize( wxSizeEvent& event )
{
	RefreshSize();
////@begin wxEVT_SIZE event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_LOCUSLISTCONTROL in wxDlgSelectLocusRange. 
}

}

}








