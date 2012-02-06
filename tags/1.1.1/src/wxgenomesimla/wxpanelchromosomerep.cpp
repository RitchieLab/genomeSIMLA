/////////////////////////////////////////////////////////////////////////////
// Name:        wxpanelchromosomerep.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 13:27:12 CDT
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

#include "wxpanelchromosomerep.h"

////@begin XPM images
////@end XPM images
#define wx_MAC_USE_CORE_GRAPHICS 1
namespace GenomeSIM {

namespace GUI {

/*!
 * wxPanelChromosomeRep type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPanelChromosomeRep, wxPanel )


/*!
 * wxPanelChromosomeRep event table definition
 */

BEGIN_EVENT_TABLE( wxPanelChromosomeRep, wxPanel )

////@begin wxPanelChromosomeRep event table entries
    EVT_SIZE( wxPanelChromosomeRep::OnSize )
    EVT_PAINT( wxPanelChromosomeRep::OnPaint )
    EVT_MOTION( wxPanelChromosomeRep::OnMotion )

////@end wxPanelChromosomeRep event table entries

END_EVENT_TABLE()


/*!
 * wxPanelChromosomeRep constructors
 */

wxPanelChromosomeRep::wxPanelChromosomeRep()
{
    Init();
}

wxPanelChromosomeRep::wxPanelChromosomeRep(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style)
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPanelChromosomeRep creator
 */

bool wxPanelChromosomeRep::Create(wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style)
{
////@begin wxPanelChromosomeRep creation
    wxPanel::Create(parent, id, pos, size, style);
    CreateControls();
////@end wxPanelChromosomeRep creation
    return true;
}


/*!
 * wxPanelChromosomeRep destructor
 */

wxPanelChromosomeRep::~wxPanelChromosomeRep()
{
////@begin wxPanelChromosomeRep destruction
////@end wxPanelChromosomeRep destruction
}


/*!
 * Member initialisation
 */

void wxPanelChromosomeRep::Init()
{
////@begin wxPanelChromosomeRep member initialisation
////@end wxPanelChromosomeRep member initialisation
}


/*!
 * Control creation for wxPanelChromosomeRep
 */

void wxPanelChromosomeRep::CreateControls()
{    
////@begin wxPanelChromosomeRep content construction
    if (ShowToolTips())
        this->SetToolTip(_("Click to Regenerate the Chromosome"));
////@end wxPanelChromosomeRep content construction
}


/*!
 * Should we show tooltips?
 */

bool wxPanelChromosomeRep::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPanelChromosomeRep::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPanelChromosomeRep bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPanelChromosomeRep bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPanelChromosomeRep::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPanelChromosomeRep icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPanelChromosomeRep icon retrieval
}

/*!
 * wxEVT_PAINT event handler for ID_CHROMOSOME_VISUAL
 */

void wxPanelChromosomeRep::OnPaint( wxPaintEvent& event )
{
	int lociCount = loci.size();
	wxCoord w, h;
    wxPaintDC dc(this);
	dc.GetSize(&w, &h);
	int x = 1, y = 1;
	int width = w -1,
		height = h - 1;

//	cout<<"Ends: "<<x<<" "<<y<<" - "<<width<<" "<<height<<"\n";

	dc.Clear();
	dc.DrawRectangle(x, y, width, height);
//	cout<<"Draw Rect: "<<x<<", "<<y<<", "<<width<<", "<<height<<"\n";

	if (lociCount > 1) {
		uint start = loci[0].GetLocation();
		uint stop  = loci[lociCount-1].GetLocation();
		float breadth = stop-start;
//		cout<<"Start: "<<start<<"\tStop: "<<stop<<"\tBreadth: "<<breadth<<"\n";

		for (int i=1; i<lociCount-1; i++) {
			uint cur = loci[i].GetLocation() - start;
			//pos is ranking from 0 to 1.0
			float pos = (float)cur / breadth;
			dc.DrawRectangle((int)(pos * width + x), y, 1, height);
		}

		wxString label;
		wxCoord wTxt, hTxt;
		label.Printf(_("%d Loci : %d kb"), lociCount, (stop-start)/1000);
		wxFont font(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD);
		dc.SetFont(font);
		dc.GetTextExtent(label, &wTxt, &hTxt);
		
		dc.SetBrush(wxBrush(wxColour(226, 226, 226), wxSOLID));
//		cout<<"Draw Rounded Rect: "<<(w-wTxt)/2-5<<","<<(h-hTxt)-2<<", "<<wTxt+10<<", "<<h<<", 10\n";
		dc.DrawRoundedRectangle((w-wTxt)/2-5, (h-hTxt)-2, wTxt+10, h, 10);
			

		//dc.SetBackgroundMode(wxSOLID);
		dc.SetTextBackground(wxColour(226, 226, 255, 16));
		dc.SetTextForeground(wxColour(255, 255, 255));
		
		dc.DrawText(label, (w-wTxt)/2-1, (h-hTxt)-1);
//		cout<<"Draw Text: "<<label<<" "<<(w-wTxt)/2-1<<", "<<(h-hTxt)-1<<"\n";
		//dc.SetBackgroundMode(wxTRANSPARENT);
			
		dc.SetTextForeground(wxColour(76, 86, 80));
		dc.DrawText(label, (w-wTxt)/2, (h-hTxt));
	}
}

void wxPanelChromosomeRep::Refresh_Image() {
	Refresh();
	Update();
}


/*!
 * wxEVT_MOTION event handler for ID_CHROMOSOME_VISUAL
 */

void wxPanelChromosomeRep::OnMotion( wxMouseEvent& event )
{
////@begin wxEVT_MOTION event handler for ID_CHROMOSOME_VISUAL in wxPanelChromosomeRep.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_MOTION event handler for ID_CHROMOSOME_VISUAL in wxPanelChromosomeRep. 
}

/*!
 * wxEVT_SIZE event handler for ID_CHROMOSOME_VISUAL
 */

void wxPanelChromosomeRep::OnSize( wxSizeEvent& event )
{
	Refresh();
////@begin wxEVT_SIZE event handler for ID_CHROMOSOME_VISUAL in wxPanelChromosomeRep.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_CHROMOSOME_VISUAL in wxPanelChromosomeRep. 
}

void wxPanelChromosomeRep::SetLoci(vector<Locus>& loci) {
	this->loci = loci;
	
	Refresh();
}

}

}







