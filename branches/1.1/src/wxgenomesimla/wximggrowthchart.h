//
// C++ Interface: wximggrowthchart
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXIMGGROWTHCHART_H
#define GENOMESIM_GUIWXIMGGROWTHCHART_H

#include <wx/panel.h>
#include <wx/bitmap.h>
#include <string>
#include "appinterface.h"

#define ID_GROWTH_CHART 10006
#define SYMBOL_WXIMGGROWTHCHART_STYLE wxSUNKEN_BORDER|wxTAB_TRAVERSAL
#define SYMBOL_WXIMGGROWTHCHART_IDNAME ID_GROWTH_CHART
#define SYMBOL_WXIMGGROWTHCHART_SIZE wxSize(300, 125)
#define SYMBOL_WXIMGGROWTHCHART_POSITION wxDefaultPosition


namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxImgGrowthChart class declaration
 */

class wxImgGrowthChart: public wxPanel, public AppInterface 	{    
    DECLARE_DYNAMIC_CLASS( wxImgGrowthChart )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxImgGrowthChart( );
    wxImgGrowthChart( wxWindow* parent, wxWindowID id = SYMBOL_WXIMGGROWTHCHART_IDNAME, const wxPoint& pos = SYMBOL_WXIMGGROWTHCHART_POSITION, const wxSize& size = SYMBOL_WXIMGGROWTHCHART_SIZE, long style = SYMBOL_WXIMGGROWTHCHART_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXIMGGROWTHCHART_IDNAME, const wxPoint& pos = SYMBOL_WXIMGGROWTHCHART_POSITION, const wxSize& size = SYMBOL_WXIMGGROWTHCHART_SIZE, long style = SYMBOL_WXIMGGROWTHCHART_STYLE );

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxImgGrowthChart event handler declarations

    /// wxEVT_SIZE event handler for ID_GROWTH_CHART
    void OnSize( wxSizeEvent& event );

    /// wxEVT_PAINT event handler for ID_GROWTH_CHART
    void OnPaint( wxPaintEvent& event );

    /// wxEVT_ERASE_BACKGROUND event handler for ID_GROWTH_CHART
    void OnEraseBackground( wxEraseEvent& event );

	void OnLeftDClick(wxMouseEvent& event);

////@end wxImgGrowthChart event handler declarations

////@begin wxImgGrowthChart member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxImgGrowthChart member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxImgGrowthChart member variables
////@end wxImgGrowthChart member variables
	
	void InitAppController(AppController *appController);
	void Commit();
	void RefreshSettings();
	

	void SetImagePng(const char *png);
	void RefreshImage();
	void SetModel(const char *modelConfiguration);
	void RefreshGrowthChart();
	
	/**
	 * @brief This one should never be used to check this condition
	 */
	bool HasChanged() { return false; }

	void SetGraphPoints(uint initGen, uint finalGen, uint period);

protected:
	wxBitmap imgData;
	string growthChartFilename;
	bool refreshGrowthChart;
	long initGen;
	long finalGen;
	long period;


};



}

}

#endif
