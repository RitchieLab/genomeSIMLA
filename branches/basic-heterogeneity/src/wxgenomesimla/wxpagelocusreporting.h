/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagelocusreporting.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Tue 11 Dec 2007 09:21:06 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPAGELOCUSREPORTING_H_
#define _WXPAGELOCUSREPORTING_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/listctrl.h"
#include "wx/toolbar.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxListCtrl;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_PAGE_LOCUS_REPORTING 10043
#define ID_LST_SELECTORS 10045
#define ID_CMD_ADD 10000
#define ID_CMD_REMOVE 10001
#define SYMBOL_WXPAGELOCUSREPORTING_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGELOCUSREPORTING_TITLE _("Locus Reporting")
#define SYMBOL_WXPAGELOCUSREPORTING_IDNAME ID_PAGE_LOCUS_REPORTING
#define SYMBOL_WXPAGELOCUSREPORTING_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGELOCUSREPORTING_POSITION wxDefaultPosition
////@end control identifiers


#include "appinterface.h"
#include "locusmanager.h"

namespace GenomeSIM {

namespace GUI {
/*!
 * wxPageLocusReporting class declaration
 */

class wxPageLocusReporting: public wxPanel, public AppInterface
{    
    DECLARE_DYNAMIC_CLASS( wxPageLocusReporting )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPageLocusReporting();
    wxPageLocusReporting( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGELOCUSREPORTING_IDNAME, const wxPoint& pos = SYMBOL_WXPAGELOCUSREPORTING_POSITION, const wxSize& size = SYMBOL_WXPAGELOCUSREPORTING_SIZE, long style = SYMBOL_WXPAGELOCUSREPORTING_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGELOCUSREPORTING_IDNAME, const wxPoint& pos = SYMBOL_WXPAGELOCUSREPORTING_POSITION, const wxSize& size = SYMBOL_WXPAGELOCUSREPORTING_SIZE, long style = SYMBOL_WXPAGELOCUSREPORTING_STYLE );

    /// Destructor
    ~wxPageLocusReporting();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxPageLocusReporting event handler declarations

    /// wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_SELECTORS
    void OnLstSelectorsSelected( wxListEvent& event );

    /// wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_SELECTORS
    void OnLstSelectorsDeselected( wxListEvent& event );

    /// wxEVT_SIZE event handler for ID_LST_SELECTORS
    void OnSize( wxSizeEvent& event );

    /// wxEVT_LEFT_DCLICK event handler for ID_LST_SELECTORS
    void OnLeftDClick( wxMouseEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD
    void OnCmdAddClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE
    void OnCmdRemoveClick( wxCommandEvent& event );

////@end wxPageLocusReporting event handler declarations

////@begin wxPageLocusReporting member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPageLocusReporting member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPageLocusReporting member variables
    wxListCtrl* lstLocusSelectors;
////@end wxPageLocusReporting member variables
	void OnLeftDouble( wxListEvent &event);

	void InitAppController(AppController *appController);
	void InitList();
	void SetSelector(long idx, LocusSelection &sel);
	void AddSelector(LocusSelection &sel);
	void RefreshSize();

	void SetChromosomes(vector<FileBasedChromosome *> *ch) { chromosomes=ch; }

	/**
	 * @brief This is called prior to saving the configuration
	 */
	void Commit();

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings();

	bool HasChanged();
	void ClearSelections();
	bool VerifyForRun() { return true; }

protected:
	vector<LocusSelection> selections;				
	vector<FileBasedChromosome *> *chromosomes;		///<The underlying data to be used for selection
	long selectedRow;								///<This is a cache to the currently selected row
};


}

}

#endif
    // _WXPAGELOCUSREPORTING_H_
