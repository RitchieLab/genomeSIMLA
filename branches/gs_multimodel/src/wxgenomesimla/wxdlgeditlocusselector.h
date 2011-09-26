/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgeditlocusselector.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 13 Dec 2007 14:03:56 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXDLGEDITLOCUSSELECTOR_H_
#define _WXDLGEDITLOCUSSELECTOR_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/grid.h"
#include "wx/toolbar.h"
////@end includes
#include <vector>
#include <map>
#include "simulation/locselection.h"
#include "locusmanager.h"
/*!
 * Forward declarations
 */

////@begin forward declarations
class wxGrid;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_EDITLOCUSSELECTOR 10046
#define ID_TXT_LABEL 10016
#define ID_TXT_DESCRIPTION 10000
#define ID_TXT_BL_TARGET 10006
#define ID_TXT_BL_MIN 10008
#define ID_TXT_BL_MAX 10009
#define ID_TXT_MAC_TARGET 10010
#define ID_TXT_MAF_MIN 10011
#define ID_TXT_MAF_MAX 10012
#define ID_GRD_REGIONS 10049
#define ID_CMD_ADD_REGION 10005
#define ID_CMD_REMOVE_REGION 10007
#define SYMBOL_EDITLOCUSSELECTOR_STYLE wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX|wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_EDITLOCUSSELECTOR_TITLE _("Edit Locus Selector")
#define SYMBOL_EDITLOCUSSELECTOR_IDNAME ID_EDITLOCUSSELECTOR
#define SYMBOL_EDITLOCUSSELECTOR_SIZE wxSize(400, 300)
#define SYMBOL_EDITLOCUSSELECTOR_POSITION wxDefaultPosition
////@end control identifiers

namespace GenomeSIM {

namespace GUI {

using namespace Simulation::Visualization;
using namespace std;



/*!
 * EditLocusSelector class declaration
 */

class EditLocusSelector: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( EditLocusSelector )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    EditLocusSelector();
    EditLocusSelector( wxWindow* parent, wxWindowID id = SYMBOL_EDITLOCUSSELECTOR_IDNAME, const wxString& caption = SYMBOL_EDITLOCUSSELECTOR_TITLE, const wxPoint& pos = SYMBOL_EDITLOCUSSELECTOR_POSITION, const wxSize& size = SYMBOL_EDITLOCUSSELECTOR_SIZE, long style = SYMBOL_EDITLOCUSSELECTOR_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_EDITLOCUSSELECTOR_IDNAME, const wxString& caption = SYMBOL_EDITLOCUSSELECTOR_TITLE, const wxPoint& pos = SYMBOL_EDITLOCUSSELECTOR_POSITION, const wxSize& size = SYMBOL_EDITLOCUSSELECTOR_SIZE, long style = SYMBOL_EDITLOCUSSELECTOR_STYLE );

    /// Destructor
    ~EditLocusSelector();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin EditLocusSelector event handler declarations

    /// wxEVT_SIZE event handler for ID_EDITLOCUSSELECTOR
    void OnSize( wxSizeEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_LABEL
    void OnTxtLabelTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_DESCRIPTION
    void OnTxtDescriptionTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_TARGET
    void OnTxtBlTargetTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MIN
    void OnTxtBlMinTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_BL_MAX
    void OnTxtBlMaxTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAC_TARGET
    void OnTxtMacTargetTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MIN
    void OnTxtMafMinTextUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAF_MAX
    void OnTxtMafMaxTextUpdated( wxCommandEvent& event );

    /// wxEVT_GRID_CELL_LEFT_DCLICK event handler for ID_GRD_REGIONS
    void OnLeftDClick( wxGridEvent& event );

    /// wxEVT_GRID_CELL_CHANGE event handler for ID_GRD_REGIONS
    void OnCellChange( wxGridEvent& event );

    /// wxEVT_GRID_EDITOR_HIDDEN event handler for ID_GRD_REGIONS
    void OnEditorHidden( wxGridEvent& event );

    /// wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_REGIONS
    void OnEditorShown( wxGridEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_ADD_REGION
    void OnCmdAddRegionClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_MENU_SELECTED event handler for ID_CMD_REMOVE_REGION
    void OnCmdRemoveRegionClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
    void OnOkClick( wxCommandEvent& event );

////@end EditLocusSelector event handler declarations

////@begin EditLocusSelector member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end EditLocusSelector member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin EditLocusSelector member variables
    wxTextCtrl* txtLabel;
    wxTextCtrl* txtDescription;
    wxTextCtrl* txtBlSizeTarget;
    wxTextCtrl* txtBlSizeMin;
    wxTextCtrl* txtBlSizeMax;
    wxTextCtrl* txtMafTarget;
    wxTextCtrl* txtMafMin;
    wxTextCtrl* txtMafMax;
    wxGrid* gridRegions;
////@end EditLocusSelector member variables
	LocusSelection locSelection;
	
	void RefreshSize();
	void InitList();
	void RenderRegion(int pos, Locus *start, Locus *stop);
	void Init(vector<FileBasedChromosome *> *chrom);
	void AddNewRange();
	void OnCombo(wxCommandEvent &event);
protected:
	map<string, FileBasedChromosome *> chromosomes;
	wxString *labels;							///<Seems to be required for the combo box
	
	int activeRow;								///<This is just used for recording which row is being edited
	void UpdateMAFSelection();	
	void UpdateBlockSelection();
	bool initializationComplete;
};






}

}

#endif
    // _WXDLGEDITLOCUSSELECTOR_H_
