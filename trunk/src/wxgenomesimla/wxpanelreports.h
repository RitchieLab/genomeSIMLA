/////////////////////////////////////////////////////////////////////////////
// Name:        wxpageselectgeneration.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 15 Feb 2008 11:08:28 CST
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXPANELREPORT_H_
#define _WXPANELREPORT_H_


/*!
 * Includes
 */

#include "simwizard.h"
#include "appinterface.h"
#include "treelistctrl.h"
#include "utility/executionlog.h"
#include <wx/utils.h>
/*!
 * Forward declarations
 */

/*!
 * Control identifiers
 */

#include <wx/filename.h>
namespace GenomeSIM {

namespace GUI {

////@begin control identifiers
#define ID_WXPANELREPORTING 11076
#define ID_TOOLBAR2 10139
#define ID_CMD_OPEN_REPORT 10141
#define ID_CMD_CONTINUE_SIMULATION 10142
#define ID_CMD_EXTRACT_DATASETS 10143
#define ID_CMD_DETAILED_ANALYSIS 10144
#define ID_REPORT_TREE_GEN_SELECTION 11077
#define ID_MENU_LAUNCH 10054
#define ID_MENU_DETAILED_ANALYSIS 10055
#define ID_MENU_EXTRACT_DATASETS_GEN_0 10056
#define SYMBOL_WXPAGESELECTGENERATION_STYLE wxDIALOG_MODAL|wxTAB_TRAVERSAL
#define SYMBOL_WXPAGESELECTGENERATION_TITLE _("Reporting")
#define SYMBOL_WXPAGESELECTGENERATION_IDNAME ID_WXPANELREPORTING
#define SYMBOL_WXPAGESELECTGENERATION_SIZE wxSize(400, 300)
#define SYMBOL_WXPAGESELECTGENERATION_POSITION wxDefaultPosition
////@end control identifiers
#define ID_CONTINUE_SIMULATION 203123
#define ID_DETAILED_ANALYSIS   203124
#define ID_OPEN_REPORT         203125
#define ID_EXTRACT_DATASETS    203126
#define ID_MENU_LAUNCH 10054

#ifdef WIN32
#define FILE_PREDICATE "file:///"
#else
#define FILE_PREDICATE "file://"
#endif

using namespace Utility;

class LogEntryData : public wxTreeItemData {
public:
	LogEntryData(ExecutionLog::LogEntry &entry, const char *reportFilename) : entry(entry), reportFilename(reportFilename) { }
	ExecutionLog::LogEntry entry;
	string reportFilename;
	bool Open() {
		if (reportFilename.length() == 0)
			return false;
		wxFileName filename(_(reportFilename.c_str()));
		filename.MakeAbsolute();

		wxString fn;
		fn=wxString(_(FILE_PREDICATE))+filename.GetFullPath();
		fn.Replace(" ", "%20", true);

		wxLaunchDefaultBrowser(fn);
			return true;
	}

	void GetMenu(wxMenu& menu) {
		if (reportFilename.length() == 0) {
			menu.Append(ID_CONTINUE_SIMULATION, _T("&Continue Simulation From Here"));
			if (entry.sampledReport)
				menu.Append(ID_DETAILED_ANALYSIS, _T("&Perform Detailed Analysis"));
			menu.Append(ID_EXTRACT_DATASETS, _T("&E&xtract Data Sets From This Pool"));
		}
		else {
			menu.Append(ID_OPEN_REPORT, _T("&Open Report"));
		}	
	}

	bool CanOpenReport() {		return reportFilename.length(); }
	bool CanContinue() { return !CanOpenReport(); }
	bool CanExtractDatasets() { return CanContinue(); }
	bool CanPerfromDetailed() { return CanContinue(); }

};

/*!
 * wxPanelReports class declaration
 */

class wxPanelReports: public wxPanel, public AppInterface
{    
    DECLARE_DYNAMIC_CLASS( wxPanelReports )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxPanelReports();
    wxPanelReports( wxWindow* parent, AppController *ctrl, wxWindowID id = SYMBOL_WXPAGESELECTGENERATION_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXPAGESELECTGENERATION_SIZE, long style = SYMBOL_WXPAGESELECTGENERATION_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXPAGESELECTGENERATION_IDNAME, const wxPoint& pos = SYMBOL_WXPAGESELECTGENERATION_POSITION, const wxSize& size = SYMBOL_WXPAGESELECTGENERATION_SIZE, long style = SYMBOL_WXPAGESELECTGENERATION_STYLE );

    /// Destructor
    ~wxPanelReports();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();


////@begin wxPanelReports event handler declarations

////@end wxPanelReports event handler declarations

////@begin wxPanelReports member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxPanelReports member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxPanelReports member variables
    wxToolBar* commandToolbar;
   	wxTreeCtrl* treeGenerationSelection;
////@end wxPanelReports member variables

	wxButton* cmdOpen;
	wxButton* cmdContinue;
	wxButton* cmdExtract;
	wxButton* cmdDetailed;
	wxButton *cmdLaunch;
	wxButton *cmdExtractFromZero;



	void OnTreeOpenReport( wxTreeEvent& event );
	void OnContextMenu(wxContextMenuEvent& event);
	void ShowContextMenu(const wxPoint& pos);
	void OnContinueSimulation(wxCommandEvent& event);
	void OnDetailedAnalysis(wxCommandEvent& event);
	void OpenReport(wxCommandEvent &event);
	void ExtractData(wxCommandEvent &event);

	void OnReportTreeItemActivated(wxTreeEvent& event);
	void OnReportTreeSelChanged(wxTreeEvent& event);


	void AddProjectEntries(const char *project, ExecutionLog::RunType *entries);

	bool ValidSelection();

	/**
	 * @brief Get the generation selected by the user. Returns t/f indicating that the value was changed
	 */
	bool GetSelectedGeneration(ExecutionLog::LogEntry &entry);

	bool TransferDataFromWindow();

	ExecutionLog::LogEntry *entry;
	LogEntryData *GetSelectedLogEntyData();
	
	void InitAppController(AppController *ctrl);
	
	/**
		* @brief This is called prior to saving the configuration
	 */
	void Commit();
	
	/**
		* @brief This is called just after a configuration has been loaded
	 */
	void RefreshSettings();
	
	bool HasChanged();
	bool VerifyForRun() { return true; }



protected:
	wxTreeItemId root;	

};

}

}
#endif
    // _WXPAGESELECTGENERATION_H_
