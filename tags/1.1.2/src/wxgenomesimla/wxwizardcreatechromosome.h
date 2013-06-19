/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizardcreatechromosome.h
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 11:20:39 CDT
// RCS-ID:      
// Copyright:   Copyright 2007 Ritchie Lab - See COPYING for License 
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WXWIZARDCREATECHROMOSOME_H_
#define _WXWIZARDCREATECHROMOSOME_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/wizard.h"
////@end includes
#include "wxwizpagebbchrom.h"
#include "simwizard.h"
/*!
 * Forward declarations
 */

////@begin forward declarations
class WizardPage;
class wxWizPageBBChrom;
class wxPanelChromosomeRep;
////@end forward declarations

namespace GenomeSIM {

namespace GUI {
class WizardPage;
/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_WXWIZARDCREATECHROMOSOME 10128
#define ID_CONFIGURE_CHROM 10131
#define ID_TXT_MIN_SNPCOUNT 10000
#define ID_TXT_MAX_SNPCOUNT 10001
#define ID_TXT_MIN_BLRECOMB 10002
#define ID_TXT_MAX_BLRECOMB 10003
#define ID_TXT_MIN_SNPRECOMB 10004
#define ID_TXT_MAX_SNPRECOMB 10005
#define ID_CMD_REGENERATE 10135
#define ID_TXT_MIN_FREQ 10006
#define ID_TXT_MAX_FREQ 10137
#define SYMBOL_WXWIZARDCREATECHROMOSOME_IDNAME ID_WXWIZARDCREATECHROMOSOME
////@end control identifiers


/*!
 * WizardPage class declaration
 */

class WizardPage: public wxWizardPageSimple, public SimWizard
{    
    DECLARE_DYNAMIC_CLASS( WizardPage )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    WizardPage();

    WizardPage( wxWizard* parent );

    /// Creation
    bool Create( wxWizard* parent );

    /// Destructor
    ~WizardPage();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin WizardPage event handler declarations

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MIN_SNPCOUNT
    void OnTxtMinBlockRecombUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAX_SNPCOUNT
    void OnTxtMaxBlockRecombUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MIN_SNPRECOMB
    void OnTxtDefMinBlInteriorUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TXT_MAX_SNPRECOMB
    void OnTxtDefMaxBlInteriorUpdated( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_REGENERATE
    void OnCmdRegenerateClick( wxCommandEvent& event );

////@end WizardPage event handler declarations

////@begin WizardPage member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end WizardPage member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin WizardPage member variables
    wxTextCtrl* txtMinSnpCount;
    wxTextCtrl* txtMaxSnpCount;
    wxTextCtrl* txtMinBlckRecomb;
    wxTextCtrl* txtMaxBlckRecomb;
    wxTextCtrl* txtMinSnpRecomb;
    wxTextCtrl* txtMaxSnpRecomb;
    wxWizPageBBChrom* panelConfigChrom;
    wxPanelChromosomeRep* panelChromosome;
    wxTextCtrl* txtMinFreq;
    wxTextCtrl* txtMaxFreq;
////@end WizardPage member variables
	int chromID;

	bool TransferDataFromWindow();
	void WriteLocusFile(const char *filename);
	string locusFilename;

	bool VerifyFreq();

	void Regenerate();
};


/*!
 * wxWizardCreateChromosome class declaration
 */

class wxWizardCreateChromosome: public wxWizard
{    
    DECLARE_DYNAMIC_CLASS( wxWizardCreateChromosome )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    wxWizardCreateChromosome();
    wxWizardCreateChromosome( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZARDCREATECHROMOSOME_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WXWIZARDCREATECHROMOSOME_IDNAME, const wxPoint& pos = wxDefaultPosition );

    /// Destructor
    ~wxWizardCreateChromosome();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin wxWizardCreateChromosome event handler declarations

    /// wxEVT_WIZARD_FINISHED event handler for ID_WXWIZARDCREATECHROMOSOME
    void OnWxwizardcreatechromosomeFinished( wxWizardEvent& event );

    /// wxEVT_LEFT_DOWN event handler for ID_CHROMOSOME_VISUAL
    void OnLeftDown( wxMouseEvent& event );

////@end wxWizardCreateChromosome event handler declarations

////@begin wxWizardCreateChromosome member function declarations

    /// Runs the wizard
    bool Run();

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end wxWizardCreateChromosome member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin wxWizardCreateChromosome member variables
    WizardPage* pageConfigureChrom;
////@end wxWizardCreateChromosome member variables
	void Initialize(int chromID);
	string GetLocusFilename() { 
		return pageConfigureChrom->locusFilename;
	}
};


}

}

#endif
    // _WXWIZARDCREATECHROMOSOME_H_
