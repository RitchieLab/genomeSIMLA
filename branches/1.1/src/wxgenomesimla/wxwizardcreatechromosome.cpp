/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizardcreatechromosome.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 04 Apr 2008 11:20:39 CDT
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
#include "wxwizpagebbchrom.h"
#include "wxpanelchromosomerep.h"
////@end includes

#include "wxwizardcreatechromosome.h"
#include "simulation/chrompool.h"

////@begin XPM images
#include "img/DnaSketch.xpm"
////@end XPM images

namespace GenomeSIM {

namespace GUI {
/*!
 * wxWizardCreateChromosome type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizardCreateChromosome, wxWizard )


/*!
 * wxWizardCreateChromosome event table definition
 */

BEGIN_EVENT_TABLE( wxWizardCreateChromosome, wxWizard )

////@begin wxWizardCreateChromosome event table entries
    EVT_WIZARD_FINISHED( ID_WXWIZARDCREATECHROMOSOME, wxWizardCreateChromosome::OnWxwizardcreatechromosomeFinished )

////@end wxWizardCreateChromosome event table entries

END_EVENT_TABLE()


/*!
 * wxWizardCreateChromosome constructors
 */

wxWizardCreateChromosome::wxWizardCreateChromosome() 
{
    Init();
}

wxWizardCreateChromosome::wxWizardCreateChromosome( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
    Init();
    Create(parent, id, pos);
}


/*!
 * wxWizardCreateChromosome creator
 */

bool wxWizardCreateChromosome::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos )
{
////@begin wxWizardCreateChromosome creation
    SetExtraStyle(wxWIZARD_EX_HELPBUTTON);
    wxBitmap wizardBitmap(GetBitmapResource(wxT("img/DnaSketch.xpm")));
    wxWizard::Create( parent, id, _("Create Chromosome"), wizardBitmap, pos, wxDEFAULT_DIALOG_STYLE|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX );

    CreateControls();
////@end wxWizardCreateChromosome creation
    return true;
}


/*!
 * wxWizardCreateChromosome destructor
 */

wxWizardCreateChromosome::~wxWizardCreateChromosome()
{
////@begin wxWizardCreateChromosome destruction
////@end wxWizardCreateChromosome destruction
}


/*!
 * Member initialisation
 */

void wxWizardCreateChromosome::Init()
{
////@begin wxWizardCreateChromosome member initialisation
    pageConfigureChrom = NULL;
////@end wxWizardCreateChromosome member initialisation
}


/*!
 * Control creation for wxWizardCreateChromosome
 */

void wxWizardCreateChromosome::CreateControls()
{    
////@begin wxWizardCreateChromosome content construction
    wxWizardCreateChromosome* itemWizard1 = this;

    pageConfigureChrom = new WizardPage( itemWizard1 );
    itemWizard1->GetPageAreaSizer()->Add(pageConfigureChrom);

    wxWizardPageSimple* lastPage = NULL;
    if (lastPage)
        wxWizardPageSimple::Chain(lastPage, pageConfigureChrom);
    lastPage = pageConfigureChrom;
////@end wxWizardCreateChromosome content construction
}


/*!
 * Runs the wizard.
 */

bool wxWizardCreateChromosome::Run()
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

bool wxWizardCreateChromosome::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizardCreateChromosome::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizardCreateChromosome bitmap retrieval
    wxUnusedVar(name);
    if (name == _T("img/DnaSketch.xpm"))
    {
        wxBitmap bitmap( DnaSketch_xpm);
        return bitmap;
    }
    return wxNullBitmap;
////@end wxWizardCreateChromosome bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizardCreateChromosome::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizardCreateChromosome icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizardCreateChromosome icon retrieval
}

void wxWizardCreateChromosome::Initialize(int chromID) {
	if (pageConfigureChrom) {
		pageConfigureChrom->panelConfigChrom->Initialize(chromID);
	}
}
/*!
 * WizardPage type definition
 */

IMPLEMENT_DYNAMIC_CLASS( WizardPage, wxWizardPageSimple )


/*!
 * WizardPage event table definition
 */

BEGIN_EVENT_TABLE( WizardPage, wxWizardPageSimple )

////@begin WizardPage event table entries
    EVT_TEXT( ID_TXT_MIN_SNPCOUNT, WizardPage::OnTxtMinBlockRecombUpdated )

    EVT_TEXT( ID_TXT_MAX_SNPCOUNT, WizardPage::OnTxtMaxBlockRecombUpdated )

    EVT_TEXT( ID_TXT_MIN_BLRECOMB, WizardPage::OnTxtMinBlockRecombUpdated )

    EVT_TEXT( ID_TXT_MAX_BLRECOMB, WizardPage::OnTxtMaxBlockRecombUpdated )

    EVT_TEXT( ID_TXT_MIN_SNPRECOMB, WizardPage::OnTxtDefMinBlInteriorUpdated )

    EVT_TEXT( ID_TXT_MAX_SNPRECOMB, WizardPage::OnTxtDefMaxBlInteriorUpdated )

    EVT_BUTTON( ID_CMD_REGENERATE, WizardPage::OnCmdRegenerateClick )

////@end WizardPage event table entries

END_EVENT_TABLE()


/*!
 * WizardPage constructors
 */

WizardPage::WizardPage() : chromID(0)
{
    Init();
}

WizardPage::WizardPage( wxWizard* parent ) : chromID(0)
{
    Init();
    Create( parent );
}


/*!
 * WizardPage creator
 */

bool WizardPage::Create( wxWizard* parent )
{
////@begin WizardPage creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageSimple::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end WizardPage creation
    return true;
}


/*!
 * WizardPage destructor
 */

WizardPage::~WizardPage()
{
////@begin WizardPage destruction
////@end WizardPage destruction
}


/*!
 * Member initialisation
 */

void WizardPage::Init()
{
////@begin WizardPage member initialisation
    txtMinSnpCount = NULL;
    txtMaxSnpCount = NULL;
    txtMinBlckRecomb = NULL;
    txtMaxBlckRecomb = NULL;
    txtMinSnpRecomb = NULL;
    txtMaxSnpRecomb = NULL;
    panelConfigChrom = NULL;
    panelChromosome = NULL;
    txtMinFreq = NULL;
    txtMaxFreq = NULL;
////@end WizardPage member initialisation
}


/*!
 * Control creation for WizardPage
 */

void WizardPage::CreateControls()
{    
////@begin WizardPage content construction
    WizardPage* itemWizardPageSimple2 = this;

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageSimple2->SetSizer(itemBoxSizer3);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer3->Add(itemBoxSizer4, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer5Static = new wxStaticBox(itemWizardPageSimple2, wxID_ANY, _("Default Block Size"));
    wxStaticBoxSizer* itemStaticBoxSizer5 = new wxStaticBoxSizer(itemStaticBoxSizer5Static, wxVERTICAL);
    itemBoxSizer4->Add(itemStaticBoxSizer5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer6 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer5->Add(itemBoxSizer6, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText7 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Min. SNP Count:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer6->Add(itemStaticText7, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMinSnpCount = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MIN_SNPCOUNT, _("2"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    txtMinSnpCount->SetMaxLength(25);
    if (WizardPage::ShowToolTips())
        txtMinSnpCount->SetToolTip(_("The minimum distance between a default block and the previous block"));
    itemBoxSizer6->Add(txtMinSnpCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer9 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer5->Add(itemBoxSizer9, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText10 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Max. SNP Count:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer9->Add(itemStaticText10, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMaxSnpCount = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MAX_SNPCOUNT, _("5"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    if (WizardPage::ShowToolTips())
        txtMaxSnpCount->SetToolTip(_("The maximum distance between a default block and the previous block"));
    itemBoxSizer9->Add(txtMaxSnpCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer12Static = new wxStaticBox(itemWizardPageSimple2, wxID_ANY, _("Default Block Separation"));
    wxStaticBoxSizer* itemStaticBoxSizer12 = new wxStaticBoxSizer(itemStaticBoxSizer12Static, wxVERTICAL);
    itemBoxSizer4->Add(itemStaticBoxSizer12, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer12->Add(itemBoxSizer13, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText14 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Min. Map Distance:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemStaticText14, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMinBlckRecomb = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MIN_BLRECOMB, _("0.00001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    if (WizardPage::ShowToolTips())
        txtMinBlckRecomb->SetToolTip(_("The minimum distance between a default block and the previous block"));
    itemBoxSizer13->Add(txtMinBlckRecomb, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer16 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer12->Add(itemBoxSizer16, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText17 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Max. Mp Distance:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer16->Add(itemStaticText17, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMaxBlckRecomb = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MAX_BLRECOMB, _("0.0001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    if (WizardPage::ShowToolTips())
        txtMaxBlckRecomb->SetToolTip(_("The maximum distance between a default block and the previous block"));
    itemBoxSizer16->Add(txtMaxBlckRecomb, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer19Static = new wxStaticBox(itemWizardPageSimple2, wxID_ANY, _("Default Block Interior"));
    wxStaticBoxSizer* itemStaticBoxSizer19 = new wxStaticBoxSizer(itemStaticBoxSizer19Static, wxVERTICAL);
    itemBoxSizer4->Add(itemStaticBoxSizer19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 0);

    wxBoxSizer* itemBoxSizer20 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer19->Add(itemBoxSizer20, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText21 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Min. Distance:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer20->Add(itemStaticText21, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMinSnpRecomb = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MIN_SNPRECOMB, _("0.00001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    if (WizardPage::ShowToolTips())
        txtMinSnpRecomb->SetToolTip(_("The minimum distance between SNPs in the local block's interior"));
    itemBoxSizer20->Add(txtMinSnpRecomb, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer23 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer19->Add(itemBoxSizer23, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 0);

    wxStaticText* itemStaticText24 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Max. Map Distance:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer23->Add(itemStaticText24, 0, wxALIGN_CENTER_VERTICAL|wxALL|wxADJUST_MINSIZE, 5);

    txtMaxSnpRecomb = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MAX_SNPRECOMB, _("0.0001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    if (WizardPage::ShowToolTips())
        txtMaxSnpRecomb->SetToolTip(_("The maximum distance between SNPs in the local block's interior"));
    itemBoxSizer23->Add(txtMaxSnpRecomb, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    panelConfigChrom = new wxWizPageBBChrom( itemWizardPageSimple2, ID_PANEL_CHROM_CFG, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
    itemBoxSizer3->Add(panelConfigChrom, 3, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer27 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer3->Add(itemBoxSizer27, 0, wxGROW|wxALL, 5);

    panelChromosome = new wxPanelChromosomeRep( itemWizardPageSimple2, ID_CHROMOSOME_VISUAL, wxDefaultPosition, wxSize(-1, 32), wxRAISED_BORDER|wxTAB_TRAVERSAL );
    if (ShowToolTips())
        panelChromosome->SetToolTip(_("Click to Regenerate the Chromosome"));
    itemBoxSizer27->Add(panelChromosome, 1, wxGROW|wxALL, 5);

    wxButton* itemButton29 = new wxButton( itemWizardPageSimple2, ID_CMD_REGENERATE, _("Regenerate"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer27->Add(itemButton29, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer30Static = new wxStaticBox(itemWizardPageSimple2, wxID_ANY, _("Minor Allele Frequency boundaries for the Loci to be Produced"));
    wxStaticBoxSizer* itemStaticBoxSizer30 = new wxStaticBoxSizer(itemStaticBoxSizer30Static, wxHORIZONTAL);
    itemBoxSizer3->Add(itemStaticBoxSizer30, 0, wxGROW|wxALL, 5);

    itemStaticBoxSizer30->Add(5, 5, 0, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer32 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer30->Add(itemBoxSizer32, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText33 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Min. Allele Freq."), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer32->Add(itemStaticText33, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMinFreq = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MIN_FREQ, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer32->Add(txtMinFreq, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemStaticBoxSizer30->Add(5, 5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer36 = new wxBoxSizer(wxHORIZONTAL);
    itemStaticBoxSizer30->Add(itemBoxSizer36, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxStaticText* itemStaticText37 = new wxStaticText( itemWizardPageSimple2, wxID_STATIC, _("Max. Allele Freq."), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer36->Add(itemStaticText37, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtMaxFreq = new wxTextCtrl( itemWizardPageSimple2, ID_TXT_MAX_FREQ, _T(""), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer36->Add(txtMaxFreq, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemStaticBoxSizer30->Add(5, 5, 1, wxGROW|wxALL, 5);

////@end WizardPage content construction

	AppConfig *config = AppConfig::GetPrimaryInstance();
	ChromPool::BlockDefinition *defBlock = config->GetDefaultBlock();
	UpdateTextField(txtMaxBlckRecomb, defBlock->maxBlckMap);
	UpdateTextField(txtMinBlckRecomb, defBlock->minBlckMap);
	UpdateTextField(txtMaxSnpRecomb, defBlock->maxSnpMap);
	UpdateTextField(txtMinSnpRecomb, defBlock->minSnpMap);
	UpdateTextField(txtMaxSnpCount, (int)defBlock->maxSnpCount);
	UpdateTextField(txtMinSnpCount, (int)defBlock->minSnpCount);
	UpdateTextField(txtMinFreq, ChromPool::defFre1);
	UpdateTextField(txtMaxFreq, ChromPool::defFre2);
}



bool WizardPage::TransferDataFromWindow() {
	if (VerifyFreq()) {
	
		int minSnpCount = ExtractInteger(txtMinSnpCount),
			maxSnpCount = ExtractInteger(txtMaxSnpCount);
	
		double minBlck = ExtractDouble(txtMinBlckRecomb),
			maxBlck = ExtractDouble(txtMaxBlckRecomb),
			minSnp = ExtractDouble(txtMinSnpRecomb),
			maxSnp = ExtractDouble(txtMaxSnpRecomb);
	
		ChromPool::BlockDefinition defBlock(minSnpCount, maxSnpCount, minBlck, maxBlck, minSnp, maxSnp, 0, 0);
		AppConfig *config = AppConfig::GetPrimaryInstance();
		config->SetDefaultBlock(defBlock);
	
		if (panelChromosome->GetLocusCount() == 0) {
			try {
				Regenerate();
			}
			catch (Utility::Exception::General &e){
				wxMessageBox(e.GetErrorMessage().c_str());
				return false;
			}
		}

		wxFileDialog fileSelect(this, _("Select a Filename for the new Locus File"), _("Locus File"), _("locus-report.loc"), _("Locus Report (*.loc)|*.loc"), wxSAVE|wxOVERWRITE_PROMPT);
		if (fileSelect.ShowModal() == wxID_OK) {
			wxFileName newFilename = fileSelect.GetPath();
			if (newFilename.GetExt().Len() == 0)
				newFilename.SetExt("loc");
			wxString path = newFilename.GetFullPath();
			WriteLocusFile(path.c_str());
		}
		else {
			cout<<"Trying to Veto a finish\n";
			return false;
		}
		return true;
	}
	return false;
}
/*!
 * Should we show tooltips?
 */

bool WizardPage::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap WizardPage::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin WizardPage bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end WizardPage bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon WizardPage::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin WizardPage icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end WizardPage icon retrieval
}

/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL2
 */

void WizardPage::OnTxtMinBlockRecombUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL2 in WizardPage.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL2 in WizardPage. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL3
 */

void WizardPage::OnTxtMaxBlockRecombUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL3 in WizardPage.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL3 in WizardPage. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL4
 */

void WizardPage::OnTxtDefMinBlInteriorUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL4 in WizardPage.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL4 in WizardPage. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL5
 */

void WizardPage::OnTxtDefMaxBlInteriorUpdated( wxCommandEvent& event )
{
////@begin wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL5 in WizardPage.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_TEXT_UPDATED event handler for ID_TEXTCTRL5 in WizardPage. 
}



/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_BUTTON1
 */

void WizardPage::OnCmdRegenerateClick( wxCommandEvent& event )
{
	try {
		Regenerate();	
	}
	catch (Utility::Exception::General& e){ 
		wxMessageBox(e.GetErrorMessage().c_str());
	}
}


bool WizardPage::VerifyFreq() {
	double min = ExtractDouble(txtMinFreq);
	double max = ExtractDouble(txtMaxFreq);

	bool isValid = min <= max;

	

	if (isValid) {
		ChromPool::defFre1 = min;
		ChromPool::defFre2 = max;
	}
	else
		wxMessageBox(_("Minor allele frequency boundaries must be less than or equal to the maximum."), _("Invalid minor allele frequency boundaries"), wxICON_WARNING|wxOK);

	if (isValid) {
		isValid = max <= 0.5;
		
		if (!isValid) 
			wxMessageBox(_("Minor allele freqency boundaries must be less than or equal to 0.5."), _("Invalid minor allele frequency boundaries"), wxICON_WARNING|wxOK);
	}
	return isValid;
}

	
void WizardPage::Regenerate() {
	if (VerifyFreq()) {
		int minSnpCount = ExtractInteger(txtMinSnpCount),
			maxSnpCount = ExtractInteger(txtMaxSnpCount);
		double minBlck = ExtractDouble(txtMinBlckRecomb),
			maxBlck = ExtractDouble(txtMaxBlckRecomb),
			minSnp = ExtractDouble(txtMinSnpRecomb),
			maxSnp = ExtractDouble(txtMaxSnpRecomb);
	
		ChromPool::BlockDefinition defBlock(minSnpCount, maxSnpCount, minBlck, maxBlck, minSnp, maxSnp, 0, 0);
		vector<Locus> loci = panelConfigChrom->RealizeChromosome(defBlock);
		panelChromosome->SetLoci(loci);
	}
}

void WizardPage::WriteLocusFile(const char *filename) {
	this->locusFilename = filename;
	panelConfigChrom->WriteMarkerInfo(filename);
}

/*!
 * wxEVT_WIZARD_FINISHED event handler for ID_WXWIZARDCREATECHROMOSOME
 */

void wxWizardCreateChromosome::OnWxwizardcreatechromosomeFinished( wxWizardEvent& event )	{
	event.Skip();
}

/*!
 * wxEVT_LEFT_DOWN event handler for ID_CHROMOSOME_VISUAL
 */

void wxWizardCreateChromosome::OnLeftDown( wxMouseEvent& event )
{
	cout<<"Clicky clicky!\n";
////@begin wxEVT_LEFT_DOWN event handler for ID_CHROMOSOME_VISUAL in wxWizardCreateChromosome.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_LEFT_DOWN event handler for ID_CHROMOSOME_VISUAL in wxWizardCreateChromosome. 
}


}

}





