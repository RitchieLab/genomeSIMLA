/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpageselectloci.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 11:10:40 CDT
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

#include "wxwizpageselectloci.h"
#include "wxdlgselectlocusrange.h"
#include <sstream>
#include <wx/colour.h>
#include "diseasemodeldetails.h"
#include "utility/exception.h"
#include "wxdialogtaskalert.h"
////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

/*!
 * wxWizPageSelectLoci type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageSelectLoci, wxWizardPageDataSim )


/*!
 * wxWizPageSelectLoci event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageSelectLoci, wxWizardPageDataSim )

////@begin wxWizPageSelectLoci event table entries
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageSelectLoci::OnPageSelectLociPageChanging )
    EVT_SIZE( wxWizPageSelectLoci::OnSize )

    EVT_GRID_CELL_LEFT_DCLICK( wxWizPageSelectLoci::OnLeftDClick )

    EVT_SPINCTRL( ID_SPINCTRL_LOCUSCOUNT, wxWizPageSelectLoci::OnSpinctrlUpdated )

////@end wxWizPageSelectLoci event table entries
END_EVENT_TABLE()


/*!
 * wxWizPageSelectLoci constructors
 */

wxWizPageSelectLoci::wxWizPageSelectLoci() : modelDetails(NULL), directionOfTransition(true)
{
    Init();
}

wxWizPageSelectLoci::wxWizPageSelectLoci( wxWizard* parent ): modelDetails(NULL), directionOfTransition(true)
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageSelectLoci creator
 */

bool wxWizPageSelectLoci::Create( wxWizard* parent )
{
////@begin wxWizPageSelectLoci creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageDataSim::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageSelectLoci creation
    return true;
}


/*!
 * wxWizPageSelectLoci destructor
 */

wxWizPageSelectLoci::~wxWizPageSelectLoci()
{
////@begin wxWizPageSelectLoci destruction
////@end wxWizPageSelectLoci destruction
	if (modelDetails)
		delete modelDetails;

}


/*!
 * Member initialisation
 */

void wxWizPageSelectLoci::Init()
{
////@begin wxWizPageSelectLoci member initialisation
    grdLocusSelection = NULL;
    spnLocusCount = NULL;
////@end wxWizPageSelectLoci member initialisation
}


/*!
 * Control creation for wxWizPageSelectLoci
 */

void wxWizPageSelectLoci::CreateControls()
{    
////@begin wxWizPageSelectLoci content construction
    wxWizPageSelectLoci* itemWizardPageDataSim1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageDataSim1->SetSizer(itemBoxSizer2);

    wxStaticText* itemStaticText3 = new wxStaticText( itemWizardPageDataSim1, wxID_STATIC, _("Select Loci associated with the disease model selected"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer2->Add(itemStaticText3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    grdLocusSelection = new wxGrid( itemWizardPageDataSim1, ID_GRD_LOCI, wxDefaultPosition, wxSize(200, 150), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdLocusSelection->SetDefaultColSize(50);
    grdLocusSelection->SetDefaultRowSize(25);
    grdLocusSelection->SetColLabelSize(25);
    grdLocusSelection->SetRowLabelSize(50);
    itemBoxSizer2->Add(grdLocusSelection, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer5, 0, wxALIGN_RIGHT|wxALL, 5);

    wxStaticText* itemStaticText6 = new wxStaticText( itemWizardPageDataSim1, wxID_STATIC, _("Locus Count (SIMPEN Only):"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer5->Add(itemStaticText6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    spnLocusCount = new wxSpinCtrl( itemWizardPageDataSim1, ID_SPINCTRL_LOCUSCOUNT, _T("0"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 100, 2 );
    itemBoxSizer5->Add(spnLocusCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxWizPageSelectLoci content construction
	grdLocusSelection->ClearGrid();
	grdLocusSelection->CreateGrid(0, 4);

	grdLocusSelection->SetRowLabelSize(20);

}

void wxWizPageSelectLoci::RefreshSize() {
	int width, height;

	if (grdLocusSelection) {	
		grdLocusSelection->GetSize(&width, &height);
	
		grdLocusSelection->SetColSize(0, (int)(width * 0.24));
		grdLocusSelection->SetColSize(1, (int)(width * 0.24));
		grdLocusSelection->SetColSize(2, (int)(width * 0.24));
		grdLocusSelection->SetColSize(3, (int)(width * 0.24));
	}
}
/*!
 * Should we show tooltips?
 */

bool wxWizPageSelectLoci::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageSelectLoci::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageSelectLoci bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageSelectLoci bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageSelectLoci::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageSelectLoci icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageSelectLoci icon retrieval
}

bool wxWizPageSelectLoci::Initialize(vector<FileBasedChromosome *> *chrom) {
	uint chromCount = chrom->size();
	vector<FileBasedChromosome *>::iterator itr = chrom->begin();
	vector<FileBasedChromosome *>::iterator end = chrom->end();
	labels = new wxString[chromCount];

	int i=0;
	while (itr != end) {
		FileBasedChromosome *cur = *itr++;
		labels[i++] = wxT(cur->label.c_str());
		chromosomes[cur->label] = cur;
	}

	//gridRegions->SetColLabelSize(40);
//	grdLocusSelection->InsertColumn(0, wxT("Locus Label"));
	grdLocusSelection->SetColLabelValue(0, wxT("Chromosome"));
	grdLocusSelection->SetColLabelValue(1, wxT("SNP ID"));
	grdLocusSelection->SetColLabelValue(2, wxT("Freq (al1)"));
	grdLocusSelection->SetColLabelValue(3, wxT("Disease Freq (al1)"));

	wxGridCellAttr *attrCombo  = new wxGridCellAttr;
	attrCombo->SetEditor(new wxGridCellChoiceEditor(chromCount, labels));
	grdLocusSelection->SetColAttr(0, attrCombo);

	RefreshSize();
	return i>0;
}


void wxWizPageSelectLoci::SetLocus(int idx, Locus *locus) {
	assert(idx < grdLocusSelection->GetNumberRows());

	stringstream ss;
	ss<<locus->Freq1();

	
	wxColour color(255, 255, 255);
	float goal, min, max;
	if (modelDetails && modelDetails->GetLocusRange(idx, min, max, goal)) {
		if (locus->Freq1() < min || locus->Freq1() > max)
			color.Set(255, 32, 32);
	}
	
	grdLocusSelection->SetCellValue(idx, 0, labels[locus->GetChromID()]);
	grdLocusSelection->SetCellBackgroundColour(idx, 0, color);
	grdLocusSelection->SetCellValue(idx, 1, wxT(locus->GetLabel().c_str()));
	grdLocusSelection->SetCellBackgroundColour(idx, 1, color);
	grdLocusSelection->SetCellValue(idx, 2, wxT(ss.str().c_str()));
	grdLocusSelection->SetCellBackgroundColour(idx, 2, color);
	loci[idx]=locus;
}

bool wxWizPageSelectLoci::TransferDataFromWindow() {
	string error;
	bool doContinue = modelDetails->ValidateCfg(error);
	if (!doContinue)
		wxMessageBox(error.c_str());
	return doContinue;
;
}


void wxWizPageSelectLoci::SetLocusCount(int locusCount) {
	cout<<"+++++++++++++++ SetLocusCount ("<<locusCount<<")\n";

	//SIMpen models accept variable number of loci
	if (locusCount < 0)
		locusCount = 2;

	spnLocusCount->SetValue(locusCount);

	grdLocusSelection->DeleteRows(0, grdLocusSelection->GetNumberRows(), true);
	loci.clear();
	
	FileBasedChromosome *chr = chromosomes.begin()->second;

	char label[]="A";
	for (int i=0; i<locusCount; i++){ 
		loci.push_back(NULL);
		grdLocusSelection->AppendRows();
		grdLocusSelection->SetRowLabelValue(i, wxT(label));
		SetLocus(i, chr->GetLocus(i));
		if (modelDetails)
			grdLocusSelection->SetCellValue(i, 3, wxT(modelDetails->GetLocusRange(i).c_str()));
		grdLocusSelection->SetReadOnly(i, 1);
		label[0]++;
	}
	grdLocusSelection->Refresh();
}

void wxWizPageSelectLoci::SetModel(ModelList::ModelFileItem &item) {
	diseaseModel = item;
	if (modelDetails)
		delete modelDetails;
	modelDetails = DiseaseModelDetails::OpenDiseaseModel(diseaseModel.modelType.c_str(), diseaseModel.filename.c_str());

	modelDetails->Load();

	string modelError;
	if (!modelDetails->ValidateCfg(modelError))	
		throw Utility::Exception::General(modelError.c_str());

	SetLocusCount(item.locusCount);
	spnLocusCount->Enable(item.locusCount<0);
}



/*!
 * wxEVT_GRID_CELL_LEFT_DCLICK event handler for ID_GRD_LOCI
 */

void wxWizPageSelectLoci::OnLeftDClick( wxGridEvent& event )
{
	int row = event.GetRow();
	int col = event.GetCol();
	wxString selection = grdLocusSelection->GetCellValue(row, 0);

	if (row >= 0 && col == 1) {
		FileBasedChromosome *chr = chromosomes[selection.c_str()];
		Locus *l = chr->GetLocus(grdLocusSelection->GetCellValue(row, 1).c_str());
		
		wxDlgSelectLocusRange dlg(NULL, wxID_ANY);
		dlg.Initialize(chr, l);
		
		if (dlg.ShowModal() == wxID_OK) {
			l = dlg.GetSelection();				
			if (l) {
				SetLocus(row, l);
			}
		}

	}	
}


/*!
 * wxEVT_SIZE event handler for ID_PAGE_SELECT_LOCI
 */

void wxWizPageSelectLoci::OnSize( wxSizeEvent& event )
{
	RefreshSize();
	event.Skip();
}

struct thArgs {
	stringstream *ss;
	vector<Locus *> *loci;
	AppController *controller;
	wxDialogTaskAlert *alert;

	bool isCompleted;
	bool isCancelled;

	thArgs(stringstream *ss, vector<Locus*>*loci, AppController*controller, wxDialogTaskAlert *alert) : ss(ss), loci(loci), controller(controller), alert(alert), isCompleted(false), isCancelled(false) {
cerr<<"Locus 0: "<<(*loci)[0]->GetLabel()<<"\n";
	}

	bool DoContinue() { return !isCancelled && isCompleted; }
};


void *wxWizPageSelectLoci::SummarizeDiseaseModel(void *argv) {
	thArgs *args = (thArgs*)argv;
	args->isCompleted = false;
	try {
cerr<<"LOCUS 0: "<<(*args->loci)[0]->GetLabel()<<"\n";
		args->controller->SummarizeDiseaseModel(*args->ss, *args->loci);
		args->isCompleted = true;
	} catch (Utility::Exception::General& e) {
		wxMessageBox(e.GetErrorMessage().c_str());
//		args->isCompleted = true;
		args->isCancelled = true;
	}

}

void wxWizPageSelectLoci::PrepSummary(wxWizardPageDataSummary* summary) {
	((wxWizardPageDataSim*)GetPrev())->PrepSummary(summary);
	((wxWizardPageDataSim*)GetNext())->SetDiseaseLoci(loci);
	((wxWizardPageDataSim*)GetNext())->SetDiseaseModel(diseaseModel);
	stringstream ss;

	int locusCount = loci.size();
	summary->AddHeader("Model Loci");
	ss<<"<TABLE border='1'><TR bgcolor=\"#dddddd\"><TH>Model Locus</TH><TH>Chromosome</TH><TH>SNP Index</TH><TH>SNP Label</TH><TH>Freq Al 1</TH></TR>";
	char label[] = "A";
	for (int i=0; i<locusCount; i++) {
		Locus *locus = loci[i];
		ss<<"<TR><TD>"<<label<<"</TD><TD>"<<labels[locus->GetChromID()]<<"</TD>"
			<<"<TD>"<<locus->GetID() + 1<<"</TD></TD><TD>"<<locus->GetLabel()<<"</TD>"
			<<"<TD>"<<locus->Freq1()<<"</TD></TR>";
		label[0]++;
	}
	ss<<"</TABLE>";
cerr<<"LOCUS ID(0) "<<loci[0]->GetLabel()<<"\n";

	string modelConfiguration;
	try {
		modelConfiguration = modelDetails->GetConfigurationDetails(loci).c_str();
	}
	catch (Utility::Exception::General& e) {
		wxMessageBox(_(e.GetErrorMessage().c_str()));
		return;
	}
	cout<<"Commiting Model: "<<modelConfiguration<<"\n";
	appController->parameters.DefineDiseaseModel(modelConfiguration.c_str());

		streambuf *buff = cout.rdbuf();
		stringstream outbuffer;
		cout.rdbuf(outbuffer.rdbuf());
		wxDialogTaskAlert *alert =  new wxDialogTaskAlert(this, true);
cerr<<"LOCUS ID(0) "<<loci[0]->GetLabel()<<"\n";
	try {
		thArgs thargs(&ss, &loci, appController, alert);
		pthread_t summaryThread;
		pthread_create(&summaryThread, NULL, &SummarizeDiseaseModel, &thargs);
		
		
		alert->SetTitle("Please Wait");
		alert->SetMessage("Setting up disease model");
	
		SetCursor(wxCursor(wxCURSOR_WAIT));
		alert->Show();

		do {
			wxMilliSleep(50);
			string buff = outbuffer.str();
			outbuffer.str("");
			
			alert->WriteLog(_T(buff.c_str()));
			alert->Update();
			wxYield();
		} while (!thargs.isCompleted && !thargs.isCancelled);
		
		alert->ShowCompleted();

		pthread_join(summaryThread, NULL);
		alert->ShowModal();

		SetCursor(wxCursor(wxCURSOR_ARROW));
		alert->Destroy();
	}
	catch (Utility::Exception::General& e) {
		wxMessageBox(_(e.GetErrorMessage().c_str()));
		cout.rdbuf(buff);
		alert->Destroy();
		return;
	}

		cout.rdbuf(buff);
		alert->Destroy();
	cout<<ss.str()<<"\n";

	summary->AddNote(ss.str().c_str());

}


void wxWizPageSelectLoci::PrepConfig(ostream& config) {
	config<<modelDetails->GetConfigurationDetails(loci);
}


/**
 * @brief This is called prior to saving the configuration
 */
void wxWizPageSelectLoci::Commit() { 

}

/**
 * @brief This is called just after a configuration has been loaded
 */
void wxWizPageSelectLoci::RefreshSettings() { 
	PenetranceModel *model = appController->parameters.GetDiseaseModel();
	string modelCfg = "";
	if (model)
		modelCfg = model->GetModelConfiguration();

	cout<<"Selected Model: "<<modelCfg<<"\n";
	//I need to figure out how to extract the loci that were associated with a previous run
}

/**
 * @brief Compares current state with that of the saved version at least against what is in memory
 * @return True indicates that something has changed 
 */
bool wxWizPageSelectLoci::HasChanged() { 
	return true;
} 


/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_PAGE_SELECT_LOCI
 */

void wxWizPageSelectLoci::OnPageSelectLociPageChanging( wxWizardEvent& event )
{
	directionOfTransition = event.GetDirection();

	if (directionOfTransition) {
		//If there is something that should cause us to refuse to move to the next page, 
		//we can veto here by returning false;
		Commit();
	
		string modelConfiguration;
		try {
			modelConfiguration = modelDetails->GetConfigurationDetails(loci).c_str();
		}
		catch (Utility::Exception::General& e) {
			wxMessageBox(_(e.GetErrorMessage().c_str()));
			event.Veto();
		}
	
		bool validSelections = true;
		float goal = 0.0, min = 0.0, max = 0.0;
	
		for (int i=0; i<loci.size(); i++) {
			if (modelDetails && modelDetails->GetLocusRange(i, min, max, goal)) {
				float curFreq1 = loci[i]->Freq1();
				validSelections = validSelections && (loci[i]->Freq1() >= min && loci[i]->Freq1() <= max);
			}
		}
	
		if (!validSelections)  {
			wxMessageDialog dlg(NULL, wxT("One or more of the selected loci fail to meet the necessary allele frequency range."), wxT("Invalid Locus Selection"), wxICON_ERROR|wxOK);
			dlg.ShowModal();
			event.Veto();
		}
		
			
		
	}
}

/*!
 * wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_SPINCTRL
 */

void wxWizPageSelectLoci::OnSpinctrlUpdated( wxSpinEvent& event )
{
	SetLocusCount(spnLocusCount->GetValue());
}


}

}



