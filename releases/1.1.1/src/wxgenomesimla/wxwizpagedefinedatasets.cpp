/////////////////////////////////////////////////////////////////////////////
// Name:        wxwizpagedefinedatasets.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Thu 27 Mar 2008 11:10:58 CDT
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

#include "wxdlgconfigureccdatasets.h"
#include "wxdlgconfigurepeddatasets.h"
#include "wxdlgconfigurediseasemodel.h"
#include "wxwizpageselectloci.h"
#include <iostream>


////@begin includes
////@end includes

#include "wxwizpagedefinedatasets.h"

////@begin XPM images
////@end XPM images

namespace GenomeSIM {

namespace GUI {

using namespace std;

/*!
 * wxWizPageDefineDatasets type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxWizPageDefineDatasets, wxWizardPageDataSim )


/*!
 * wxWizPageDefineDatasets event table definition
 */

BEGIN_EVENT_TABLE( wxWizPageDefineDatasets, wxWizardPageDataSim )

////@begin wxWizPageDefineDatasets event table entries
    EVT_WIZARD_PAGE_CHANGING( -1, wxWizPageDefineDatasets::OnPageDefineDatasetsPageChanging )
    EVT_SIZE( wxWizPageDefineDatasets::OnSize )

    EVT_LIST_ITEM_SELECTED( ID_LST_DATASETS, wxWizPageDefineDatasets::OnLstDatasetsSelected )
    EVT_LIST_ITEM_DESELECTED( ID_LST_DATASETS, wxWizPageDefineDatasets::OnLstDatasetsDeselected )
    EVT_LIST_ITEM_ACTIVATED( ID_LST_DATASETS, wxWizPageDefineDatasets::OnLstDatasetsItemActivated )

    EVT_BUTTON( ID_CMD_ADD_CC, wxWizPageDefineDatasets::OnCmdAddCcClick )

    EVT_BUTTON( ID_CMD_ADD_PED, wxWizPageDefineDatasets::OnCmdAddPedClick )

    EVT_BUTTON( ID_DEL_DATASET, wxWizPageDefineDatasets::OnDelDatasetClick )

    EVT_CHOICE( ID_CMB_MODEL, wxWizPageDefineDatasets::OnCmbModelSelected )

    EVT_BUTTON( ID_CMD_CONFIG_MODEL, wxWizPageDefineDatasets::OnCmdConfigModelClick )

////@end wxWizPageDefineDatasets event table entries

END_EVENT_TABLE()


/*!
 * wxWizPageDefineDatasets constructors
 */

wxWizPageDefineDatasets::wxWizPageDefineDatasets() : modelLibrary(NULL)
{
    Init();
}

wxWizPageDefineDatasets::wxWizPageDefineDatasets( wxWizard* parent ) : modelLibrary(NULL)
{
    Init();
    Create( parent );
}


/*!
 * wxWizPageDefineDatasets creator
 */

bool wxWizPageDefineDatasets::Create( wxWizard* parent )
{
////@begin wxWizPageDefineDatasets creation
    wxBitmap wizardBitmap(wxNullBitmap);
    wxWizardPageDataSim::Create( parent, NULL, NULL, wizardBitmap );

    CreateControls();
    if (GetSizer())
        GetSizer()->Fit(this);
////@end wxWizPageDefineDatasets creation
    return true;
}


/*!
 * wxWizPageDefineDatasets destructor
 */

wxWizPageDefineDatasets::~wxWizPageDefineDatasets()
{
////@begin wxWizPageDefineDatasets destruction
////@end wxWizPageDefineDatasets destruction
}


/*!
 * Member initialisation
 */

void wxWizPageDefineDatasets::Init()
{
////@begin wxWizPageDefineDatasets member initialisation
    txtDatasetCount = NULL;
    lstDatasets = NULL;
    chkStandardPedigreeHeader = NULL;
    chkWriteBinary = NULL;
    lstDiseaseModels = NULL;
////@end wxWizPageDefineDatasets member initialisation
}


/*!
 * Control creation for wxWizPageDefineDatasets
 */

void wxWizPageDefineDatasets::CreateControls()
{    
////@begin wxWizPageDefineDatasets content construction
    wxWizPageDefineDatasets* itemWizardPageDataSim1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemWizardPageDataSim1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer4 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer4, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer4->Add(itemBoxSizer5, 0, wxALIGN_RIGHT|wxALL, 0);

    wxStaticText* itemStaticText6 = new wxStaticText( itemWizardPageDataSim1, wxID_STATIC, _("Number of Data Sets:"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    itemBoxSizer5->Add(itemStaticText6, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtDatasetCount = new wxTextCtrl( itemWizardPageDataSim1, ID_TXT_DATASET_COUNT, _("1000"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer5->Add(txtDatasetCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    lstDatasets = new wxListCtrl( itemWizardPageDataSim1, ID_LST_DATASETS, wxDefaultPosition, wxSize(300, 100), wxLC_REPORT );
    itemBoxSizer4->Add(lstDatasets, 1, wxGROW|wxALL, 5);

    chkStandardPedigreeHeader = new wxCheckBox( itemWizardPageDataSim1, ID_CHK_STD_HEADER, _("Use Standard Pedigree Header"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    chkStandardPedigreeHeader->SetValue(true);
    itemBoxSizer4->Add(chkStandardPedigreeHeader, 0, wxALIGN_RIGHT|wxALL, 5);

    chkWriteBinary = new wxCheckBox( itemWizardPageDataSim1, ID_CHK_BINARY_DATASETS, _("Write Binary Datasets"), wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT );
    chkWriteBinary->SetValue(false);
    chkWriteBinary->SetHelpText(_("Save space with binary datasets, but only if your application can read them"));
    if (wxWizPageDefineDatasets::ShowToolTips())
        chkWriteBinary->SetToolTip(_("Save space with binary datasets, but only if your application can read them"));
    itemBoxSizer4->Add(chkWriteBinary, 0, wxALIGN_RIGHT|wxALL, 5);

    wxBoxSizer* itemBoxSizer11 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer4->Add(itemBoxSizer11, 0, wxGROW|wxALL, 0);

    itemBoxSizer4->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer13 = new wxBoxSizer(wxVERTICAL);
    itemBoxSizer3->Add(itemBoxSizer13, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton14 = new wxButton( itemWizardPageDataSim1, ID_CMD_ADD_CC, _("Add &Case/Control"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton14, 0, wxGROW|wxALL, 5);

    wxButton* itemButton15 = new wxButton( itemWizardPageDataSim1, ID_CMD_ADD_PED, _("Add &Pedigree"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton15, 0, wxGROW|wxALL, 5);

    wxButton* itemButton16 = new wxButton( itemWizardPageDataSim1, ID_DEL_DATASET, _("&Delete Dataset"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer13->Add(itemButton16, 0, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer17Static = new wxStaticBox(itemWizardPageDataSim1, wxID_ANY, _("Status Model"));
    wxStaticBoxSizer* itemStaticBoxSizer17 = new wxStaticBoxSizer(itemStaticBoxSizer17Static, wxHORIZONTAL);
    itemBoxSizer2->Add(itemStaticBoxSizer17, 0, wxGROW|wxALL, 5);

    wxArrayString lstDiseaseModelsStrings;
    lstDiseaseModels = new wxChoice( itemWizardPageDataSim1, ID_CMB_MODEL, wxDefaultPosition, wxDefaultSize, lstDiseaseModelsStrings, 0 );
    itemStaticBoxSizer17->Add(lstDiseaseModels, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxButton* itemButton19 = new wxButton( itemWizardPageDataSim1, ID_CMD_CONFIG_MODEL, _("Confi&gure"), wxDefaultPosition, wxDefaultSize, 0 );
    itemStaticBoxSizer17->Add(itemButton19, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

////@end wxWizPageDefineDatasets content construction
	InitListBox();
	RefreshSize();
	InitializeDiseaseModels();
}

void wxWizPageDefineDatasets::InitListBox() {
	lstDatasets->ClearAll();
	lstDatasets->InsertColumn(0, wxT("Label"));
	lstDatasets->InsertColumn(1, wxT("Type"));
	lstDatasets->InsertColumn(2, wxT("Details"));
}

void wxWizPageDefineDatasets::InitializeDiseaseModels() {
	if (modelLibrary == NULL) 
		modelLibrary = ModelList::Instance();
	int fileCount = modelLibrary->GetFileCount();
	lstDiseaseModels->Insert(_("New Disease"), 0);
	lstDiseaseModels->SetSelection(0);

	for (int i=0; i<fileCount; i++) {
		ModelList::ModelFileItem &item = modelLibrary->GetFile(i);
		lstDiseaseModels->Insert(item.label.c_str(), i+1);
	}
}

void wxWizPageDefineDatasets::RefreshSize() {
	int width, height;
	if (lstDatasets) {
		lstDatasets->GetSize(&width, &height);
		lstDatasets->SetColumnWidth(0, (int)(width * 0.2));
		lstDatasets->SetColumnWidth(1, (int)(width * 0.2));
		lstDatasets->SetColumnWidth(2,(int)(width * 0.58));
	}
}	

/*!
 * Should we show tooltips?
 */

bool wxWizPageDefineDatasets::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxWizPageDefineDatasets::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxWizPageDefineDatasets bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxWizPageDefineDatasets bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxWizPageDefineDatasets::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxWizPageDefineDatasets icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxWizPageDefineDatasets icon retrieval
}




/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_CC
 */

void wxWizPageDefineDatasets::OnCmdAddCcClick( wxCommandEvent& event )
{
	wxDlgConfigureCCDatasets dlg(this, wxID_ANY, _T("Configure Case Control Dataset"));
	
	if (dlg.ShowModal()) {
		string label = dlg.GetModelLabel();
		string modelDetails = dlg.GetModelDetails();
		cout<<"Details: "<<modelDetails<<"\n";	
		string modelSummary = dlg.GetModelSummary();
		cout<<"Summary: "<<modelSummary<<"\n";
		long idx = lstDatasets->InsertItem(datasets.size(), wxT(label.c_str()));
		lstDatasets->SetItem(idx, 0, wxT(dlg.GetModelLabel().c_str()));
		lstDatasets->SetItem(idx, 1, wxT("C/C"));
		lstDatasets->SetItem(idx, 2, wxT(modelSummary.c_str()));
		datasets.push_back(DatasetDetails(label.c_str(), modelDetails.c_str()));
	}
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_SELECTED event handler for ID_LST_DATASETS
 */

void wxWizPageDefineDatasets::OnLstDatasetsSelected( wxListEvent& event )
{
	selectedModel = event.GetIndex();
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_DESELECTED event handler for ID_LST_DATASETS
 */

void wxWizPageDefineDatasets::OnLstDatasetsDeselected( wxListEvent& event )
{
	selectedModel = -1; 
}


/*!
 * wxEVT_COMMAND_LIST_ITEM_ACTIVATED event handler for ID_LST_DATASETS
 */

void wxWizPageDefineDatasets::OnLstDatasetsItemActivated( wxListEvent& event )
{
	selectedModel = event.GetIndex();

	wxListItem info;
	info.m_itemId = selectedModel;
	info.m_col = 1;
	info.m_mask=wxLIST_MASK_TEXT;

	if (lstDatasets->GetItem(info)) {
		string type = info.m_text.c_str();
	
		if (type == "PED") {
			wxDlgConfigurePedDatasets dlg(this, wxID_ANY, _T("Edit Pedigree Dataset"));
			dlg.ConfigureModelDetails(datasets[selectedModel].details.c_str());
			cout<<"Model Details: "<<datasets[selectedModel].details.c_str()<<"\n";
			if (dlg.ShowModal()) {
				string label = dlg.GetModelLabel();
				string modelDetails = dlg.GetModelDetails();
				string modelSummary = dlg.GetModelSummary();
				lstDatasets->SetItem(selectedModel, 0, wxT(label.c_str()));
				lstDatasets->SetItem(selectedModel, 2, wxT(modelSummary.c_str()));
				datasets[selectedModel] =  DatasetDetails(label.c_str(), modelDetails.c_str());
			}
			
		}
		else if (type == "C/C") 	{
			wxDlgConfigureCCDatasets dlg(this, wxID_ANY, _T("Edit Case/Control Dataset"));
			cout<<"Model Details: "<<datasets[selectedModel].details.c_str()<<"\n";
			dlg.ConfigureModelDetails(datasets[selectedModel].details.c_str());
			if (dlg.ShowModal()) {
				string label = dlg.GetModelLabel();
				string modelDetails = dlg.GetModelDetails();
				string modelSummary = dlg.GetModelSummary();
				lstDatasets->SetItem(selectedModel, 0, wxT(label.c_str()));
				lstDatasets->SetItem(selectedModel, 2, wxT(modelSummary.c_str()));
				datasets[selectedModel] =  DatasetDetails(label.c_str(), modelDetails.c_str());
			}
		}
		else {
			wxMessageBox(_("Ooops!"));
			cout<<"Type: "<<type<<"\n";
		}
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_ADD_PED
 */

void wxWizPageDefineDatasets::OnCmdAddPedClick( wxCommandEvent& event )
{
	wxDlgConfigurePedDatasets dlg(this, wxID_ANY, _T("Configure Pedigree Dataset"));
	if (dlg.ShowModal()) {
		string label = dlg.GetModelLabel();
		string modelDetails = dlg.GetModelDetails();
		string modelSummary = dlg.GetModelSummary();
		long idx = lstDatasets->InsertItem(datasets.size(), wxT(label.c_str()));
		lstDatasets->SetItem(idx, 0, wxT(dlg.GetModelLabel().c_str()));
		lstDatasets->SetItem(idx, 1, wxT("PED"));
		lstDatasets->SetItem(idx, 2, wxT(modelSummary.c_str()));
		datasets.push_back(DatasetDetails(label.c_str(), modelDetails.c_str()));
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_DEL_DATASET
 */

void wxWizPageDefineDatasets::OnDelDatasetClick( wxCommandEvent& event )
{
	int selection = selectedModel;

	//We can't let them delete "New Disease model"
	if (selection >= 0) {
		lstDatasets->DeleteItem(selection);
		vector<DatasetDetails>::iterator itr = datasets.begin();
		vector<DatasetDetails>::iterator end = datasets.end();

		while (selection-- > 0 && itr != end) 
			itr++;
		if (itr!=end)
			datasets.erase(itr);
	
		//Drop the selected value
		selectedModel--;
	}
}


/*!
 * wxEVT_COMMAND_CHOICE_SELECTED event handler for ID_CMB_MODEL
 */
void wxWizPageDefineDatasets::OnCmbModelSelected( wxCommandEvent& event )
{
	wxDlgConfigureDiseaseModel dlg(this, wxID_ANY, _("Model Configuration"));
	int selection = lstDiseaseModels->GetSelection();
	if (selection == 0) {
		if (dlg.ShowModal()) {
			//Check to see if we need to add this to the list
			if (dlg.newFile) {
				string label = dlg.GetLabel().c_str();
				selection = modelLibrary->AddFile(label.c_str(), dlg.filetype.c_str(), dlg.filename.c_str(), true, dlg.GetLocusCount());
				lstDiseaseModels->Insert(_(label.c_str()), selection);
				lstDiseaseModels->SetSelection(selection);
			}
		}
	}
	else {
		string filename = modelLibrary->GetFilename(selection - 1).c_str();
		cout<<"Filename: "<<filename<<"\n";
		if (!wxFileExists(_T(filename.c_str()))) {
				ModelList::ModelFileItem &item = modelLibrary->GetFile(selection-1);
				wxString fileType = _T("Penetrance Files (*.pen)|*.pen;*.PEN|All Files (*.*)|*.*");
				if (item.modelType == "SIMPEN")
					fileType = _T("SIMPEN Files (*.simpen)|*.simpen|All Files (*.*)|*.*");
				else if (item.modelType == "SIMLA")
					fileType = _T("SIMLA Files (*.simla)|*.simla|All Files (*.*)|*.*");

				wxFileDialog modelFilename(this, _(item.modelType.c_str()), _("The file describing the disease model selected could not be found. Please select the appropriate file"), _(item.modelType.c_str()), fileType, wxOPEN);
				if (modelFilename.ShowModal() == wxID_OK) {
					wxString newFilename = modelFilename.GetPath().c_str();
					cout<<"You selected "<<newFilename<<"\n";
					item.filename = newFilename.c_str();
					modelLibrary->Update(selection-1, item);
					//I need to fix the model's reference in the library
				}

		}
	}
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ID_CMD_CONFIG_MODEL
 */

void wxWizPageDefineDatasets::OnCmdConfigModelClick( wxCommandEvent& event )
{
	wxDlgConfigureDiseaseModel dlg(this, wxID_ANY, _("Model Configuration"));
	int selection = lstDiseaseModels->GetSelection();
	if (selection > 0) {
		ModelList::ModelFileItem &item = modelLibrary->GetFile(selection-1);
		dlg.LoadModel(item.modelType.c_str(), item.filename.c_str(), item.editable);
	}

	if (dlg.ShowModal()) {
		//Check to see if we need to add this to the list
		if (dlg.newFile) {
			string label = dlg.GetLabel().c_str();
			selection = modelLibrary->AddFile(label.c_str(), dlg.filetype.c_str(), dlg.filename.c_str(), true, dlg.GetLocusCount());
			lstDiseaseModels->Insert(_(label.c_str()), selection);
			lstDiseaseModels->SetSelection(selection);
		}
	}
}

/*!
 * wxEVT_SIZE event handler for ID_PAGE_DEFINE_DATASETS
 */

void wxWizPageDefineDatasets::OnSize( wxSizeEvent& event )
{
////@begin wxEVT_SIZE event handler for ID_PAGE_DEFINE_DATASETS in wxWizPageDefineDatasets.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_SIZE event handler for ID_PAGE_DEFINE_DATASETS in wxWizPageDefineDatasets. 
}


bool wxWizPageDefineDatasets::TransferDataFromWindow() {
	if (lstDiseaseModels->GetSelection() <= 0) {
		wxMessageBox(wxT("You must choose your model first"), wxT("No Model Selected"), wxICON_WARNING | wxOK, this);
		return false;
	}
	if (datasets.size() <= 0) {
		wxMessageBox(wxT("You must define one or more datasets"), wxT("No Datasets Configured"), wxICON_WARNING | wxOK, this);
		return false;
	}
	Commit();
	return true;
}

/*!
 * wxEVT_WIZARD_PAGE_CHANGING event handler for ID_PAGE_DEFINE_DATASETS
 */

void wxWizPageDefineDatasets::OnPageDefineDatasetsPageChanging( wxWizardEvent& event )
{
	if (event.GetDirection()) {
		int selection = lstDiseaseModels->GetSelection();
		ModelList::ModelFileItem &item = modelLibrary->GetFile(selection-1);
		bool modelLoaded = false;

		while (!modelLoaded ){
			try {
				((wxWizPageSelectLoci*)GetNext())->SetModel(item);
				modelLoaded = true;
			}
 			catch (Utility::Exception::FileNotFound& e) {
				//Ask for the user to find the file
				wxFileDialog modelFilename(this, _(item.modelType.c_str()), _("The file describing the disease model selected could not be found. Please select the appropriate file"), _(item.modelType.c_str()), _("Penetrance Files (*.pen)|*.pen;*.PEN|SIMPEN Files (*.simpen)|*.simpen|SIMLA Files (*.simla)|*.simla|All Files (*.*)|*.*"), wxOPEN);
				if (modelFilename.ShowModal() == wxID_OK) {
					wxString newFilename = modelFilename.GetPath();
					cout<<"You selected "<<newFilename<<"\n";
					item.filename = newFilename.c_str();
					((wxWizPageSelectLoci*)GetNext())->SetModel(item);
					modelLibrary->Update(selection-1, item);
					//I need to fix the model's reference in the library
				}
				else
					return;
			} catch (Utility::Exception::General& e) {
				wxMessageBox(e.GetErrorMessage().c_str());
				event.Veto();
				return;
			}		
		}


	}
}


void wxWizPageDefineDatasets::PrepSummary(wxWizardPageDataSummary *summary) {
	static string yesNo[]={"NO","YES"};
	summary->Clear("Dataset Generation Summary");
	summary->AddHeader("Dataset Definitions");
	stringstream ss;
	int detailsCount = datasets.size();
	ss<<"<TABLE border=1><TR bgcolor=\"#dddddd\"><TH>Dataset Label</TH><TH>Type</TH><TH>Details</TH></TR>";
	for (int i=0; i<detailsCount; i++) {
		DatasetDetails& curDataset = datasets[i];
		ss<<"<TR><TD>"<<curDataset.label<<"</TD><TD> </TD><TD><PRE>"<<curDataset.details<<"</PRE></TD></TR>";
	}
	ss<<"</TABLE>";
	ss<<"<P><TABLE border=1><TR><TH  bgcolor=\"#dddddd\">Datasets Per Type</TH> <TD>"<<ExtractInteger(txtDatasetCount)<<"</TD> </TR>";
	ss<<"          <TR><TH bgcolor=\"#dddddd\">Standard Pedigree Header</TH><TD>"<<yesNo[chkStandardPedigreeHeader->GetValue()]<<"</TD></TR>";
	ss<<"<TR><TH bgcolor=\"#dddddd\">Binary Datasets</TH><TD>"<<yesNo[chkWriteBinary->GetValue()]<<"</TD></TR>";
	ss<<"<TR><TH bgcolor=\"#dddddd\">Disease Model Selected</TH><TD>";
	int selection = lstDiseaseModels->GetSelection();
	assert(selection > 0);
	ModelList::ModelFileItem &item = modelLibrary->GetFile(selection-1);
	ss<<item.label<<"</TD></TR>";
	ss<<"<TR><TH bgcolor=\"#dddddd\">Model Filename</TH><TD>"<<item.filename<<"</TD></TR></TABLE>";
	
	summary->AddNote(ss.str().c_str());

}

void wxWizPageDefineDatasets::PrepConfig(ostream& config) {
	static string yesNo[]={"NO","YES"};
	config<<"\n#Dataset Definitions\n";
	config<<"DATASET_COUNT"<<ExtractInteger(txtDatasetCount)<<"\n";
	int detailsCount = datasets.size();
	for (int i=0; i<detailsCount; i++) 
		config<<datasets[i].details<<"\n";

	config<<"\n#Standard Pedigree Header. Yes (10) columns. No(6)\n";
	config<<"USE_STD_PEDIGREE_HEADER\t"<<yesNo[chkStandardPedigreeHeader->GetValue()]<<"\n";

	config<<"\n#Binary Datasets are currently only supported by plato\n";
	config<<"BINARY_DATASETS        \t"<<yesNo[chkWriteBinary->GetValue()]<<"\n";
	
}


/**
 * @brief This is called prior to saving the configuration
 */
void wxWizPageDefineDatasets::Commit() {
	appController->parameters.configuration->ClearDatasets();
	int detailsCount = datasets.size();
	for (int i=0; i<detailsCount; i++) {
		stringstream lines(datasets[i].details.c_str());
		char line[4096];
		while (!lines.eof()) {
			lines.getline(line, 4096);
			appController->parameters.configuration->AddDataset(line);
		}
	}	
	appController->parameters.configuration->datasetSettings.simSets = ExtractInteger(txtDatasetCount);
	appController->parameters.configuration->datasetSettings.binaryDatasets = chkWriteBinary->GetValue();
	Individual::StandardPedigreeHeader = chkStandardPedigreeHeader->GetValue();

}

/**
 * @brief This is called just after a configuration has been loaded
 */
void wxWizPageDefineDatasets::RefreshSettings() {
	chkWriteBinary->SetValue(appController->parameters.configuration->datasetSettings.binaryDatasets);
	chkStandardPedigreeHeader->SetValue(Individual::StandardPedigreeHeader);
	
	UpdateTextField(txtDatasetCount, (int)appController->parameters.configuration->datasetSettings.simSets);

	vector<Sample*>::iterator itr = appController->parameters.configuration->datasetSettings.samples.begin();
	vector<Sample*>::iterator end = appController->parameters.configuration->datasetSettings.samples.end();

	while (itr != end) {
		Sample *sample = *itr;
		long idx = lstDatasets->InsertItem(datasets.size(), wxT(sample->GetLabel().c_str()));
		//long idx = lstDatasets->InsertItem(datasets.size(), wxT(itr->GetLabel().c_str()));
		lstDatasets->SetItem(idx, 0, wxT(sample->GetLabel().c_str()));
		lstDatasets->SetItem(idx, 1, wxT(sample->GetType().c_str()));
		lstDatasets->SetItem(idx, 2, wxT(sample->GetSummary().c_str()));
		datasets.push_back(DatasetDetails(sample->GetLabel().c_str(), sample->GetDetails().c_str()));
		itr++;
	}
}
/**
 * @brief Compares current state with that of the saved version at least against what is in memory
 * @return True indicates that something has changed 
 */
bool wxWizPageDefineDatasets::HasChanged() {
	return true;
}

}

}





