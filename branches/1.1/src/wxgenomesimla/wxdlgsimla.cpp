/////////////////////////////////////////////////////////////////////////////
// Name:        wxdlgsimla.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 15:55:49 CST
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

#include "wxdlgsimla.h"
#include <fstream>

////@begin XPM images
////@end XPM images
namespace GenomeSIM {

namespace GUI {


/*!
 * wxDlgSIMLA type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxDlgSIMLA, wxPanel )


/*!
 * wxDlgSIMLA event table definition
 */

BEGIN_EVENT_TABLE( wxDlgSIMLA, wxPanel )

////@begin wxDlgSIMLA event table entries
    EVT_SIZE( wxDlgSIMLA::OnSize )

////@end wxDlgSIMLA event table entries

	EVT_SPINCTRL(ID_SPINCTRL, wxDlgSIMLA::OnDialLocusSpin)
END_EVENT_TABLE()



/*!
 * wxDlgSIMLA constructors
 */
wxDlgSIMLA::wxDlgSIMLA() : locGrdMaster(NULL), intGrdMaster(NULL)
{
    Init();
}

wxDlgSIMLA::wxDlgSIMLA( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : locGrdMaster(NULL), intGrdMaster(NULL)
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxDlgSIMLA creator
 */

bool wxDlgSIMLA::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxDlgSIMLA creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxDlgSIMLA creation
    return true;
}


/*!
 * wxDlgSIMLA destructor
 */

wxDlgSIMLA::~wxDlgSIMLA()
{
////@begin wxDlgSIMLA destruction
////@end wxDlgSIMLA destruction
}


/*!
 * Member initialisation
 */

void wxDlgSIMLA::Init()
{
////@begin wxDlgSIMLA member initialisation
    dialLocusCount = NULL;
    txtTargetPrevalence = NULL;
    grdLocusDefinition = NULL;
    grdInteractions = NULL;
////@end wxDlgSIMLA member initialisation
}


/*!
 * Control creation for wxDlgSIMLA
 */

void wxDlgSIMLA::CreateControls()
{    
////@begin wxDlgSIMLA content construction
    wxDlgSIMLA* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 0, wxGROW|wxALL, 0);

    wxStaticText* itemStaticText4 = new wxStaticText( itemPanel1, wxID_STATIC, _("Locus Count:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer3->Add(itemStaticText4, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    dialLocusCount = new wxSpinCtrl( itemPanel1, ID_SPINCTRL, _T("2"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 6, 2 );
    itemBoxSizer3->Add(dialLocusCount, 0, wxGROW|wxALL, 5);

    itemBoxSizer3->Add(5, 5, 1, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText7 = new wxStaticText( itemPanel1, wxID_STATIC, _("Disease Prevalence:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer3->Add(itemStaticText7, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtTargetPrevalence = new wxTextCtrl( itemPanel1, ID_TARGETPREV, _("0.001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer3->Add(txtTargetPrevalence, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer2->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxBoxSizer* itemBoxSizer10 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer10, 1, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer11Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Main Effects And Locus Definitions:"));
    wxStaticBoxSizer* itemStaticBoxSizer11 = new wxStaticBoxSizer(itemStaticBoxSizer11Static, wxHORIZONTAL);
    itemBoxSizer10->Add(itemStaticBoxSizer11, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer12 = new wxBoxSizer(wxVERTICAL);
    itemStaticBoxSizer11->Add(itemBoxSizer12, 1, wxGROW|wxALL, 5);

    grdLocusDefinition = new wxGrid( itemPanel1, ID_GRID, wxDefaultPosition, wxSize(200, 100), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdLocusDefinition->SetDefaultColSize(50);
    grdLocusDefinition->SetDefaultRowSize(25);
    grdLocusDefinition->SetColLabelSize(25);
    grdLocusDefinition->SetRowLabelSize(50);
    grdLocusDefinition->CreateGrid(2, 2, wxGrid::wxGridSelectCells);
    itemBoxSizer12->Add(grdLocusDefinition, 1, wxGROW|wxALL, 5);

    itemBoxSizer10->Add(5, 5, 0, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer15Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Interactions:"));
    wxStaticBoxSizer* itemStaticBoxSizer15 = new wxStaticBoxSizer(itemStaticBoxSizer15Static, wxVERTICAL);
    itemBoxSizer10->Add(itemStaticBoxSizer15, 1, wxGROW|wxALL, 5);

    grdInteractions = new wxGrid( itemPanel1, ID_GRID1, wxDefaultPosition, wxSize(150, 100), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdInteractions->SetDefaultColSize(50);
    grdInteractions->SetDefaultRowSize(25);
    grdInteractions->SetColLabelSize(25);
    grdInteractions->SetRowLabelSize(50);
    grdInteractions->CreateGrid(2, 2, wxGrid::wxGridSelectCells);
    itemStaticBoxSizer15->Add(grdInteractions, 1, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText17 = new wxStaticText( itemPanel1, wxID_STATIC, _("Values of 1.0 are allowed and indicate that no\n effect is to be modelled between these loci"), wxDefaultPosition, wxDefaultSize, wxST_NO_AUTORESIZE );
    itemStaticBoxSizer15->Add(itemStaticText17, 0, wxGROW|wxALL, 0);

////@end wxDlgSIMLA content construction
	
	InitGrids();
}


/*!
 * Should we show tooltips?
 */

bool wxDlgSIMLA::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxDlgSIMLA::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxDlgSIMLA bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxDlgSIMLA bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxDlgSIMLA::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxDlgSIMLA icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxDlgSIMLA icon retrieval
}

bool wxDlgSIMLA::Evaluate() {
	bool success = false;
/*	stringstream ss(modelCfg.c_str());

	string cmd;
	
	ss>>cmd;			//>>percentage;

	if (model)  {
		cout<<"Overwriting previous model definition. (";
		model->GenerateReport(cout, 0);
		cout<<")\n";
		delete model;
	}

	//Right now, all models are based on penetranceModels- which is a factory
	model = Simulation::StatusModel::PenetranceModel::GetModel(ss, poolMgr);

	if (requireModel && model==NULL) {
		if (modelCfg=="")
			cout<<"There was no model specified. Please set up a disease model in the configuration file\n";
		else
			cout<<"Invalid model definition: "<<modelCfg<<". Does this file really exist?\n";
		abort();
	}
	return model;*/
	return success;
}

void wxDlgSIMLA::InitGrids() {
	grdLocusDefinition->SetDefaultColSize(75);
	grdLocusDefinition->SetColLabelSize(50);
	grdLocusDefinition->SetRowLabelSize(35);
	
	grdInteractions->SetDefaultColSize(75);
	grdInteractions->SetColLabelSize(25);
	grdInteractions->SetRowLabelSize(80);
	

	if (locGrdMaster)
		delete locGrdMaster;
	locGrdMaster = new wxGridSimlaLoci(dialLocusCount->GetValue());
	grdLocusDefinition->SetTable(locGrdMaster, true);

	locGrdMaster->SetLocusCount(dialLocusCount->GetValue());

	if (intGrdMaster)
		delete intGrdMaster;
	intGrdMaster = new wxGridSimlaInteractions(dialLocusCount->GetValue());
	grdInteractions->SetTable(intGrdMaster, true);

	intGrdMaster->SetLocusCount(dialLocusCount->GetValue());

	
}

bool wxDlgSIMLA::Save(const char *filename, const char *desc) {
	//Save interactions
	if (locGrdMaster->Verify() && intGrdMaster->Verify()) {
		ofstream file(filename);
	
		file<<SIMLA_PREVALENCE<<"\t"<<ExtractDouble(txtTargetPrevalence)<<"\n";
		locGrdMaster->Save(file);
		intGrdMaster->Save(file);	
		return true;
	}
	
	//intGrdMaster->WriteToStream(file);
	return false;
}

void wxDlgSIMLA::OnDialLocusSpin(wxSpinEvent& event) {
	locGrdMaster->SetLocusCount(dialLocusCount->GetValue());
//	locGrdMaster->Refresh();
	intGrdMaster->SetLocusCount(dialLocusCount->GetValue());
//	intGrdMaster->Refresh();
}

void wxGridBetaWise::SetLocusCount(int locusCount) {
	this->locusCount=locusCount;
	ResizeContents(locusCount);
}


int wxGridBetaWise::GetLocusCount() {
	return locusCount;
}


void wxGridSimlaInteractions::Load(FileToMap& file, int locCount) {
	SetLocusCount(locCount);
	vector<string> interactionDetails;
	
	if (file.GetLines(SIMLA_INTERACTION, interactionDetails)) {
		assert(interactionDetails.size() <= interactions.size() * 2);

		for (size_t i =0; i<interactionDetails.size(); ) {
			string id;
			double beta;
			id = interactionDetails[i++];
			beta = atof(interactionDetails[i++].c_str());
			if (interactions.find(id) != interactions.end())
				interactions[id] = Interactions(id, beta);
		}
	}
}

void wxGridSimlaInteractions::LoadFromStream(istream &stream) {
	//Clear previous values from the map
	map<string, Interactions>::iterator itr = interactions.begin();
	map<string, Interactions>::iterator end = interactions.end();
	
	while (itr != end) {
		itr->second.beta = 0.0;
		itr++;
	}	

	while (!stream.eof()) {
		string id;
		double beta;

		stream>>id>>beta;
		if (interactions.find(id) != interactions.end())
			interactions[id] = Interactions(id, beta);
	}
}

bool wxGridSimlaInteractions::Verify() {
	map<string, Interactions>::iterator itr = interactions.begin();
	map<string, Interactions>::iterator end = interactions.end();

	while (itr != end) {
		if (itr->second.beta < 0.00001) {
			wxMessageBox(_("Rel. Risk much be greater than 0"), _("Invalid Input"), wxICON_WARNING|wxOK);
			return false;
		}
		itr++;	
	}
	return true;
}

void wxGridSimlaInteractions::Save(ostream& os) {
	map<string, Interactions>::iterator itr = interactions.begin();
	map<string, Interactions>::iterator end = interactions.end();

	while (itr != end) {
		if (itr->second.beta != 1.0) 
			os<<SIMLA_INTERACTION<<"\t"<<itr->second.ID<<"\t"<<itr->second.beta<<"\n";
		itr++;	
	}
}
	
void wxGridSimlaInteractions::WriteToStream(ostream &stream) {
	map<string, Interactions>::iterator itr = interactions.begin();
	map<string, Interactions>::iterator end = interactions.end();

	while (itr != end) {
		if (itr->second.beta > 0.00001) 
			stream<<itr->second.ID<<"\t"<<itr->second.beta<<"\n";
		itr++;	
	}
}

void wxGridSimlaLoci::InitCheckBoxes() {
	int curSize = loci.size();
	wxGrid *grid = GetView();

	if (grid) {
		for (int i=0; i<curSize; i++){
			grid->SetCellRenderer(i, 0, new wxGridCellBoolRenderer);
			grid->SetCellEditor(i, 0, new wxGridCellBoolEditor);
		}
	}
	
}


void wxGridSimlaLoci::ResizeContents(int locusCount) {
	int curSize = loci.size();

	if (GetView() && curSize > 0) {
		wxGridTableMessage msg1(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, curSize);
		GetView()->ProcessTableMessage(msg1);	
	}

	loci.resize(locusCount);
	if (curSize < locusCount) {
		char lbl[]="A";
		for (int i=0; i<locusCount; i++) {
			loci[i].ID = lbl;
			lbl[0]++;
		}
	}
	this->locusCount = locusCount;
	
	
	if (GetView()) {
		wxGridTableMessage msg2(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, loci.size());
		GetView()->ProcessTableMessage(msg2);
		InitCheckBoxes();
	}
}



void wxGridSimlaLoci::Refresh() {
	wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 0);
	GetView()->ProcessTableMessage(msg);
}

void wxGridSimlaLoci::Reset() {
	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, locusCount);
		GetView()->ProcessTableMessage(msg);
	}
	loci.clear();
}
double wxGridSimlaLoci::GetValueAsDouble(int row, int col) {
	double value = 0.0;
	if (row < locusCount)  {
		if (col == 1)
			value = loci[row].beta;
		else if (col == 2)
			value = loci[row].modelType;
	}
	return value;
}

bool wxGridSimlaLoci::GetValueAsBool(int row, int col) {
	bool value = false;
	if (row < locusCount && col == 0)
		value = loci[row].diseaseAtMajor;
	return value;
}

void wxGridSimlaLoci::SetValueAsBool(int row, int col, bool isSet) {
	if (row < locusCount && col == 0)
		loci[row].diseaseAtMajor = isSet;
}


void wxGridSimlaLoci::SetValueAsDouble(int row, int col, double value) {
	if (row < locusCount && col == 1) 
		loci[row].beta = value;
	else if (row < locusCount && col == 2)
		loci[row].modelType = value;
	
		
}


bool wxGridSimlaLoci::CanGetValueAs(int row, int col, const wxString& typeName) {
	if (typeName == wxGRID_VALUE_BOOL && col == 0)
		return true;
	if (typeName == wxGRID_VALUE_FLOAT && col == 1)
		return true;
	if (typeName == wxGRID_VALUE_FLOAT && col == 2)
		return true;

	return false;
}

bool wxGridSimlaLoci::CanSetValueAs(int row, int col, const wxString& typeName) {
	return CanGetValueAs(row, col, typeName);
}

wxString wxGridSimlaLoci::GetTypeName(int row, int col) {
	return wxGRID_VALUE_FLOAT;
}

wxString wxGridSimlaLoci::GetRowLabelValue(int row) {
	if (row<locusCount)
		return loci[row].ID.c_str();
	return "";
}

wxString wxGridSimlaLoci::GetColLabelValue(int col) {
	if (col == 0)	
		return "Disease\nAt Major";
	else if (col==1)
		return "Rel.\nRisk";
	else if (col==2)
		return "Model\nType";
	else
		return "";
}

void wxGridSimlaLoci::SetValue(int row, int col, const wxString& value) {
	if (row<locusCount) {
		LocusDetails &details = loci[row];
		if (col == 0)	
			details.diseaseAtMajor = atoi(value.c_str());
		else if (col == 1) 
			details.beta = atof(value.c_str());
		else if (col == 2)
			details.modelType = atof(value.c_str());
	}
}

wxString wxGridSimlaLoci::GetValue(int row, int col) {
	LocusDetails& details = loci[row];

	wxString val;
	if (col == 0)
		val.Printf("%d", (int)details.diseaseAtMajor);
	else if (col == 1)
		val.Printf("%f", details.beta);
	else if (col == 2)
		val.Printf("%.2f", details.modelType);

	return val;
}

bool wxGridSimlaLoci::IsEmptyCell(int row, int col) {
	if (row < locusCount && col < 3)
		return true;

	return false;
}

int wxGridSimlaLoci::GetNumberCols() {
	return 3;
}

int wxGridSimlaLoci::GetNumberRows() {
	return locusCount;
}


void wxGridSimlaLoci::Load(FileToMap &file) {
	vector<string> locusData;
	int locCount = atoi(file.GetLine("SIMLA_LOC_COUNT").c_str());
	SetLocusCount(locCount);

	file.GetLines(SIMLA_LOCUS, locusData);
	int fragCount = locusData.size();
	assert(fragCount * 0.25 == (size_t)locCount);

	loci.clear();
	for (int i=0; i<fragCount; ) {
		string ID = locusData[i++];
		bool diseaseAtMajor = atoi(locusData[i++].c_str());
		double beta = atof(locusData[i++].c_str());
		double type = atof(locusData[i++].c_str());
		loci.push_back(LocusDetails(ID, diseaseAtMajor, beta, type));
	}
	
}

bool wxGridSimlaLoci::Verify() {
	vector<LocusDetails>::iterator itr = loci.begin();
	vector<LocusDetails>::iterator end = loci.end();

	while (itr != end) {
		if (itr->beta < 0.00001) {
			wxMessageBox(_("Rel. Risk muist be greater than 0"), _("Invalid Input"), wxICON_WARNING|wxOK);
			return false;
		}
		else if (itr->modelType > 1.0 || itr->modelType < 0.0) {
			wxMessageBox(_("The Model Type must be between 0.0 and 1.0 (inclusive)"), _("Invalid Input"), wxICON_WARNING|wxOK);
			return false;
		}
		itr++;	
	}
	return true;
}

void wxGridSimlaLoci::Save(ostream& os) {
	os<<"SIMLA_LOC_COUNT\t"<<loci.size()<<"\n";
		
	vector<LocusDetails>::iterator itr = loci.begin();
	vector<LocusDetails>::iterator end = loci.end();

	while (itr != end) {
		os<<SIMLA_LOCUS<<"\t"<<itr->ID<<"\t"<<itr->diseaseAtMajor<<"\t"<<itr->beta<<"\t"<<itr->modelType<<"\n";	
		itr++;
	}
	
}

bool wxDlgSIMLA::Import(const char *filename, string &desc) {
	FileToMap file;
	file.Parse(filename);
	txtTargetPrevalence->SetValue(file.GetLine(SIMLA_PREVALENCE).c_str());

	locGrdMaster->Load(file);	
	dialLocusCount->SetValue(locGrdMaster->GetLocusCount());
	intGrdMaster->Load(file, locGrdMaster->GetLocusCount());
	
	return true;
}

void wxGridSimlaInteractions::Refresh() {
	if (GetView()){
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 0);
		GetView()->ProcessTableMessage(msg);
	}
}

void wxGridSimlaInteractions::Reset() {
	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, interactions.size());
		GetView()->ProcessTableMessage(msg);
	}
	interactions.clear();
}
double wxGridSimlaInteractions::GetValueAsDouble(int row, int col) {
	double value = 0.0;
	if (row < locusCount && col == 0)
		value = interactions[interactionNames[row]].beta;
	return value;
}

void wxGridSimlaInteractions::SetValueAsDouble(int row, int col, double value) {
	if (row < locusCount && col == 0) 
		interactions[interactionNames[row]].beta = value;
		
}


bool wxGridSimlaInteractions::CanGetValueAs(int row, int col, const wxString& typeName) {
	if (typeName == wxGRID_VALUE_FLOAT && col == 01)
		return true;

	return false;
}

bool wxGridSimlaInteractions::CanSetValueAs(int row, int col, const wxString& typeName) {
	return CanGetValueAs(row, col, typeName);
}

wxString wxGridSimlaInteractions::GetTypeName(int row, int col) {
	return wxGRID_VALUE_FLOAT;
}

wxString wxGridSimlaInteractions::GetColLabelValue(int col) {
	if (col == 0)	
		return "Relative Risk";
	else
		return "";
}
wxString wxGridSimlaInteractions::GetRowLabelValue(int row) {
	if ((size_t)row<interactions.size())
		return interactions[interactionNames[row]].ID.c_str();
	return "";
}

void wxGridSimlaInteractions::SetValue(int row, int col, const wxString& value) {
	if ((size_t)row<interactions.size()) {
		Interactions &details = interactions[interactionNames[row]];
		if (col == 0) 
			details.beta = atof(value.c_str());
	}
}

wxString wxGridSimlaInteractions::GetValue(int row, int col) {
	Interactions& details = interactions[interactionNames[row]];

	wxString val;
	if (col == 0)
		val.Printf("%f", details.beta);
		
	return val;
}

bool wxGridSimlaInteractions::IsEmptyCell(int row, int col) {
	if ((size_t)row<interactions.size() && col < 1)
		return true;

	return false;
}

int wxGridSimlaInteractions::GetNumberCols() {
	return 1;
}

int wxGridSimlaInteractions::GetNumberRows() {
	return interactions.size();
}


void wxGridSimlaInteractions::Append(string prev, char letter, int curIdx) {
	//Recursive solution, if we haven't gone to far, we'll add some more
	if (curIdx < locusCount) {
		if (prev.length() > 0) {
			prev.append("x");
		}
		
		prev.append(1, letter);
		//Add local node first

		Append(prev, letter+1, curIdx+1);

		if (prev.length() > 1 && interactions.find(prev) == interactions.end())  {
			interactionNames.push_back(prev);
			interactions[prev] = Interactions(prev, 1.0);
		}

		string local;
		local.append(1, letter);

		for (int i=curIdx+1; i<locusCount; i++) {
			Append(local, letter+i, curIdx+i);
			Append("", letter+i, curIdx+i);
		}

	}
}




void wxGridSimlaInteractions::ResizeContents(int locusCount) {
	//int curSize = interactions.size();
	Reset();
	//interactions.clear();
	string base;
	interactions.clear();
	interactionNames.clear();
	Append(base, 'A', 0);
	sort(interactionNames.begin(), interactionNames.end());
	

	if (GetView() ) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, interactions.size());
		GetView()->ProcessTableMessage(msg);
	}
}



/*!
 * wxEVT_SIZE event handler for ID_WXDLGSIMLA
 */

void wxDlgSIMLA::OnSize( wxSizeEvent& event )
{
	if (grdLocusDefinition) {
		int labelWidth = grdLocusDefinition->GetColLabelSize();
		int w,h;
	
		grdLocusDefinition->GetSize(&w, &h);
		int availableWidth = w - labelWidth - 25 - 75;
		
		grdLocusDefinition->SetColSize(0, 75);
		grdLocusDefinition->SetColSize(1,  (int) (availableWidth * 0.60));
		grdLocusDefinition->SetColSize(2,  (int) (availableWidth * 0.40));
	
		grdLocusDefinition->SetColFormatFloat(1, 5, 3);
		grdLocusDefinition->SetColFormatFloat(2, 4, 2);
	
		labelWidth = grdInteractions->GetColLabelSize();
		grdInteractions->GetSize(&w, &h);
		availableWidth = w - labelWidth - 80;
	
		grdInteractions->SetColSize(0, availableWidth);
	
		event.Skip();
	}	
}
}

}




