/////////////////////////////////////////////////////////////////////////////
// Name:        wxpagepenetrancemodel.cpp
// Purpose:     
// Author:      Eric Torstenson
// Modified by: 
// Created:     Fri 22 Feb 2008 13:00:53 CST
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

#include "wxpagepenetrancemodel.h"
#include <sstream>
#include <fstream>
////@begin XPM images
////@end XPM images
#include "wx/artprov.h"
#include "utility/lineparser.h"
#include "simulation/penetranceeval.h"
#include "wxdlgbasichtmlreport.h"


namespace GenomeSIM {

namespace GUI {
#define ID_INVERT_ALLELES 151245

using namespace Utility;
using namespace Simulation;



/*!
 * wxPagePenetranceModel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( wxPagePenetranceModel, wxPanel )


/*!
 * wxPagePenetranceModel event table definition
 */

BEGIN_EVENT_TABLE( wxPagePenetranceModel, wxPanel )

////@begin wxPagePenetranceModel event table entries
    EVT_SPINCTRL( ID_LOCUS_COUNT, wxPagePenetranceModel::OnLocusCountUpdated )
    EVT_TEXT( ID_LOCUS_COUNT, wxPagePenetranceModel::OnLocusCountTextUpdated )

    EVT_GRID_LABEL_RIGHT_CLICK( wxPagePenetranceModel::OnLabelRightClick )
    EVT_GRID_EDITOR_SHOWN( wxPagePenetranceModel::OnEditorShown )
    EVT_GRID_CMD_LABEL_LEFT_CLICK( ID_GRD_ALLELE_FREQ, wxPagePenetranceModel::OnGrdAlleleFreqLabelLeftClick )
    EVT_GRID_CMD_LABEL_RIGHT_CLICK( ID_GRD_ALLELE_FREQ, wxPagePenetranceModel::OnGrdAlleleFreqLabelRightClick )

    EVT_GRID_CMD_LABEL_LEFT_CLICK( ID_GRD_PENETRANCE, wxPagePenetranceModel::OnGrdPenetranceLabelLeftClick )
    EVT_GRID_CMD_LABEL_RIGHT_CLICK( ID_GRD_PENETRANCE, wxPagePenetranceModel::OnGrdPenetranceLabelRightClick )

////@end wxPagePenetranceModel event table entries
	EVT_MENU(ID_INVERT_ALLELES, wxPagePenetranceModel::InvertAlleles )
END_EVENT_TABLE()


/*!
 * wxPagePenetranceModel constructors
 */

wxPagePenetranceModel::wxPagePenetranceModel() : penGridMaster(NULL)
{
    Init();
}

wxPagePenetranceModel::wxPagePenetranceModel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ): freqGridMaster(NULL), penGridMaster(NULL)
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * wxPagePenetranceModel creator
 */

bool wxPagePenetranceModel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin wxPagePenetranceModel creation
    wxPanel::Create( parent, id, pos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end wxPagePenetranceModel creation
    return true;
}


/*!
 * wxPagePenetranceModel destructor
 */

wxPagePenetranceModel::~wxPagePenetranceModel()
{
////@begin wxPagePenetranceModel destruction
////@end wxPagePenetranceModel destruction
}


/*!
 * Member initialisation
 */

void wxPagePenetranceModel::Init()
{
////@begin wxPagePenetranceModel member initialisation
    txtFreqThreshold = NULL;
    dialLocusCount = NULL;
    grdAlleleFrequencies = NULL;
    grdPenetrances = NULL;
////@end wxPagePenetranceModel member initialisation
	isLocked = false;
}


/*!
 * Control creation for wxPagePenetranceModel
 */

void wxPagePenetranceModel::CreateControls()
{  
	wxBitmap SaveIcon = wxArtProvider::GetBitmap(wxART_FILE_SAVE, wxART_TOOLBAR, wxSize(16,16));
  	wxBitmap SaveAsIcon = wxArtProvider::GetBitmap(wxART_FILE_SAVE_AS, wxART_TOOLBAR, wxSize(16, 16));
	wxBitmap FolderIcon = wxArtProvider::GetBitmap(wxART_FOLDER, wxART_TOOLBAR, wxSize(16, 16));
////@begin wxPagePenetranceModel content construction
    wxPagePenetranceModel* itemPanel1 = this;

    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemPanel1->SetSizer(itemBoxSizer2);

    wxBoxSizer* itemBoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer2->Add(itemBoxSizer3, 1, wxGROW|wxALL, 5);

    wxStaticBox* itemStaticBoxSizer4Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Allele Frequencies:"));
    wxStaticBoxSizer* itemStaticBoxSizer4 = new wxStaticBoxSizer(itemStaticBoxSizer4Static, wxHORIZONTAL);
    itemBoxSizer3->Add(itemStaticBoxSizer4, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer5 = new wxBoxSizer(wxVERTICAL);
    itemStaticBoxSizer4->Add(itemBoxSizer5, 1, wxGROW|wxALL, 5);

    wxBoxSizer* itemBoxSizer6 = new wxBoxSizer(wxHORIZONTAL);
    itemBoxSizer5->Add(itemBoxSizer6, 0, wxGROW|wxALL, 0);

    wxStaticText* itemStaticText7 = new wxStaticText( itemPanel1, wxID_STATIC, _("Threshold:"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer6->Add(itemStaticText7, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    txtFreqThreshold = new wxTextCtrl( itemPanel1, ID_TXT_FREQ_THRESH, _("0.001"), wxDefaultPosition, wxDefaultSize, 0, wxTextValidator(wxFILTER_NUMERIC) );
    itemBoxSizer6->Add(txtFreqThreshold, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    itemBoxSizer6->Add(5, 5, 1, wxGROW|wxALL, 5);

    dialLocusCount = new wxSpinCtrl( itemPanel1, ID_LOCUS_COUNT, _T("2"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 5, 2 );
    itemBoxSizer6->Add(dialLocusCount, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    grdAlleleFrequencies = new wxGrid( itemPanel1, ID_GRD_ALLELE_FREQ, wxDefaultPosition, wxSize(300, 100), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdAlleleFrequencies->SetDefaultColSize(50);
    grdAlleleFrequencies->SetDefaultRowSize(25);
    grdAlleleFrequencies->SetColLabelSize(25);
    grdAlleleFrequencies->SetRowLabelSize(50);
    grdAlleleFrequencies->CreateGrid(2, 2, wxGrid::wxGridSelectCells);
    itemBoxSizer5->Add(grdAlleleFrequencies, 1, wxGROW|wxALL, 5);

    wxStaticText* itemStaticText12 = new wxStaticText( itemPanel1, wxID_STATIC, _("*Frequencies for allele 1 are specifically associated with the first allele and may be either major or minor, depending on the underlying population. "), wxDefaultPosition, wxSize(25, -1), wxST_NO_AUTORESIZE );
    itemStaticText12->SetFont(wxFont(8, wxSWISS, wxITALIC, wxNORMAL, false, wxT("Sans")));
    itemBoxSizer5->Add(itemStaticText12, 0, wxGROW|wxALL, 7);

    itemBoxSizer3->Add(5, 5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 10);

    wxStaticBox* itemStaticBoxSizer14Static = new wxStaticBox(itemPanel1, wxID_ANY, _("Penetrance Values:"));
    wxStaticBoxSizer* itemStaticBoxSizer14 = new wxStaticBoxSizer(itemStaticBoxSizer14Static, wxVERTICAL);
    itemBoxSizer3->Add(itemStaticBoxSizer14, 0, wxGROW|wxALL, 5);

    grdPenetrances = new wxGrid( itemPanel1, ID_GRD_PENETRANCE, wxDefaultPosition, wxSize(250, 150), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL );
    grdPenetrances->SetDefaultColSize(50);
    grdPenetrances->SetDefaultRowSize(25);
    grdPenetrances->SetColLabelSize(25);
    grdPenetrances->SetRowLabelSize(50);
    grdPenetrances->CreateGrid(5, 5, wxGrid::wxGridSelectCells);
    itemStaticBoxSizer14->Add(grdPenetrances, 1, wxGROW|wxALL, 5);

    // Connect events and objects
    grdAlleleFrequencies->Connect(ID_GRD_ALLELE_FREQ, wxEVT_SIZE, wxSizeEventHandler(wxPagePenetranceModel::OnSize), NULL, this);
////@end wxPagePenetranceModel content construction

	InitGrids();
}


/*!
 * Should we show tooltips?
 */

bool wxPagePenetranceModel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap wxPagePenetranceModel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin wxPagePenetranceModel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end wxPagePenetranceModel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon wxPagePenetranceModel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin wxPagePenetranceModel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end wxPagePenetranceModel icon retrieval
}


/*!
 * wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_SPINCTRL1
 */

void wxPagePenetranceModel::OnLocusCountUpdated( wxSpinEvent& event )
{
////@begin wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_SPINCTRL1 in wxPagePenetranceModel.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_COMMAND_SPINCTRL_UPDATED event handler for ID_SPINCTRL1 in wxPagePenetranceModel. 
}


/*!
 * wxEVT_COMMAND_TEXT_UPDATED event handler for ID_SPINCTRL1
 */

void wxPagePenetranceModel::OnLocusCountTextUpdated( wxCommandEvent& event )
{
	if (dialLocusCount) {
		freqGridMaster->SetLocusCount(dialLocusCount->GetValue());
		penGridMaster->SetLocusCount(dialLocusCount->GetValue());
	}
}



void wxPagePenetranceModel::InitGrids() {
    grdAlleleFrequencies->SetDefaultColSize(75);
    grdAlleleFrequencies->SetColLabelSize(25);
    grdAlleleFrequencies->SetRowLabelSize(80);
	grdPenetrances->SetDefaultColSize(100);
	grdPenetrances->SetRowLabelSize(125);


	if (freqGridMaster)
		delete freqGridMaster;
	freqGridMaster = new wxFreqGrid(2);
	grdAlleleFrequencies->SetTable(freqGridMaster, true);

	if (penGridMaster)
		delete penGridMaster;
	penGridMaster = new wxPenGrid(2);
	grdPenetrances->SetTable(penGridMaster, true);
	
}



wxFreqGrid::wxFreqGrid(int modelSize) : frequencies(2), lociCount(modelSize) { 
	
}

int wxFreqGrid::GetNumberCols() {
	return 2;
}

int wxFreqGrid::GetNumberRows() { 
	return lociCount;
}


bool wxFreqGrid::IsEmptyCell( int row, int col ) {
	return false;
}

wxString wxFreqGrid::GetValue( int row, int col ) {	
	wxString rVal = wxT("");


	if ((size_t)row < frequencies.size()) {
		Freqs &l = frequencies[row];
	
//		cout<<"Fetching	"<<row<<"x"<<col<<" from SNP\n";

		switch (col) {
			case 0:
				rVal.Printf("%f", l.al1);
				break;
			case 1:			
				rVal.Printf("%f", l.al2);
				break;
			default:
				cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
				assert(0);
		}
	}
	return rVal;
}
void wxFreqGrid::Reset() {

	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, frequencies.size());
		GetView()->ProcessTableMessage(msg);
	}
	frequencies.clear();
}


void wxFreqGrid::Load(vector<string>& freqs) {
	size_t allCount = freqs.size() / 2;
	lociCount = allCount / 2;

	/*cout<<"Allele Count: "<<allCount<<"\n- Locus Count: "<<lociCount<<"\n";
	for (int i=0; i<freqs.size(); i++) {
		cout<<"* "<<freqs[i]<<"\n";
	}*/
	
	

	map<string, double> freqCache;

	for (size_t i=0; i<allCount; i++) {
//		cout<<"--"<<freqs[i*2]<<" : "<<freqs[(i*2)+1]<<"\n";
		freqCache[freqs[(i*2)]]=atof(freqs[(i*2)+1].c_str());
	}


	char al1[]="A";
	char al2[]="a";

	Reset();

	for (int i=0; i<lociCount; i++) {
		double f1=freqCache[al1];
		double f2=freqCache[al2];
		Freqs f(f1, f2);
		frequencies.push_back(f);
		al1[0]++;
		al2[0]++;
	}	

	wxGridTableMessage msg2(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, frequencies.size());
	GetView()->ProcessTableMessage(msg2);

}


void wxFreqGrid::SetValue( int row, int col, const wxString& value ) {
	Freqs &l = frequencies[row];
	double temp;
	switch (col) {
		case 0:
			if (value.ToDouble(&temp))	{
				l.al1 = temp;
				l.al2 = 1.0 - temp;
			}
			break;
		case 1:
			if (value.ToDouble(&temp)) {
				l.al1 = 1.0 - temp;
				l.al2 = temp;
			}
			break;
		default:
			cout<<"Invalid index for wxGridLoci: "<<col<<"\n";
			assert(0);
	}
	Refresh();					
}

wxString wxFreqGrid::GetRowLabelValue( int row ) {
	wxString rVal;
	rVal.Printf("Locus %d", row + 1);
	return rVal;
}

wxString wxFreqGrid::GetColLabelValue( int col ) {
	if (col == 0)
		return _T("All. 1");
	else 
		return _T("All. 2");
}

wxString wxFreqGrid::GetTypeName( int row, int col ) {
	return wxGRID_VALUE_FLOAT;
}


/*wxGridCellAttr *wxGridLoci::GetAttr(int row, int col) {
	return attr;
}
*/
bool wxFreqGrid::CanGetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)col >= frequencies.size())
		return false;

	if (typeName == wxGRID_VALUE_FLOAT)
		return true;
	else 
		return false;
}

bool wxFreqGrid::CanSetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)col >= frequencies.size())
		return CanGetValueAs(row, col, typeName);
	else
		return false;
}



double wxFreqGrid::GetValueAsDouble( int row, int col ) {
	if (col == 0) 
		return frequencies[row].al1;
	else
		return frequencies[row].al2;
}



void wxFreqGrid::SetValueAsDouble( int row, int col, double value ) {
	if (col == 0) {
		frequencies[row].al1 = value;
		frequencies[row].al2 = 1.0 - value;
	}
	else	{
		frequencies[row].al1 = 1.0 - value;
		frequencies[row].al2 = value;
	}
}

void wxFreqGrid::InvertFrequencies(int row) {
	row--;
	//cout<<"Inverting Frequencies for row: "<<row<<"\n";
	double temp = frequencies[row].al1;
	frequencies[row].al1 = frequencies[row].al2;
	frequencies[row].al2 = temp;
	Refresh();
}
void wxFreqGrid::Refresh() {
	//size_t rowCount = GetNumberRows();
	wxGridTableMessage msg2(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 0);//rowCount);// frequencies.size());
	GetView()->ProcessTableMessage(msg2);
}

size_t wxFreqGrid::GetLocusCount() {
	return lociCount;
}
void wxFreqGrid::SetLocusCount(size_t locusCount) {
	if (locusCount > frequencies.size())
		while (frequencies.size() < locusCount)
			frequencies.push_back(Freqs());

	//std::cout<<"Updating locus count: "<<locusCount<<" : "<<frequencies.size()<<"\n";

	if (GetView()) {
		int difference = locusCount - lociCount;
		if (difference > 0)	{
			wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, difference);
			GetView()->ProcessTableMessage(msg);
		}
		else if (difference < 0) {
			wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, locusCount, difference * -1);
			GetView()->ProcessTableMessage(msg);
		}
	}
	
	lociCount = locusCount;

}



wxPenGrid::wxPenGrid(int modelSize) : lociCount(0) { 
	SetLocusCount(modelSize);
}

int wxPenGrid::GetNumberCols() {
	return 1;
}

int wxPenGrid::GetNumberRows() { 
	int penSize = (int)pow((float)3, (float)lociCount);
	return penSize;
}


bool wxPenGrid::IsEmptyCell( int row, int col ) {
	return false;
}

wxString wxPenGrid::GetValue( int row, int col ) {	
	wxString rVal = wxT("");
	assert(col == 0);

	if ((size_t)row < penetrances.size()) {
		double value = penetrances[row];
	
//		cout<<"Fetching	"<<row<<"x"<<col<<" from SNP\n";
		rVal.Printf("%f",value);
	}
	return rVal;
}



void wxPenGrid::SetValue( int row, int col, const wxString& value ) {
	assert(col == 0);
	double temp;

	if (value.ToDouble(&temp))
		penetrances[row] = temp;
	Refresh();					
}

wxString wxPenGrid::GetRowLabelValue( int row ) {
	return _T(penLabels[row].c_str());

}

wxString wxPenGrid::GetTypeName( int row, int col ) {
	return wxGRID_VALUE_FLOAT;
}


/*wxGridCellAttr *wxGridLoci::GetAttr(int row, int col) {
	return attr;
}
*/
bool wxPenGrid::CanGetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)row >= penetrances.size())
		return false;

	if (typeName == wxGRID_VALUE_FLOAT)
		return true;
	else 
		return false;
}

bool wxPenGrid::CanSetValueAs( int row, int col, const wxString& typeName ) {
	if ((size_t)row >= penetrances.size())
		return CanGetValueAs(row, col, typeName);
	else
		return false;
}



double wxPenGrid::GetValueAsDouble( int row, int col ) {
	return penetrances[row];
}



void wxPenGrid::SetValueAsDouble( int row, int col, double value ) {
	penetrances[row] = value;
}
void wxPenGrid::Reset(bool clearPenetrances) {

	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_DELETED, 0, (int)pow((float)3.0, (float)lociCount));
		GetView()->ProcessTableMessage(msg);
	}
	if (clearPenetrances)
		penetrances.clear();

}


void wxPenGrid::Refresh() {
	wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, 0);//, penetrances.size());
	GetView()->ProcessTableMessage(msg);
}

void wxPenGrid::SetLocusCount(int locusCount) {
	BuildGenotypeLabels(locusCount);

	size_t penCount = (size_t)pow((float)3.0, (float)locusCount);

	Reset(false);
	while (penetrances.size() < penCount)
		penetrances.push_back(0.0);

	if (GetView()) {
		wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, penCount);
		GetView()->ProcessTableMessage(msg);
	}
	
	lociCount = locusCount;

}


string wxPenGrid::BuildGenotypeLabel(uint genotype, uint position) {
	char A='A'+position;
	char a='a'+position;
	
	char *label = new char[3];
	if (genotype == 0)
		sprintf(label, "%c%c", A, A);
	else if (genotype == 1)
		sprintf(label, "%c%c", A, a);
	else if (genotype == 2)
		sprintf(label, "%c%c", a, a);
	string finalLabel = label;
	delete[] label;
	return finalLabel;
}


string wxPenGrid::BuildGenotypeLabel(uint *genotypes, uint modelSize) {
	stringstream ss;
	for (uint i=0; i<modelSize; i++) {
		ss<<BuildGenotypeLabel(genotypes[i], i);
		//cout<<genotypes[i]<<" ";
	}
	//cout<<"= "<<ss.str()<<"\n";

	return ss.str();
}


void wxPenGrid::BuildGenotypeLabels(uint modelSize) {
	uint *genotypes = new uint[modelSize+1];
	uint position = 0;
	penLabels.clear();

//	for (uint i=0; i<modelSize; i++) 
//		genotypes[i]=0;
	memset((void*)genotypes, 0, (modelSize + 1 )*sizeof(uint));

	penLabels.clear();	
	position = modelSize - 1;

	while (genotypes[0]<3) {	
		string label = BuildGenotypeLabel(genotypes, modelSize);
		
		penLabels.push_back(label);
		//cout<<id++<<"\t"<<label<<"\n";

		if (++genotypes[position]>2 && position > 0) {
			//Find the highest position of rollover
			while (position-- > 0 && ++genotypes[position] > 2) {}

			while (position < modelSize - 1) 
				genotypes[++position] = 0;
		}
	}	
}




/*!
 * wxEVT_GRID_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE
 */

void wxPagePenetranceModel::OnLabelRightClick( wxGridEvent& event )
{
////@begin wxEVT_GRID_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE in wxPagePenetranceModel.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE in wxPagePenetranceModel. 
}


/*!
 * wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_PENETRANCE
 */

void wxPagePenetranceModel::OnGrdPenetranceLabelLeftClick( wxGridEvent& event )
{
////@begin wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_PENETRANCE in wxPagePenetranceModel.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_PENETRANCE in wxPagePenetranceModel. 
}


bool wxPagePenetranceModel::Import(const char *filename, string &comment) {
	FileToMap file;
	bool success = file.Parse(filename);
	
	if (success) {
		comment = file.comments;

//		cout<<"Original Comment: "<<comment<<"\n";
		if (comment.length() > 0) {
			size_t loc = comment.rfind("\n#", comment.length() - 1);
				
			while (loc != string::npos) {
//				cout<<"Found at: "<<loc<<"\n";
				comment = comment.erase(loc + 1, 1);
				loc = comment.rfind("\n#", loc - 1);	
			}
		}		

		if (comment[0] == '#')
			comment=comment.erase(0,1);
		//cout<<"Loaded Comment:   "<<comment<<"\n";

		vector<string> values;

		txtFreqThreshold->SetValue(_(file.GetLine("FREQ_THRESHOLD").c_str()));

		vector<string> lines;
		if (file.GetLines("FREQ", lines)) {
			freqGridMaster->Load(lines);
			//freqGridMaster->Refresh();
		}
		penGridMaster->Load(freqGridMaster->GetLocusCount(), file);
		//penGridMaster->Refresh();
		
		dialLocusCount->SetValue(freqGridMaster->GetLocusCount());
	
		
	}

	return success;
}

bool wxPagePenetranceModel::Save(const char *filename, const char *description) { 
	bool success=true;
		
	string desc(description);

	if (desc.length() > 0) {
		size_t loc = desc.rfind('\n', desc.length() - 1);
	
		while (loc != string::npos) {
			desc.insert(loc + 1, 1, '#');
			loc = desc.rfind('\n', loc - 1);	
		}
	}
	
	cout<<"Saving "<<filename<<" Desc: "<<description<<"\n";
	ofstream file(filename);
	file<<"#"<<desc<<"\n\n";
	file<<"# The Threshold will help protect genomeSIMLA from using a penetrance\n";
	file<<"# table with data whose allele frequencies vary too greatly from the population\n";
	file<<"# for which the table was designed.\n";
	file<<"FREQ_THRESHOLD "<<ExtractDouble(txtFreqThreshold)<<"\n";
	file<<"\n\n#Allele Frequencies Associated with this Penetrance Table\n";

	char al1 = 'A';
	char al2 = 'a';
	int locCount = dialLocusCount->GetValue();
	for (int i=0; i<locCount; i++){
		file<<"FREQ "<<al1++<<" "<<freqGridMaster->GetValue(i, 0)<<"\n";
		file<<"FREQ "<<al2++<<" "<<freqGridMaster->GetValue(i, 1)<<"\n";
	}
	
	file<<"\n\n# The Penetrances are listed below\n";
	int gtCount = penGridMaster->GetNumberRows();
	for (int i=0; i<gtCount; i++) {
		file<<penGridMaster->GetRowLabelValue(i)<<" "<<penGridMaster->GetValue(i, 0)<<"\n";
	}

	file.close();
		
	
	return success;
}


/*!
 * wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_PENETRANCE
 */

void wxPagePenetranceModel::OnGrdPenetranceLabelRightClick( wxGridEvent& event )
{
	event.Skip();
	//ShowContextMenu(pos, row);
}


void wxPagePenetranceModel::ShowContextMenu(const wxPoint& pos, int SNP) {
	wxMenu menu;
	selectedSnp = SNP;
	menu.Append(ID_INVERT_ALLELES, _("Invert Alleles"), _("Swap allele frequencies in this table and in the penetrance table"), wxITEM_NORMAL);
	PopupMenu(&menu, pos.x, pos.y);

}


bool wxPagePenetranceModel::Evaluate() {
	bool success = true;

	stringstream ss;
	PenetranceEval eval;
	int locusCount = freqGridMaster->GetLocusCount();
	
	char al1[]="A";
	char al2[]="a";

	for (int i=0; i<locusCount; i++) {
		eval.AddLocus(freqGridMaster->GetValueAsDouble(i, 0), freqGridMaster->GetValueAsDouble(i, 1), al1[0], al2[0]);
		al1[0]++;
		al2[0]++;
	}
	eval.InitializeCells();
	int cellCount = penGridMaster->GetNumberRows();
	
	for (int i=0; i<cellCount; i++) {
		string cellID = penGridMaster->GetRowLabelValue(i).c_str();
		eval.GetCell(cellID.c_str()).penetrance = penGridMaster->GetValueAsDouble(i,0);
	}

	eval.Evaluate(ss);
	wxDlgBasicHtmlReport rpt(this);
	rpt.Append(ss.str().c_str());
	//cout<<ss.str()<<"\n";
	rpt.ShowModal();
	return success;
}


void wxPagePenetranceModel::InvertAlleles(wxCommandEvent& event) {
	cout<<"Inverting SNP "<<selectedSnp<<"\n";
	penGridMaster->InvertAllele(selectedSnp);
	
	//Swap the frequencies
	freqGridMaster->InvertFrequencies(selectedSnp);
 	penGridMaster->Refresh();
}


/*!
 * wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_ALLELE_FREQ
 */

void wxPagePenetranceModel::OnGrdAlleleFreqLabelLeftClick( wxGridEvent& event )
{
////@begin wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_ALLELE_FREQ in wxPagePenetranceModel.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_CMD_LABEL_LEFT_CLICK event handler for ID_GRD_ALLELE_FREQ in wxPagePenetranceModel. 
}


/*!
 * wxEVT_GRID_CMD_LABEL_RIGHT_CLICK event handler for ID_GRD_ALLELE_FREQ
 */
void wxPagePenetranceModel::OnGrdAlleleFreqLabelRightClick( wxGridEvent& event )
{
	int row = event.GetRow();

	wxPoint point = event.GetPosition();

	point = ScreenToClient(::wxGetMousePosition());
	if (row >= 0) 
		ShowContextMenu(point, row + 1);

}

/*!
 * wxEVT_SIZE event handler for ID_GRD_ALLELE_FREQ
 */

void wxPagePenetranceModel::OnSize( wxSizeEvent& event )
{
	int labelWidth = grdAlleleFrequencies->GetColLabelSize();
	int availableWidth = event.GetSize().GetWidth() - labelWidth - 60;
	
	grdAlleleFrequencies->SetColSize(0, availableWidth/2);
	grdAlleleFrequencies->SetColSize(1, availableWidth/2);

	event.Skip();
	
}


bool wxPenGrid::InvertAllele(int locusToInvert) {
	vector<double> tableCopy(penetrances);

	int totalGtCount = (int)pow((float)3.0, (float)lociCount);
	vector<int> genotypes(lociCount);

	int idx = locusToInvert-1;
	for (int i=0; i<totalGtCount; i++) {
		ConvertToGenotypes(genotypes, lociCount, i);
		if (genotypes[idx] == 0) {
			genotypes[idx] = 2;
		}
		else if (genotypes[idx] == 2) {
			genotypes[idx] = 0;
		}
	 
		int destGt;
		ConvertFromGenotypes(genotypes, lociCount, destGt);
		penetrances[destGt] = tableCopy[i];
	}
	return true;
}

void wxPenGrid::ConvertFromGenotypes(vector<int>& genotypes, int genotypeCount, int& mlGenotype) {
	int currOffset = 1;
	
	mlGenotype = 0;
	for (int i=0; i<genotypeCount; i++) {
		mlGenotype+=(currOffset*genotypes[i]);
		currOffset*=3;
	}
}

void wxPenGrid::ConvertToGenotypes(vector<int>& genotypes, int genotypeCount, int mlGenotype) {
	//mlGenotype++;
	genotypes.clear();
	genotypes.resize(genotypeCount);
	int curOffset = (int)pow(3.0, (int)genotypeCount-1);
	int remainderGT = mlGenotype;
	for (int i=genotypeCount - 1; i>=0; i--) {
		int gt = remainderGT/curOffset;
		remainderGT = remainderGT%curOffset;
		//genotypes.push_back(gt);
		genotypes[i]=gt;
		curOffset /= 3;
	}	
}





/*
void wxPagePenetranceModel::Save(const char *filename) {
	ofstream file(filename);

	file<<"# genomeSIMLA Penetrance filename\n";
	file<<"#\n# Created using wxGenomeSIMLAs\n";
	file<<"-"<<	
}
*/


/*
	vector<double> penetrances;
	int lociCount;
	vector<string> penLabels;
*/

bool wxPenGrid::Load(int lociCount, FileToMap &file) {	
	Reset(false);	
	
	SetLocusCount(lociCount);
	
	size_t gtCount = GetNumberRows();


	for (size_t i=0; i<gtCount; i++) {
		string key = GetRowLabelValue(i).c_str();
		string line = file.GetLine(GetRowLabelValue(i).c_str());
		penetrances[i]=atof(line.c_str());
	}

	//Refresh();

	wxGridTableMessage msg(this, wxGRIDTABLE_NOTIFY_ROWS_APPENDED, penetrances.size());
	GetView()->ProcessTableMessage(msg);

	return true;
	
	
}

bool	wxPenGrid::Save(const char *file) {
	//We aren't doing anything here, which seems wrong
	assert(0);
	return true;
	
}




/*!
 * wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_ALLELE_FREQ
 */

void wxPagePenetranceModel::OnEditorShown( wxGridEvent& event )
{
////@begin wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_ALLELE_FREQ in wxPagePenetranceModel.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_GRID_EDITOR_SHOWN event handler for ID_GRD_ALLELE_FREQ in wxPagePenetranceModel. 
}


}

}



