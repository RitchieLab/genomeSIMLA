//
// C++ Interface: locusmanager
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLOCUSMANAGER_H
#define SIMULATIONLOCUSMANAGER_H
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "locus.h"
#include "utility/exception.h"
#include "utility/rbtree.h"
#include "utility/stat.h"

#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif
namespace Simulation {

#ifdef CPPUNIT
void WriteLocusFile(const char *filename);
void WriteLocusXY(const char *filename);
void WriteLocusFile(const char *filename, int locusCount);
#endif


using namespace std;




template <class T>
struct LocusLT { 
	bool operator()(const T *l, const T *r) const {
		return l->GetLocation() < r->GetLocation();
	}

};

struct floatLT { 
	int operator()(const float l, const float r) const {
		if (l<r) return -1;
		if (l>r) return 1;
		return 0;
	}
};

/**
@Brief Represents a single chromosome's loci. 

	@author Eric Torstenson
*/
template <class T>
class LocusManager {
public:
typedef Utility::RBTree<float, T *, floatLT> PositionLookup;
typedef Utility::RBTreeNode<float, T *, floatLT> PositionLookupNode;
	/**
	 * @Brief basc construction: project is used to create filenames
	 * @param chromID Integer ID (generally an index into an array of chromosomes)
	 * @param lambdaDev This is the deviation of the poison lambda from 1 for males and females. This is the value added or subtracted from 1.0 depending on male/females and the resulting value is used to adjust XO in males or females. XY chromosomes should have a lambdaDev or 0.0 since these two chromosomes should have already been adjusted for gender XO Rates
	 */
    LocusManager(int chromID, float lamdbaDev=0.25) : isSuspended(false), minThreshold(0.0), chromID(chromID), lambdaDev(lambdaDev), gdStart(0.0), gdLength(0.0) { }						
    ~LocusManager();

	int LocusCount();

	/**
	 * @Brief Returns locus at the specified index (sorted by position)
	 * @Note Index is zero based
	 */
	T *operator[](int idx);

	/**
	 * @Brief Returns locus at the specified index (sorted by position)
	 * @Note Index is zero based
	 */
	T *At(int idx);
	/**
	 * @Brief Returns locus with specified lable (if it exists)
	 */
	T *operator[](const char *label);

	T *At(const char *label);

	int ChromID();
	void ChromID(int id);
	bool BuildXOEventList(Utility::Random& rnd, std::vector<size_t>& events, bool isXX);
	/**
	 * @Brief Return the locus based on the genetic position on the chromosome
	 * @Note The child class is responsible for building this up
	 */
	T *operator[](float geneticPosition);
	T *At(float geneticPosition);

	void Show(float geneticPosition);

	/**
	 * @Brief Build the array of loci usable for LD calculations
	 * @param start The first SNP to be added to the array
	 * @param stop The last SNP
	 * @param loci the locus array to be populated
	 */
	virtual int GetLocusArray(int start, int stop, vector<Locus>& loci, int type);

	virtual void BuildLocusMap(map<string, Locus>& locusMap);

	void SetPoissonLambda(float lambda);

	
	void ForceAlleleFrequency(uint locusID, float al1, float al2);
protected:
	vector<T*> sortedLoci;				///<This is just a linear version of the map, sorted by the <
	std::map<std::string, T*> loci;		///<The main "database" of loci for this "array"
	bool isSuspended;					///<Indicate that the data needs to be "resourced"
	float minThreshold;
	int chromID;						///<ID associated with chromosome index
	PositionLookup recombIndexLookup;	///<Find a locus based on it's genetic position
	float poissonLambda;				///<Lambda associated with this report
	float lambdaDev;					///<This number will be added/subtracted to 1.0 to adjust the poisson lambda based on the gender of the individual in which the XO is occurring	
	float gdStart;						///<The genetic distance of the first snp from the beginning of the chromosome
	float gdLength;						///<genetic distance between first and last SNP on the chromosome
};

template <class T>
inline
void LocusManager<T>::SetPoissonLambda(float lambda) {
	poissonLambda = lambda;
}
template <class T>
inline
int LocusManager<T>::GetLocusArray(int start, int stop, vector<Locus>& loci, int t) {
	assert (start < stop);
	assert (stop <= LocusCount());
	int fixedCount = 0;
	for (int i=start; i<stop; i++) {
		Locus l;
		if (sortedLoci.at(i)->Clone(l, t, fixedCount)) {
			loci.push_back(l);
		}
	}
	return fixedCount;
}

template <class T>
inline
void LocusManager<T>::BuildLocusMap(map<string, Locus>& locusMap) {
	for (int i=0; i<sortedLoci.size(); i++) {
		Locus l;
		if (sortedLoci.at(i)->Distill(l))
			locusMap[l.GetLabel()] = l;
	}
}


template <class T>
inline
bool LocusManager<T>::BuildXOEventList(Utility::Random& rnd, std::vector<size_t>& events, bool isFemale) {
	float lambdaAdjustment = 1.0 - lambdaDev;
	if (isFemale)
		lambdaAdjustment = 1.0 + lambdaDev;

	float lambda = poissonLambda * lambdaAdjustment;
	size_t eventCount = PoissonEventCount(rnd, lambda);

	for (size_t i=0; i<eventCount; i++) {
		float loc = gdStart +  rnd(gdLength);
		Locus *locus = At(loc);
		int idx = locus->GetID();
		vector<size_t>::iterator other = find(events.begin(), events.end(), idx);
		if (other != events.end())
			events.erase(other);
		events.push_back(idx);
	}
	sort(events.begin(), events.end());
	return events.size()%2;
}

template <class T>
class LocusManagerFileBased  : public LocusManager<T> {
public:
	LocusManagerFileBased(const char *project, int chromID);
	~LocusManagerFileBased();
	

	/**
	 * @Brief verify the array object.
	 * @Param If generation is -1, we will call the virtual version of the function, otherwise, we are checking for the present of filenames
	 * @note Returns appropriate error Message
	 */
	void Verify(int generation=-1);

	/**
	 * @brief Saves the current state to file 
	 */
	void Save(int generation, const char *label);
	/**
	 * @Brief Loads the specified filename. 
	 * @note This should only be used by the client during initialization-subsequent subsequent calls should be made to refresh
	 */
	void Load(const char *filename);

	/**
	 * @Brief Refreshes locus details, like allele frequency, based on a pre-existing file which is named based on the local "project" name and the specified generation
	 */
	void Refresh(int generation);
	/**
	 * @brief Write locus details to file
	 * @param label Used inside locus report for informational purposes....probably should be the chromosome name
	 * @param filename The name of the file to be written
	 */
	void WriteLocusReport(const char *label, const char *filename);
	void WriteLoci(ostream& os);

	void WriteMarkerInfo(const char *filename);
	/**
	 * @brief constructs filename based on pattern and generation
	 */
	std::string GetFilename(int generation);

	void SetFilename(const char *filename);

	void SetPrefix(const char *prefix);

	std::string Filename();
protected:
	std::string sourceFilename;			///<This is recorded mainly for reporting purposes
	std::string prefix;					///<Used to construct the filename
};

template<class T>
inline
void LocusManagerFileBased<T>::SetPrefix(const char *prefix) {
	this->prefix = prefix;
}


template<class T>
inline
void LocusManager<T>::ForceAlleleFrequency(uint locusID, float al1, float al2) {
	assert(locusID < sortedLoci.size());
	T *locus = sortedLoci[locusID];
	locus->SetFreq1(al1);
	locus->SetFreq2(al2);
}

template<class T>
inline
int LocusManager<T>::ChromID() {
	return chromID;
}

template<class T>
inline
void LocusManager<T>::ChromID(int id) {
	chromID = id;
}


template <class T>
inline
T *LocusManager<T>::operator[](float geneticPosition) {
	return At(geneticPosition);
}

template <class T>
inline
T *LocusManager<T>::At(float geneticPosition) {
	T *locus = NULL;
	PositionLookupNode *node = recombIndexLookup.FindNearestMin(geneticPosition);
	if (node)
		locus = node->GetData();
	return locus;
}

template <class T>
inline
void LocusManager<T>::Show(float geneticPosition) {
	T *locus = NULL;
	int idx = 0;
	PositionLookupNode *node = recombIndexLookup.GetFirst();
	while (node) {
		if (geneticPosition == node->GetData()->GetID())
			cerr<<"*";
		cerr<<"\t"<<++idx<<"\t"<<node->GetKey()<<"\t"<<node->GetData()->GetID()<<"\n";
		node = node->GetNext();
	}
}


template<class T>
inline
T *LocusManager<T>::operator[](const char *label) {
	return At(label);
}

template<class T>
inline
T* LocusManager<T>::At(const char *label) {
	T *loc = NULL;
	typename map<string, T* >::iterator itr = loci.find(label);
	if (itr != loci.end() )
		loc = itr->second;
	return loc;	
}

template<class T>
inline
T *LocusManager<T>::At(int idx) {
	assert(idx < loci.size());
	return sortedLoci[idx];
}

template<class T>
inline
T *LocusManager<T>::operator[](int idx) {
	return At(idx);
}

template<class T>
inline
LocusManager<T>::~LocusManager() {
	typename std::map<std::string, T*>::iterator itr = loci.begin();
	typename std::map<std::string, T*>::iterator end = loci.end();
	
	while (itr != end) {
		delete (itr++)->second;
	}
}

template<class T>
inline
int LocusManager<T>::LocusCount() {
	return loci.size();
}










/******************************************************************************
 ***             Locus Array File Based
 ******************************************************************************/
template<class T>
inline
LocusManagerFileBased<T>::LocusManagerFileBased(const char *project, int chromID) : LocusManager<T>(chromID), prefix(project) { }

template<class T>
inline
LocusManagerFileBased<T>::~LocusManagerFileBased() { }

template<class T>
inline
void LocusManagerFileBased<T>::SetFilename(const char *filename) {
	this->sourceFilename = filename;
}

template<class T>
inline
string LocusManagerFileBased<T>::Filename() {
	return sourceFilename;
}

template<class T>
inline
void LocusManagerFileBased<T>::WriteLoci(ostream& os) {
	int lociCount = LocusManager<T>::loci.size();
	if (lociCount > 0) {
		if (LocusManager<T>::sortedLoci[0]->GetMinAlleleFreq() >= LocusManager<T>::minThreshold)
		//os << setw(width)<<"Locus ID" << " ";
			LocusManager<T>::sortedLoci[0]->WriteHeader(os, 20);
	}
	for(uint i=0; i<lociCount; i++){
		//os << setw(width) << i+1 << " ";
		if (LocusManager<T>::sortedLoci[i]->GetMinAlleleFreq() >= LocusManager<T>::minThreshold)
			LocusManager<T>::sortedLoci[i]->WriteFormatted(os, 20);
	}
}


template <class T>
inline
void LocusManagerFileBased<T>::WriteMarkerInfo(const char *filename) {
	ofstream os(filename);
	if (os.is_open()) {
		int width=16;
		int lociCount = LocusManager<T>::loci.size();	
		for(uint i=0; i<lociCount; i++){
			//os << setw(width) << i+1 << " ";
			if (LocusManager<T>::sortedLoci[i]->GetMinAlleleFreq() >= LocusManager<T>::minThreshold)
				LocusManager<T>::sortedLoci[i]->WriteMarkerInfo(os);
		}
	}
	else {
		cout<<"FileIOError: "<<filename<<"\n";
		throw Utility::Exception::FileNotWritable(filename);
	}
}
template<class T>
inline
void LocusManagerFileBased<T>::WriteLocusReport(const char *label, const char *filename) {
	ofstream os(filename);
	if (os.is_open()) {
		int width=16;
		int lociCount = LocusManager<T>::loci.size();
		os << "Locus Log for "<<label<<"\n";
		os << lociCount<<" Loci\n";
	
		WriteLoci(os);
	}
	else {
		cout<<"FileIOError: "<<filename<<"\n";
		throw Utility::Exception::FileNotWritable(filename);
	}
}

template<class T>
inline
void LocusManagerFileBased<T>::Save(int generation, const char *label) {
	string filename = GetFilename(generation);
	WriteLocusReport(label, filename.c_str());
}

template<class T>
inline
void LocusManagerFileBased<T>::Refresh(int generation) {	
	//Attempting to "reload" a locus file is a waste of time. If we want to load a new pool, we need to make sure to get rid of the old stuff
	string filename = GetFilename(generation);
	ifstream file(filename.c_str(), ios_base::in);
	if (file.is_open()) {
		char line[4096];
		file.getline(line, 4096);			///<The chromosome id
		file.getline(line, 4096);			///<The count of loci
		file.getline(line, 4096);			///<The header information
		int lociCount = LocusManager<T>::loci.size();
		float geneticPosition = 0.0;
		while (!file.eof()) {
			T *loc = new T(LocusManager<T>::chromID, lociCount);

			file>>*loc;
			loc->MapPosition(geneticPosition);
			geneticPosition+=loc->MapDistance();

			if (loc->Valid()) {
				if (LocusManager<T>::loci.find(loc->GetLabel())==LocusManager<T>::loci.end()) {
					delete loc;
					string errMsg = "The locus file, " + filename + ", had one or more loci that was missing from the previous definition (" + loc->GetLabel() + ").";
					throw Utility::Exception::General(errMsg.c_str());
				}
				Locus *oldLocus = LocusManager<T>::loci[loc->GetLabel()];
				oldLocus->SetFreq1(loc->Freq1());
				oldLocus->SetFreq2(loc->Freq2());
				//We just want to update the contents of the locus-not replace any pointers 
				//*(LocusManager<T>::loci[loc->GetLabel()]) = *loc;
				
			}
			delete loc;
		}
	}
	else {
		throw Utility::Exception::FileNotFound(filename.c_str());
	}
}

template<class T>
inline
void LocusManagerFileBased<T>::Load(const char *filename) {
	assert(LocusManager<T>::loci.size()==0);
	sourceFilename = filename;
	ifstream file(filename, ios_base::in);
	if (file.is_open()) {
		char line[4096];
		file.getline(line, 4096);			///<The chromosome id
		file.getline(line, 4096);			///<The count of loci
		file.getline(line, 4096);			///<The header information
	
		//We need to make sure we are starting with an empty locus pool
		LocusManager<T>::loci.clear();
		while (!file.eof()) {
			int lociCount = LocusManager<T>::loci.size();
			T *loc = new T(LocusManager<T>::chromID, lociCount + 1);

			file>>*loc;

			if (loc->Valid()) {
				if (LocusManager<T>::loci.find(loc->GetLabel())!=LocusManager<T>::loci.end()) {
					string errMsg = "The locus file, " + string(filename) + ", had one or more duplicate loci (" + loc->GetLabel() + ").";
					delete loc;
					throw Utility::Exception::General(errMsg.c_str());
				}
				//we want to store these pointers for others to have access to 
				LocusManager<T>::loci[loc->GetLabel()] = loc;
				LocusManager<T>::sortedLoci.push_back(loc);
			}
			else {
				delete loc;
			}
		}

		//Sort the loci based on position- We can't sort by genetic position, due to the fact that there might not be any genetic positions yet (Calibration should take care of that)
		sort(LocusManager<T>::sortedLoci.begin(), LocusManager<T>::sortedLoci.end(), LocusLT<T>());
		float mapDistance = 0.0;
		typename T::PositionCalibration calibration;
		typename vector<T*>::iterator itr = LocusManager<T>::sortedLoci.begin();
		typename vector<T*>::iterator end = LocusManager<T>::sortedLoci.end();
		while (itr != end) {
			T *locus = *itr;
			calibration.Calibrate(*locus);
			LocusManager<T>::recombIndexLookup.Set(locus->MapPosition(), locus);
			itr++;
		}
	}
	else 
		throw Utility::Exception::FileNotFound(filename);
	if (LocusManager<T>::sortedLoci.size() == 0) {
		cerr<<"Locus File, "<<filename<<", seemed to be empty. This chromosome will have no SNPs in it\n";
		return;
	}
	LocusManager<T>::gdStart = LocusManager<T>::sortedLoci[0]->MapPosition();
	LocusManager<T>::gdLength = LocusManager<T>::sortedLoci[LocusManager<T>::sortedLoci.size()-1]->MapPosition() - LocusManager<T>::gdStart;
	LocusManager<T>::poissonLambda = LocusManager<T>::gdLength * 0.01;
}

template<class T>
inline
std::string LocusManagerFileBased<T>::GetFilename(int generation) {
	stringstream filename;
	if (generation == -1)
		filename << this->sourceFilename;
	else
		filename<<prefix<<"."<<generation<<".loc";
	return filename.str();
}

template<class T>
inline
void LocusManagerFileBased<T>::Verify(int generation) {
	string filename = GetFilename(generation);
	//Check that the file exists
	ifstream file(filename.c_str);
	if (!file.is_open()) {
		string errMsg = "An error was encountered opening the locus file, " + filename +".";
		throw Utility::Exception::General(errMsg.c_str());
	}
}






#ifdef CPPUNIT
class LocusManagerTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( LocusManagerTest );
	CPPUNIT_TEST( TestRegularLocus );
	CPPUNIT_TEST( TestLocusXY );
	CPPUNIT_TEST_SUITE_END();
public:
	LocusManagerTest();
	~LocusManagerTest();

	void setUp();
	void tearDown();


	void TestRegularLocus();
	void TestLocusXY();
};

#endif //CPPUNIT





}

#endif
