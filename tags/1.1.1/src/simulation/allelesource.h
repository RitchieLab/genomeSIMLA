//
// C++ Interface: allelesource
//
// Description: 
//
//
// Author:  <Eric Torstenson>, Marylyn Ritchie (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONALLELESOURCE_H
#define SIMULATIONALLELESOURCE_H
#include "par_region.h"
#include <iostream>
#include "locusmanager.h"
#include <boost/dynamic_bitset.hpp>
#include "utility/rbtree.h"
#include "kinshipcalculator.h"
#include <iomanip>



#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif
namespace Simulation {

extern bool DoSuspend;			///<Determines if wake/suspend actually need to conserve memory

template <class T> class AlleleSource;
/**
 * @brief streaming operator...requires pointer, since we are assuming some magic occurs with polymorphism within
 * @note This is only useful for doing text output in phased format
*/
template <class T>
inline
std::ostream& operator<<(std::ostream & os, AlleleSource<T>* source) {
	for(unsigned int i=0; i<source->locusCount; i++)
    	os << source->At(1)+1 << " ";
	os<<"\n";
  	return os;
}
/**
@Brief Functional details associated with various ways of representing allelic information inside the chromosome

	@author Eric Torstenson
*/
template<class T>
class AlleleSource{
public:

   AlleleSource(LocusManager<T>* loci, int sourceID) : loci(loci), sourceID(sourceID)  {
		locusCount = loci->LocusCount();
	}
    virtual ~AlleleSource() {}

	/**
	 * @Brief depending on the type this is, the output could be anything from 0-1 to 0-N
	 * @note We can't use the operator[] because some forms might not have a single
	 * 		source for the data to return a reference to. As a result, we need two
	 * 		distint functions for the different behaviors.
	 */
	virtual int At(uint locusIndex)=0;
	virtual void At(uint locusIndex, int val)=0;
	/**
	 * @brief Performs cross-over with "other"
	 * @Note Unlike the previous chromosome, we are assuming this is the phased chromosome
	 * @note If the user wants this to be phased, but start on the other chromosome for 
	 * 			efficiency reasons, the first XO point should be 0
	 */
	AlleleSource<T>* Cross(AlleleSource<T>* other, vector<size_t>& xoEvents);
	/**
	 * @brief Crosses 2 alleles and returns a new chromosome
	 * @param father The data passed from father. If X/Y, this is assumed to be a Y
	 * @param par Only used for X/Y cross-overs. If is Autosomal, then use NULL
	 * @param returnX only used when par!=NULL. This indicates the return type is X or Y
	 * @param rnd the random number generator
	 * @note if T is LocusXY, returnX is used to determine X/Y of the returned chromosome only if
	 * 			par is non-null (pass NULL if crossing 2 Xs)
	 */
	AlleleSource<T>* Cross(AlleleSource<T>* father, Utility::Random& rnd, PAR_Region<T> *par=NULL, bool returnX=true);	

	virtual int GetMaxGenotypeCount( )=0;

	/**
	 * @Brief write to alleles to binary pool
	 */
	virtual void WriteBinary(std::ostream& file)=0;
	virtual void WriteCompleteBinary(std::ostream& file, int sourceID)=0;
	virtual void PrintCrossOvers(std::ostream& os)=0;
	virtual int CountCrossOvers()=0;
	/**
	 * @brief read alleles from binary
	 */
	virtual void ReadBinary(std::istream& file, Utility::RBTree<size_t, AlleleSource<T>*>* sourceAlleles)=0;

	/**
	 * @brief Resets the internal representation to hold ONLY the model loci
	 */
	virtual void ResetToMinimal(vector<uint>& modelLoci) { }

	virtual size_t GetLocusCount();		///<Count the loci associated with this object
	
	/**
	 * @brief Returns a vector containing the allele sources for a region of an object
	 * @note The lower bound should occur at all times unless lowerBound==upperBound. 
	 * 		The upper boud should NEVER be inserted
	 *		These indices are the point of switching, so the local chromosome would not contribute data
	 * 		at a given locus, but switch phase at that point. The first entry at 5 would indicate that
	 * 		the primary phased source contributes alleles 0-4.
	 */
	virtual void GetSource(Utility::RBTree<size_t, AlleleSource<T>*>&, int lowerBound, int upperBound);

	LocusManager<T> *GetLoci() { return loci; }

	int GetSourceID() { return sourceID; }
	void SetSourceID(int sourceID) { this->sourceID=sourceID; }
	/**
	 * @brief Clone the local object using polymorphism
	 */
	virtual AlleleSource *Clone()=0;

	/**
	 * @brief Attempt to reconstruct a complete chromosome based on local data
	 */
	virtual AlleleSource *Realize(int sourceID) = 0;

	virtual void Suspend(vector<uint>& modelLoci) = 0;
	virtual void Wake(istream& infile) = 0;

	/**
	 * @brief Attempt to reconstruct a complete chromosome based on other
	 */
	virtual AlleleSource *Realize(int sourceID, AlleleSource*other) = 0;

	/**
	 * @brief Resets the information rendering a blank slate
	 */
	virtual void Reset()=0;

	virtual bool HasCompleteGenotypes()=0;

	virtual int EvaluateKinship(AlleleSource<T> *other);

	/**
	 * @brief streaming operator...requires pointer, since we are assuming some magic occurs with polymorphism within
	 * @note This is only useful for doing text output in phased format
	 */	
	friend std::ostream& operator<<<>(std::ostream & os, AlleleSource<T>* source);
	virtual void GetSourceData(Utility::RBTree<size_t, AlleleSource<T>*>& srcs)=0;
	void ShowGenotypes(std::ostream& os, int count, char div);
	
	/**	
	 * @brief this is intended for debugging, not really comparing in an efficient way
	 */
	bool Compare(AlleleSource<T>* other, int idx);
protected:
	int locusCount;						///<number of loci
	LocusManager<T> *loci;				///<Shared locus array (don't ever delete here) 
	int sourceID;						///<ID for use in descendants
};



/** 
 * @brief bi-allelic SNPs. 
 */
template <class T>
class AlleleSourceSingle : public AlleleSource<T> {
public:
    AlleleSourceSingle(LocusManager<T>* loci, int sourceID);
	AlleleSourceSingle(const AlleleSourceSingle& other);
	~AlleleSourceSingle();

	int At(uint locusIndex);			///<Returns value at given locus
	void At(uint locusIndex, int val);	///<Sets value at given locus


//	int EvaluateKinship(AlleleSource<T>* other);
	/**
	 * @brief Performs cross-over with "other"
	 * @Note Unlike the previous chromosome, we are assuming this is the phased chromosome
	 * @note If the user wants this to be phased, but start on the other chromosome for efficiency reasons,
	 * 			the first XO point should be 0
	 */
//	AlleleSource<T>* Cross(AlleleSource<T>* other, vector<size_t>& xoEvents);
	/**
	 * @brief Crosses 2 alleles and returns a new chromosome
	 * @param father The data passed from father. If X/Y, this is assumed to be a Y
	 * @param par Only used for X/Y cross-overs. If is Autosomal, then use NULL
	 * @param returnX only used when par!=NULL. This indicates the return type is X or Y
	 * @param rnd the random number generator
	 * @note if T is LocusXY, returnX is used to determine X/Y of the returned chromosome only if
	 * 			par is non-null (pass NULL if crossing 2 Xs)
	 */
//	AlleleSource<T>* Cross(AlleleSource<T>* father, Utility::Random& rnd, PAR_Region<T> *par=NULL, bool returnX=true);

	/**
	 * @Brief write to alleles to binary pool using the most "compressed" manner
	 */
	void WriteBinary(std::ostream& file);

	/**
	 * @Brief Read data from the pool
	 * @param file Source being read
	 * @param sourceAlleles Means nothing for this object
	 */
	void ReadBinary(std::istream& file, Utility::RBTree<size_t, AlleleSource<T>*>* sourceAlleles);

	/**	
	 * @brief Write binary data for each locus
	 */
	void WriteCompleteBinary(std::ostream& file, int sourceID);

	void PrintCrossOvers(std::ostream& os);

	int CountCrossOvers();
	/**
	 * @brief Resets the internal representation to hold ONLY the model loci
	 */
	void ResetToMinimal(vector<uint>& modelLoci);

	/**
	 * @brief returns the maximum number of genotypes possible
	 */
	int GetMaxGenotypeCount();
	/**
	 * @brief Produce a copy
	 */
	AlleleSource<T> *Clone();

	void Wake(istream& file);
	void Suspend(vector<uint>& modelLoci);
	/**
	 * @brief Attempt to reconstruct a complete chromosome based on local data
	 */
	AlleleSource<T> *Realize(int sourceID) { AlleleSource<T> *newSource = new AlleleSourceSingle(*this); return newSource; }
	/**
	 * @brief Attempt to reconstruct a complete chromosome based on other
	 */
	AlleleSource<T> *Realize(int sourceID, AlleleSource<T>*other);
	/**
	 * @brief Return object to clean state
	 */
	void Reset();
	
	bool HasCompleteGenotypes() { return alleles!=NULL; }

	void GetSourceData(Utility::RBTree<size_t, AlleleSource<T>*>& srcs);
	void InitLoci(Utility::Random& rnd, bool isY);
protected:
	boost::dynamic_bitset<> *alleles;		///<real data
	map<size_t, int> modelOnlyLoci;		///<"Minimal" binary support
};



/**
 * @brief specialized version of the source data, in which the data lies in 2 or more sources.
 * @note This type of source is the result of a mating of two chromosomes
 * Compressed Binary format:
 * 	This is a special format where the source stores ONLY the XO points and the source allele at each point
 *  3      0 64		100 1    176 64     251 1
 *  might be a single inherited allele where 3 events occured at loci: 100, 176 and 251. Genetic information is
 * 		derived in the following manner: 
 * 		0-99 		data comes from allele source 64
 * 		100-175 	data comes from allele source 1
 * 		etc....
 */
template <class T>
class AlleleSourceInh: public AlleleSource<T> {
public:
	AlleleSourceInh(AlleleSource<T>* phased, AlleleSource<T>* unphased, vector<size_t>& xoEvents);
	AlleleSourceInh(LocusManager<T>* loci, int sourceID);
	std::string LocusAt(uint locusIndex);
	int At(uint locusIndex);			///<Returns value at given locus
	void At(uint locusIndex, int val);	///<Sets value at given locus

	/**
	 * @Brief write to alleles to binary pool using the most "compressed" manner
	 */
	virtual void WriteBinary(std::ostream& file);

	/**	
	 * @brief Write binary data for each locus
	 */
	virtual void WriteCompleteBinary(std::ostream& file, int sourceID);
	virtual void PrintCrossOvers(std::ostream& os);
	int CountCrossOvers();
	/**
	 * @brief read alleles from binary
	 * @param file Source of the data
	 * @param sourceAlleles living pool to grab pointers from to populate XO lookup
	 */
	virtual void ReadBinary(std::istream& file, Utility::RBTree<size_t, AlleleSource<T>*>* sourceAlleles);

	/**
	 * @brief Resets the internal representation to hold ONLY the model loci
	 */
//	virtual void ResetToMinimal(vector<uint>& modelLoci);

	/**
	 * @Brief For now, we aren't worrying about this function. 
	 * @note Once we have done tests on memory consumption, we may find we have to drop these values too and reload them from wake()
	 */
	void Suspend(vector<uint>& modelLoci) { }
	void Wake(istream& infile) { }
	void ReportSources(std::ostream& os);
	/**
	 * @brief Returns a vector containing the allele sources for a region of an object
	 * @note The lower bound should occur at all times. The upper bound should NEVER be inserted
	 */
	void GetSource(Utility::RBTree<size_t, AlleleSource<T>*>& sources, int lowerBound, int upperBound);

	/**
	 * @Brief Returns a copy of the local object 
	 */
	AlleleSource<T> *Clone();
	/**
	 * @brief Attempt to reconstruct a complete chromosome based on local data
	 */
	AlleleSource<T> *Realize(int sourceID);
	/**
	 * @brief Attempt to reconstruct a complete chromosome based on other
	 */
	AlleleSource<T> *Realize(int sourceID, AlleleSource<T>*other);
	void Reset();
	bool HasCompleteGenotypes();


//	AlleleSource<T>* Cross(AlleleSource<T>* other, vector<size_t>& xoEvents);
	/**
	 * @brief Crosses 2 alleles and returns a new chromosome
	 * @param father The data passed from father. If X/Y, this is assumed to be a Y
	 * @param par Only used for X/Y cross-overs. If is Autosomal, then use NULL
	 * @param returnX only used when par!=NULL. This indicates the return type is X or Y
	 * @param rnd the random number generator
	 * @note if T is LocusXY, returnX is used to determine X/Y of the returned chromosome only if
	 * 			par is non-null (pass NULL if crossing 2 Xs)
	 */
//	AlleleSource<T>* Cross(AlleleSource<T>* other, Utility::Random& rnd, PAR_Region<T> *par=NULL, bool isFemale=true);
	int GetMaxGenotypeCount( );
	void GetSourceData(Utility::RBTree<size_t, AlleleSource<T>*>& srcs);
protected:

	/**
	 * @Brief Indices pointing to the original source for a given allele
	 * @note This MUST be the index into the chromosome, because the XY chromosomes
	 * 			might have different genetic positions for the different chromosomes, and 
	 * 			we can't assume that regions within the PAR are using either X or Y genetic
	 * 			mapping....
	 */
	Utility::RBTree<size_t, AlleleSource<T>*> alleles;
};

template <class T>
void AlleleSource<T>::GetSource(Utility::RBTree<size_t, AlleleSource<T>*>& sources, int lowerBound, int upperBound) {
	if (lowerBound < upperBound) {
assert(lowerBound >= 0);
		sources.Set(lowerBound, this); 
	}
}

/****************************************************************************
 ***      AlleleSource
 ***
 ****************************************************************************/
template <class T>
inline
bool AlleleSource<T>::Compare(AlleleSource<T>* other, int idx) {
	bool identical = locusCount == other->locusCount;
	if (!identical) {
		cerr<<"Ugh, these aren't even the same size!\n";
		return false;
	}
	int width = 5;
	stringstream a, b, c, d;
	a<<idx<<")\t";
	b<<"\t";
	c<<"\t";
	d<<"\t";
	for (int i=0; i<locusCount; i++) {
		a<<setw(width)<<i;
		b<<setw(width)<<At(i);
		d<<setw(width)<<other->At(i);
		bool identicalLoci = At(i) == other->At(i);
		if (!identicalLoci) {
			c<<"    X";
			identical=false;
		}	
		else
			c<<"     ";
	}
	if (!identical)
		cerr<<"\n"<<a.str()<<"\n"
			<<b.str()<<"\n"
			<<c.str()<<"\n"
			<<d.str()<<"\n";
	else
		cerr<<idx<<b.str()<<"\n";
	return identical;
}

template <class T>
inline
size_t AlleleSource<T>::GetLocusCount() {
	return locusCount;
}


template <class T>
inline
int AlleleSourceSingle<T>::GetMaxGenotypeCount( ){
	return 3;
}
template <class T>
inline
void AlleleSource<T>::ShowGenotypes(std::ostream& os, int count, char div){
	if (count > locusCount)
		count = locusCount;
	os<<" "<<div;

	os<<" "<<locusCount<<"  ";
	for (uint i=0; i<count; i++) 
		os<<At(i)<<" ";

	if (count < locusCount) {
		os<<"...";
		for (uint i=locusCount - count; i<locusCount; i++) 
			os<<At(i)<<" ";
	}
}

template <class T>
inline
int AlleleSource<T>::EvaluateKinship(AlleleSource<T> *other) {
	KinshipCalculator<AlleleSource<T>*> k;
	RBTree<size_t, AlleleSource<T>*> localSources;
	RBTree<size_t, AlleleSource<T>*> otherSources;
	GetSourceData(localSources);
	other->GetSourceData(otherSources);
	return k.EvaluateKinship(localSources, otherSources);
}
/****************************************************************************
 ***      AlleleSourceSingle 
 ***
 ****************************************************************************/
template <class T>
inline
AlleleSourceSingle<T>::AlleleSourceSingle(LocusManager<T>* loci, int sourceID) 
		: AlleleSource<T>(loci, sourceID), alleles(NULL) {

}
	

template <class T>
inline
AlleleSourceSingle<T>::AlleleSourceSingle(const AlleleSourceSingle& other) 
	: AlleleSource<T>(other.loci, other.sourceID), alleles(NULL), modelOnlyLoci(other.modelOnlyLoci) {
	if (other.alleles) 
		alleles = new boost::dynamic_bitset<>(*other.alleles);
}

template <class T>
inline
AlleleSourceSingle<T>::~AlleleSourceSingle(){
	if (alleles)
		delete alleles;

}

template <class T>
int AlleleSourceSingle<T>::At(uint locusIndex) {
	if (alleles) {
		return (*alleles)[locusIndex];
	}
	else {
		map<size_t, int>::iterator itr = modelOnlyLoci.find(locusIndex);
		if (itr != modelOnlyLoci.end())
			return itr->second;
		/** This should probably be an exception **/
		cerr<<"AlleleSourceSingle::At("<<locusIndex<<") called with no allelic information. \n";
		assert(0);
		return -1;
	}
}
template <class T>
void AlleleSourceSingle<T>::Wake(istream& file) {
	if (DoSuspend) {
		ReadBinary(file, NULL);
	}
}
template <class T>
void AlleleSourceSingle<T>::Suspend(vector<uint>& modelLoci) {
	if (DoSuspend) {
		ResetToMinimal(modelLoci);
	}
}

template<class T>
void AlleleSourceSingle<T>::GetSourceData(Utility::RBTree<size_t, AlleleSource<T>*>& src) {
	src.Add(0, this);
}
template <class T>
inline
void AlleleSourceSingle<T>::PrintCrossOvers(std::ostream& os) {
	os<<"-\t0) "<<AlleleSource<T>::GetSourceID()<<"\n";
}

template <class T>
inline
int AlleleSourceSingle<T>::CountCrossOvers() {
	return 1;
}

template <class T>
inline
void AlleleSourceInh<T>::GetSourceData(Utility::RBTree<size_t, AlleleSource<T>*>& srcs) {
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	while (node) {
		float pos = node->GetKey();
		int idx = AlleleSource<T>::loci->At(pos)->GetID();
		srcs.Add(idx, node->GetData());
		node = node->GetNext();
	}
	
}

template <class T>
inline
int AlleleSourceInh<T>::CountCrossOvers() {
	return alleles.GetCount();
}
template <class T>
inline
void AlleleSourceInh<T>::PrintCrossOvers(std::ostream& os) {
	int xoCount = 0;
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	os<<"*\t";
	while (node) {
		int pos = node->GetKey();
		AlleleSource<T>* chrom = node->GetData();
		os<<xoCount<<") "<<AlleleSource<T>::loci->At(pos)->GetID()<<","<<chrom->GetSourceID()<<"\t";
		node = node->GetNext();
	}
	os<<"\n";
}

template <class T>
inline
void AlleleSourceSingle<T>::InitLoci(Utility::Random& rnd, bool isY) {
	int count = AlleleSource<T>::loci->LocusCount();
	if (alleles == NULL) 
		alleles = new boost::dynamic_bitset<>(count);
	
	for (int i=0; i<count; i++) {
		T *locus = AlleleSource<T>::loci->At(i);
		alleles->set(i, locus->GetInitialValue(rnd, isY) == 1);
	}
}

template <class T>
inline
void AlleleSourceSingle<T>::At(uint locusIndex, int val) {
	if (alleles == NULL) 
		alleles = new boost::dynamic_bitset<>(AlleleSource<T>::loci->LocusCount());
	assert(alleles&&locusIndex < AlleleSource<T>::loci->LocusCount());
	assert(val < 2 && val >=0);

	(*alleles)[locusIndex] = val > 0;
}

template <class T>
inline
AlleleSource<T> *AlleleSourceSingle<T>::Realize(int sourceID, AlleleSource<T>* other) {
	AlleleSourceSingle *newSource = new AlleleSourceSingle(other->GetLoci(), sourceID);
	for (int i=0; i<AlleleSource<T>::locusCount; i++) {
		newSource->At(i, other->At(i));
	}
	return newSource;
}
/**
 * @brief Performs cross-over with "other".
 */
template <class T>
inline
AlleleSource<T>* AlleleSource<T>::Cross(AlleleSource<T>* other, vector<size_t>& xoEvents) {
	return new AlleleSourceInh<T>(this, other, xoEvents);
}

template <class T>
inline
AlleleSource<T> *AlleleSource<T>::Cross(AlleleSource<T>* father, Utility::Random& rnd, PAR_Region<T>* par, bool isFemale) {
	vector<size_t> xoEvents;
	bool phase;

	//X Y -> ??
	//We are assuming local chromosome is an X
	if (par) {
		phase = par->BuildXOEventList(rnd, xoEvents);
		if (!isFemale)
			phase=!phase;
	}
		
	//X X -> ??
	else {
		phase = AlleleSource<T>::loci->BuildXOEventList(rnd, xoEvents, false);
		if (rnd.drand() < 0.5)
			phase = !phase;
	}
	


	AlleleSource<T> *newChrom = NULL;
	if (phase)
		newChrom = new AlleleSourceInh<T>(father, this, xoEvents);
	else
		newChrom = new AlleleSourceInh<T>(this, father, xoEvents);
	return newChrom;
}

template <class T>
inline
void AlleleSourceSingle<T>::ReadBinary(std::istream& file, Utility::RBTree<size_t, AlleleSource<T>*>* sourceAlleles) {
	uint totalBlocks = (AlleleSource<T>::locusCount / boost::dynamic_bitset<>::bits_per_block) + 1;
	boost::dynamic_bitset<>::block_type raw[totalBlocks];
	
	if (alleles)
		delete alleles;

	assert(!file.eof());
	file.read((char*)&raw, (sizeof(boost::dynamic_bitset<>::block_type)*totalBlocks));
	alleles = new boost::dynamic_bitset<>(&raw[0], &raw[totalBlocks]);
	alleles->resize(AlleleSource<T>::locusCount);
}
/**
 * @Brief write to alleles to binary pool using the most "compressed" manner
 */
template <class T>
inline
void AlleleSourceSingle<T>::WriteBinary(std::ostream& file) {
	assert(HasCompleteGenotypes());
	uint totalBlocks = (alleles->size() / boost::dynamic_bitset<>::bits_per_block) + 1;

	boost::dynamic_bitset<>::block_type raw[totalBlocks];

	//Convert the bitset to raw data
	to_block_range(*alleles, raw);
	file.write((char*)&raw, (sizeof(boost::dynamic_bitset<>::block_type)*totalBlocks));	
}

/**	
 * This is the same as the other write, so no need to change anything
 */
template <class T>
inline
void AlleleSourceSingle<T>::WriteCompleteBinary(std::ostream& file, int sourceID) {
	WriteBinary(file);
}



template <class T>
inline
AlleleSource<T> *AlleleSourceSingle<T>::Clone() {
	AlleleSource<T> *newSource = new AlleleSourceSingle(*this);
	return newSource;
}

template <class T>
inline
void AlleleSourceSingle<T>::Reset() {
	if (alleles)
		delete alleles;
	modelOnlyLoci.clear();
}

/**
 * @brief Resets the internal representation to hold ONLY the model loci
 */
template <class T>
inline
void AlleleSourceSingle<T>::ResetToMinimal(vector<uint>& modelLoci) {
	vector<uint>::iterator itr = modelLoci.begin();
	vector<uint>::iterator end = modelLoci.end();

	modelLoci.clear();
	while (itr != end) {
		modelLoci[*itr] = At(*itr);
		itr++;
	}

	delete alleles;
	alleles=NULL;
}

/****************************************************************************
 ***      AlleleSourceSingle 
 ***
 ****************************************************************************/
template <class T>
inline
AlleleSourceInh<T>::AlleleSourceInh(AlleleSource<T>* phased, AlleleSource<T>* unphased, vector<size_t>& xoEvents) : AlleleSource<T>(phased->GetLoci(), -1){
	vector<size_t>::iterator itr = xoEvents.begin();
	vector<size_t>::iterator end = xoEvents.end();
	AlleleSource<T> *phases[2];
	phases[0]=phased;
	phases[1]=unphased;
	bool phase=false;
	size_t curEvent=0;
	alleles.Clear();
//	alleles.Add(0, phased);
	while (itr != end) {
		int event = *itr++;
		phases[phase]->GetSource(alleles, curEvent, event);
		curEvent = event;
		phase = !phase;
	}
//	if (curEvent != 0)
		phases[phase]->GetSource(alleles, curEvent, AlleleSource<T>::locusCount);


}

template <class T>
inline
AlleleSourceInh<T>::AlleleSourceInh(LocusManager<T>* loci, int sourceID) : AlleleSource<T>(loci, sourceID) {
	
}
/*
template <class T>
inline
AlleleSource<T> *AlleleSourceInh<T>::Cross(AlleleSource<T>* other, Utility::Random& rnd, PAR_Region<T>* par, bool isFemale) {
	vector<size_t> xoEvents;
	bool phase;
	if (par) 
		phase = par->BuildXOEventList(rnd, xoEvents);
	else 
		phase = AlleleSource<T>::loci->BuildXOEventList(rnd, xoEvents, isFemale);
	

	if (rnd.drand() < 0.5)
		phase = !phase;

	AlleleSource<T> *newChrom = NULL;
	if (phase)
		newChrom = new AlleleSourceInh(other, this, xoEvents);
	else
		newChrom = new AlleleSourceInh(this, other, xoEvents);
	return newChrom;
}


template <class T>
inline
AlleleSource<T> *AlleleSourceInh<T>::Cross(AlleleSource<T>* other, vector<size_t>& xoEvents) {
	AlleleSourceInh *newAlleleSource = new AlleleSourceInh(this, other, xoEvents);
	return newAlleleSource;
}
*/
template <class T>
inline
int AlleleSourceInh<T>::GetMaxGenotypeCount( ){
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	if (node)
		return node->GetData()->GetMaxGenotypeCount();
	else
		return -1;
}

template <class T>
inline
void AlleleSourceInh<T>::ReadBinary(std::istream& file, Utility::RBTree<size_t, AlleleSource<T>*> *sourceAlleles) {
	size_t xoCount = 0; 
	unsigned int poolSize = sourceAlleles->GetCount();
	file.read((char*)&xoCount, 4);

	int xoEvent;
	int source;
	for (int i=0; i<xoCount; i++) {
		file.read((char*)&xoEvent, 4);
		file.read((char*)&source, 4);

//		Utlity::RBTreeNode<int, AlleleSource<T>*> *node = alleles.Find(xoEvent);
		if (alleles.Find(xoEvent)==NULL) {
			Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = sourceAlleles->Find(source);
assert(node);
			if (node)
				alleles.Add(xoEvent, node->GetData());
		} else {
			cout<<"ERROR: A problem was encountered trying to load the inherited pool\n";
			exit(1);
		}
	}

}
	

/**
 * @Brief write to alleles to binary pool using the most "compressed" manner
 */
template <class T>
inline
void AlleleSourceInh<T>::WriteBinary(std::ostream& file) {
	size_t xoCount = alleles.GetCount();
	file.write((char*)&xoCount, 4);

	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
/*	if (node->GetKey() != 0) {
		cout<<"An error was ecountered trying to write binary for inherited allele, "<<AlleleSource<T>::sourceID<<". The first XO event occurs at: "<<node->GetKey()<<". We need to know who do we start with!\n";
		exit(1);
	}*/
	while (node) {
		size_t location = node->GetKey();
		int source = node->GetData()->GetSourceID();
		file.write((char*)&location, 4);
		file.write((char*)&source, 4);

		node = node->GetNext();
	}

}
	
template <class T>
inline
void AlleleSourceInh<T>::Reset() {
	alleles.Clear();
//	modelOnlyLoci.clear();
}
/**
 * Realize a new chromosome based on other
 */
template <class T>
inline
AlleleSource<T> *AlleleSourceInh<T>::Realize(int sourceID, AlleleSource<T>* other) {
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	AlleleSource<T> *newSource = NULL, *source = NULL;
	if (node) {
		source = node->GetData();
		newSource = source->Realize(sourceID, other);
	}
	return newSource;
}	

template <class T>
inline
AlleleSource<T> *AlleleSourceInh<T>::Realize(int sourceID) {
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	AlleleSource<T> *newSource = NULL, *source = NULL;
	if (node) {
		source = node->GetData();
		newSource = source->Realize(sourceID, this);
	}
	return newSource;
}
/**	
 * This is the same as the other write, so no need to change anything
 */
template <class T>
inline
void AlleleSourceInh<T>::WriteCompleteBinary(std::ostream& file, int sourceID) {
	//I need to convert to a realized chromosome (i.e. AlleleSourceSingle) and call it's Write
	AlleleSource<T> *realized = Realize(sourceID);
	realized->WriteBinary(file);
	delete realized;
}

template <class T>
inline
AlleleSource<T> *AlleleSourceInh<T>::Clone() {
	AlleleSource<T> *newSource = new AlleleSourceInh(*this);
	assert(((AlleleSourceInh<T>*)newSource)->alleles.GetCount() > 0);
	return newSource;
}


template <class T>
inline
bool AlleleSourceInh<T>::HasCompleteGenotypes() {
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	if (node) {	
		AlleleSource<T> *source = node->GetData();

		if (source)
			return source->HasCompleteGenotypes();
	}
	return false;
}


/**
 * @brief Returns a vector containing the allele sources for a region of an object
 * @note The lower bound should occur at all times. The upper boud should NEVER be inserted
 */

template <class T>
inline
void AlleleSourceInh<T>::GetSource(
		Utility::RBTree<size_t, AlleleSource<T>*>& sources, 
		int lowerBound, 
		int upperBound) {
	typename Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	
	assert(lowerBound <=upperBound);
	while (node) {
		size_t cur = node->GetKey();
assert(cur >= 0);
		if (cur >= lowerBound && cur <= upperBound) {
			sources.Set(cur,node->GetData());
		}
		else if (cur < lowerBound && cur < upperBound) {
			sources.Set(lowerBound, node->GetData());
		}
		node = node->GetNext();
	}
}



template <class T>
inline
std::string AlleleSourceInh<T>::LocusAt(uint locusIndex) {
	T *loc = (*AlleleSource<T>::loci)[locusIndex];
	return loc->GetLabel();
}
///<Returns value at given locus
template <class T>
inline
int AlleleSourceInh<T>::At(uint locusIndex) {
	assert(alleles.GetCount() > 0);	
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.FindNearestMin(locusIndex);

	AlleleSource<T> *source = node->GetData();
	assert(source);
//cerr<<"At("<<locusIndex<<":"<<loc->GetLabel()<<") "<<source->At(locusIndex)<<"\n";
	return source->At(locusIndex);	
}	

template <class T>
inline
void AlleleSourceInh<T>::ReportSources(std::ostream& os) {
	int size = alleles.GetCount();
	Utility::RBTreeNode<size_t, AlleleSource<T>*> *node = alleles.GetFirst();
	os<<size<<" : ";
	while (node) {
		os<<"("<<node->GetKey()<<") "<<node->GetData()->GetSourceID()<<" ";
		node = node->GetNext();
	}
	os<<"\n";
}

template<class T>
inline
void AlleleSourceInh<T>::At(uint locusIndex, int val) {
	cerr<<"Error: Attempt to set value at inherited allele source (allelesource.h:600)\n";
	exit(1);
}	




#ifdef CPPUNIT
class AlleleSourceTest : public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( AlleleSourceTest );
	CPPUNIT_TEST( TestInitialization );
	CPPUNIT_TEST( TestCrossOver );
	CPPUNIT_TEST( TestCrossOverXY );
	CPPUNIT_TEST_SUITE_END();
public:
	AlleleSourceTest();
	~AlleleSourceTest();

	void setUp();
	void tearDown();

	void TestCrossOverXY();
	void TestInitialization();
	void TestCrossOver();
};
#endif

}




#endif
