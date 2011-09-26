//
// C++ Interface: cpair
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONCPAIR_H
#define SIMULATIONCPAIR_H
#include "allelesource.h"
#include "locusxy.h"

namespace Simulation {

/**
@brief Represents a pair of chromosomes (or more) which combine to for a single genotye at each locus position
		An individual will hold an array or vector of these objects to represent individuals with 1 or more chromosomes

	@author Eric Torstenson
*/
template<class T>
class CPair{
public:
	CPair();
	CPair(const CPair& other);
    CPair(AlleleSource<T>* ch1, AlleleSource<T>* ch2, LocusManager<T>* locusSouce);
    virtual ~CPair();

	void Init(AlleleSource<T>* ch1, AlleleSource<T>* chr2, LocusManager<T>* locusSource);

	virtual bool operator=(const CPair& other);

	/**
	 * @Brief assembles genotype 
	 * @note -1 Missing, 0 aa, 1 Aa, 2 AA (for duplex)
	 */
	virtual int GetGenotype(int locusIdx);

	/**
	 * @brief mark a given locus as missing
	 */
	virtual bool SetMissing(int locus, bool isMissing=false);
	/**
	 * @Brief Writes genotypes to stream
	 * @param os The stream to be written to
	 * @param minAlleleThresh If maf exceeds this threshold, it will be written
	 * @param count The number of SNPs written (first N). -1 will write them all
	 * @param div How the SNPs are divided (Space, tab, comma, etc)
	 */
	virtual void ShowGenotypes(std::ostream& os, float minAlleleThresh=0.0, int count=-1, char div=' ', uint *genotypes= NULL);

	/**
	 * @Brief Writes genotypes to stream in pedigree format
	 * @param os The stream to be written to
	 * @param minAlleleThresh If maf exceeds this threshold, it will be written
	 * @param count The number of SNPs written (first N). -1 will write them all
	 * @param div How the SNPs are divided (Space, tab, comma, etc)
	 */
	virtual void ShowGenotypesPedigree(std::ostream& os, float minAlleleThresh=0.0, int count=-1, char div=' ', uint *genotypes= NULL);	

	AlleleSource<T> *chrom1;				///<First Chromosome
	AlleleSource<T> *chrom2;				///<
	LocusManager<T> *locusSource;
	Utility::BitSetType missingData;

	//draw a new alleleSource based on the local pair
	AlleleSource<T> *Draw(Utility::Random& rnd, PAR_Region<T> *par, bool isFemale);				
};


class CPairXY : public CPair<LocusXY> {
public:
	CPairXY();
	CPairXY(const CPairXY& other);
	CPairXY(AlleleSource<LocusXY>*y, AlleleSource<LocusXY>*y, LocusManager<LocusXY>* locusSource, bool isXX);
	~CPairXY();

	virtual bool operator=(const CPairXY& other);

	/**
	 * @Brief basically, this differs by having to handle X or Y linked traits....
	 */
	virtual int GetGenotype(int locusIdx);


	bool IsXX() { return isXX; }
	void IsXX(bool isXX) { this->isXX = isXX; }

	/**
	 * @brief Returns X/Y based on returnX derived from the two local chromosomes
	 * @note X/X should pass NULL as par
	 */
	AlleleSource<LocusXY> *Draw(Utility::Random& rnd, PAR_Region<LocusXY> *par, bool returnX);				

protected:
	bool isXX;								///<Used in determining XX or XY
};

template <class T>
inline
AlleleSource<T> *CPair<T>::Draw(Utility::Random& rnd, PAR_Region<T> *par, bool isFemale) {
	return chrom1->Cross(chrom2, NULL, isFemale);
}

template <class T>
inline
void CPair<T>::ShowGenotypes(std::ostream& os, float minAlleleThresh, int count, char div, uint *genotypeCounts) {
	int locusCount = chrom1->GetLocusCount();
	LocusManager<T> *loci = chrom1->GetLoci();					
	int maxGenotypes = chrom1->GetMaxGenotypeCount() + 1;	//Extra one for missing data
	if (count < 0 || count > locusCount)
		count  = locusCount;

	for(int i=0; i<count; i++) 
		if (loci->At(i)->PassThresholdMAF(minAlleleThresh) ) {
			int gt = chrom1->At(i) + chrom2->At(i) + 1;
    		os<<gt<<div;
			if (genotypeCounts) {
				genotypeCounts[gt]++;
				//For now, we are assuming 4 genotypes (3 + missing). That might have to change for sequence data
				genotypeCounts+=maxGenotypes;
			}
		}
}

template <class T>
inline
void CPair<T>::ShowGenotypesPedigree(std::ostream& os, float minAlleleThresh, int count, char div, uint *genotypeCounts) {

	static string genotypes[] = { " 0 0", " 1 1", " 1 2", " 2 2" };
	int locusCount = locusSource->LocusCount();
	int maxGenotypes = chrom1->GetMaxGenotypeCount() + 1;	//Extra one for missing data
	if (count < 0 || count > locusCount)
		count  = locusCount;

	for(int i=0; i<count; i++) 
		if (locusSource->At(i)->PassThresholdMAF(minAlleleThresh) ) {
			int gt = chrom1->At(i) + chrom2->At(i) + 1;
    		os<<genotypes[gt]<<div;
			if (genotypeCounts) {
				genotypeCounts[gt]++;
				//For now, we are assuming 4 genotypes (3 + missing). That might have to change for sequence data
				genotypeCounts+=maxGenotypes;
			}
		}
}

template <class T>
inline
CPair<T>::CPair() : chrom1(NULL), chrom2(NULL), locusSource(NULL) {

}

template <class T>
inline
CPair<T>::CPair(AlleleSource<T>* ch1, AlleleSource<T> *ch2, LocusManager<T> *locusSource) : chrom1(ch1), chrom2(ch2), locusSource(locusSource), missingData(locusSource->LocusCount(), false) { 
	if (locusSource == NULL)
		this->locusSource = ch1->GetLoci();
}

template <class T>
inline
void CPair<T>::Init(AlleleSource<T>* ch1, AlleleSource<T>* ch2, LocusManager<T> *locSource) {
	chrom1 = ch1; 
	chrom2 = ch2;
	if (locSource)
		locusSource = locSource;
	else
		locusSource = ch1->GetLoci();

	missingData.resize(locusSource->LocusCount(), false);
}

template <class T>
inline
int CPair<T>::GetGenotype(int locusIdx) {
	assert(locusIdx < chrom1->GetLocusCount());
	if (missingData[locusIdx])
		return -1;
	return chrom1->At(locusIdx) + chrom2->At(locusIdx);
}

template <class T>
inline

bool CPair<T>::SetMissing(int locus, bool isMissing) {
	assert(locus<missingData.size());
	missingData[locus] = isMissing;
}

template <class T>
inline
CPair<T>::~CPair() {
	if (chrom1)
		delete chrom1;
	if (chrom2)
		delete chrom2;
}

template <class T>
inline
CPair<T>::CPair(const CPair<T>& other) : locusSource(other.locusSource) {
	chrom1 = other.chrom1->Clone();
	chrom2 = other.chrom2->Clone();
}

template <class T>
inline
bool CPair<T>::operator=(const CPair<T>& other) {
	chrom1 = other.chrom1->Clone();
	chrom2 = other.chrom2->Clone();
	locusSource = other.locusSource;
	missingData = other.missingData;
}

inline
CPairXY::CPairXY(const CPairXY& other) : CPair<LocusXY>(other.chrom1, other.chrom2, other.locusSource) {
	isXX = other.isXX;
}

inline
bool CPairXY::operator=(const CPairXY& other) {
	chrom1 = other.chrom1->Clone();
	chrom2 = other.chrom2->Clone();
	locusSource = other.locusSource;
	missingData = other.missingData;
	isXX = other.isXX;
}

inline
CPairXY::CPairXY() : isXX(false) { }

inline
CPairXY::CPairXY(AlleleSource<LocusXY>* x, AlleleSource<LocusXY>* y, LocusManager<LocusXY>* locusSource, bool isXX) : 
		CPair<LocusXY>(x, y, locusSource) {
	this->isXX=isXX;
}

inline
CPairXY::~CPairXY() { }




inline
AlleleSource<LocusXY> *CPairXY::Draw(Utility::Random& rnd, PAR_Region<LocusXY> *par, bool returnX) {
	return chrom1->Cross(chrom2, rnd, par, returnX);
}
inline
int CPairXY::GetGenotype(int locusIdx) {
	LocusXY *locus = CPair<LocusXY>::locusSource->At(locusIdx);
	if (missingData[locusIdx])
		return -1;
	if (locus) {
		if (isXX) {
			if (locus->type == LocusXY::Y_Only)
				return chrom2->At(locusIdx)*2;				///<Just let it be homozygous whatever is at the Y
			else
				return CPair<LocusXY>::GetGenotype(locusIdx);
		} else {	
			switch (locus->type) {
				case LocusXY::X_Only:
					return chrom1->At(locusIdx)*2;
				case LocusXY::Y_Only:
					return chrom2->At(locusIdx)*2;
				default:
					return CPair<LocusXY>::GetGenotype(locusIdx);
			}
		}
	}
	return -1;
}

#ifdef CPPUNIT
class CPairTest : public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( CPairTest );
	CPPUNIT_TEST( Test );
	CPPUNIT_TEST( TestXY );
	CPPUNIT_TEST_SUITE_END();
public:
	CPairTest();
	~CPairTest();

	void setUp();
	void tearDown();

	void Test();
	void TestXY();
};
#endif

}

#endif
