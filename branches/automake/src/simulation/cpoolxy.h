//
// C++ Interface: cpoolxy
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C)Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONCPOOLXY_H
#define SIMULATIONCPOOLXY_H

#include "cpool.h"
#include "locusxy.h"
#include "par_region.h"
#include "utility/random.h"
#include "individual.h"

namespace Simulation {
using namespace Utility;

#define OOOOPS(c) std::cerr<<"Ooops, the code for this section has not been written. Please send a message to genomeSIMLA@chgr.mc.vanderbilt.edu with the subject: genomeSIMLA Runtime error:\n"<<c<<"\n"
/**
@brief Specialized pool for the XY chromosome(s)
	@author 
*/
class CPoolXY : public CPool<LocusXY>
{
public:
	friend class CPoolXYTest;

	/**
	 * @brief Constructor...
	 * @param chromID Index of this chromosome
	 * @param project Prefix used to write/read all files
	 * @param label Used to help user identify a given chromosome in reports and filenames
	 */
    CPoolXY(int chromID, const char *project, const char *label);
	/**
	 * @brief Pool takes ownership of the locus array...so beware
	 */
    ~CPoolXY();

	/**
	 * @Brief Stores locus array and materializes the PAR
	 */
	void Init(LocusManager<LocusXY>* loci, Random& rnd);

	bool ContinueGrowingY();
	bool ContinueGrowingX();

	int GetPopulationSize();
	int GetPopulationSizeX();
	int GetPopulationSizeY();

	void Refresh(int generation);
	void Bake();
	int Bake(Population<LocusXY>::Type *source, Population<LocusXY>::Type *cur, int &sourceID);
	size_t CalculateAlleleFrequencies(uint maxIndividualCount);
	size_t CalculateAlleleFrequencies(Population<LocusXY>::Type::iterator start, Population<LocusXY>::Type::iterator end, bool isY);
	void Open(LocusManager<LocusXY>* loci, Random& rnd);
	void Wake();
	void Suspend();
	void Save(int generation);

	uint AdvanceGenerations(uint generations, PopulationGrowth::GrowthRate* f, Utility::Random& rnd, uint tCount);
	/**
	 * @Brief Populates the pool based on source pool using cross-over
	 */
	virtual void PopulatePool(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *destX, Population<LocusXY>::Type *destY, Utility::Random& rnd);
	/**
	 * @Brief Threaded pool population (generational advancement)
	 */
	static void *PopulatePool(void *args);

	void DrawY(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *dest, Utility::Random& rnd);
	void DrawX(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *dest, Utility::Random& rnd);
	AlleleSource<LocusXY> *DrawX(Utility::Random& rnd);
	AlleleSource<LocusXY> *DrawY(Utility::Random& rnd);

	void DrawIndividual(Individual& ind, bool isXX, Utility::Random& rnd);

	/**
	 * @brief Builds initial population based on random draws against the minor allele frequencies
	 */
	virtual void BuildInitialPopulation(Random& rnd, int expressionCount);

	void WriteLdData(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh );	
	void WriteLdDataX(const char *basefilename, ostream &summary, Visualization::LocusReport& locReport, float thresh);
	void WriteLdDataY(const char *basefilename, ostream &summary, Visualization::LocusReport& locReport, float thresh);

	void DebugPrint(ostream& os, int start, int stop, Population<LocusXY>::Type* popX=NULL, Population<LocusXY>::Type* popY=NULL);

	bool VerifyAlleleFrequencies();

	/**
	 * @Sanity check to make sure that all non X alleles are fixed
	 */
	bool VerifyXPool(Population<LocusXY>::Type* pop);
	/**
	 * @Sanity check to make sure that all non Y alleles are fixed
	 */
	bool VerifyYPool(Population<LocusXY>::Type* pop);
	void WriteHaploview(ostream& os, Population<LocusXY>::Type *pop, int count);

	bool VerifyChromosomeIs(AlleleSource<LocusXY>* chromosome, bool isX);
protected:
	void BuildPAR();								///<Scan the locus array object for PAR regions
	bool BuildXOEvents(vector<size_t>& events, Random& rnd, float poisson=0);

	/** @brief Builds XO points where lambda is adjusted for gender */
	bool BuildXOEvents(vector<size_t>& events, Random& rnd, bool isMale);	

	/** @brief The founding population (Y chromosomes) */
	vector<AlleleSource<LocusXY> *> *sourcePopulationY;		

	/** @brief Y population at current generation */
	vector<AlleleSource<LocusXY> *> *curPopulationY;

	PAR_Region<LocusXY> par;						///<For XY only, this is the region where XO occurs

	boost::dynamic_bitset<> defaultBits;
	int remainingChromsY;							///<Number of chromosomes to be added to a growing pool


};


#ifdef CPPUNIT
class CPoolXYTest : public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( CPoolXYTest );
	CPPUNIT_TEST( TestInitialization );
	CPPUNIT_TEST_SUITE_END();
public:
	CPoolXYTest();
	~CPoolXYTest();

	void setUp();
	void tearDown();

	void TestInitialization();
};
#endif

}

#endif
