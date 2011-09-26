//
// C++ Implementation: cpoolxy
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "cpoolxy.h"
#include "lineargrowth.h"				//For unit tests
#include "blocklistnodefourgammetes.h"
#include "ldcalculator.h"
#include "ldplotmanager.h"
#include "chrompool.h"					///Just for the static members
#include <boost/timer.hpp>
#include "ldhtmlreport.h"


using namespace Utility;
using namespace Simulation::Visualization;
	

namespace Simulation {
	enum LocusType {
		Normal = 0, 
		X_Only,
		Y_Only,
		X_Homolog,
		Y_Homolog,
		PAR
	};

CPoolXY::CPoolXY(int chromID, const char *project, const char *label) : CPool<LocusXY>(chromID, project, label), curPopulationY(NULL), sourcePopulationY(NULL) {
	
}

CPoolXY::~CPoolXY() {
	if (sourcePopulationY) {
		Purge(sourcePopulationY);
		delete sourcePopulationY;
	}

	if (curPopulationY) {
		Purge(curPopulationY);
		delete curPopulationY;
	}
}
void CPoolXY::Init(LocusManager<LocusXY>* loci, Random& rnd) {
	this->loci = loci;
	par.Init(loci);
	BuildPAR();
	int count = loci->LocusCount();

	//Build out the default loci. I'm using the opposite allele chromosome's allele frequency, because I want to 
	//make sure it's not 0.0 or 1.0. As far as I'm aware, the actual value is irrelavent
	defaultBits.resize(count, false);
	for (int i=0; i<count; i++) {
		LocusXY *loc = loci->At(i);
		if (loc->type == LocusXY::X_Only || loc->type == LocusXY::X_Homolog)
			defaultBits[i] = rnd.drand() < loc->Freq2X();
		else
			defaultBits[i] = rnd.drand() < loc->Freq2Y();
	}
	gdStart = loci->At(0)->MapPosition();
	int locusCount = loci->LocusCount();
	gdLength = (loci->At(locusCount -1)->MapPosition()-gdStart);
	poissonLambda = gdLength * 0.01;
}

void CPoolXY::BuildPAR( ) {
	LocusXY *prev = NULL;
	LocusXY *cur = NULL;

	int count = loci->LocusCount();

	for (int i=0; i<count; i++) {
		LocusXY *loc = loci->At(i);
			prev = cur;
			cur = loc;
		if (loc->type == LocusXY::PAR) {
			if (prev && prev->type == LocusXY::PAR)
				par.InsertRegion(*prev, *cur);
		}
	}
}

bool CPoolXY::ContinueGrowingY() {
	bool doContinue = false;
	LOCKPOOL;
	doContinue = remainingChromsY-- > 0;
	UNLOCKPOOL;
	return doContinue;
}
bool CPoolXY::ContinueGrowingX() {
	bool doContinue = false;
	LOCKPOOL;
	doContinue =  remainingChroms-- > 0;
	UNLOCKPOOL;
	return doContinue;
}

bool CPoolXY::BuildXOEvents(vector<size_t>& events, Utility::Random& rnd, float lambda) {
	if (lambda == 0.0)
		lambda = poissonLambda;
	size_t eventCount = PoissonEventCount(rnd, lambda);

	for (size_t i=0; i<eventCount; i++) {
		float loc = gdStart +  rnd(gdLength);
		Locus *locus = loci->At(loc);
		int idx = locus->GetID();
		vector<size_t>::iterator other = find(events.begin(), events.end(), idx);
		if (other != events.end())
			events.erase(other);
		events.push_back(idx);
	}
	sort(events.begin(), events.end());
	return events.size()%2;
}
/**
 * @brief when we are performing XO when gender is relevant, we adjust it according to (f=1.25N, m=0.75N) where N is the average lambda
 */
bool CPoolXY::BuildXOEvents(vector<size_t>& events, Utility::Random& rnd, bool isMale) {
	float lambda;
	if (isMale)
		lambda = poissonLambda * 0.75;
	else
		lambda = poissonLambda * 1.25;
	return BuildXOEvents(events, rnd, lambda);
}	


void CPoolXY::DrawIndividual(Individual& ind, bool isXX, Utility::Random& rnd) {
	if (isXX) {
		ind.SetXX(DrawX(rnd), DrawX(rnd), &par);
	} else {
		ind.SetXY(DrawX(rnd), DrawY(rnd), &par);
	}
}

AlleleSource<LocusXY> *CPoolXY::DrawX(Utility::Random& rnd) {
	return CPool<LocusXY>::Draw(rnd);
}

AlleleSource<LocusXY> *CPoolXY::DrawY(Utility::Random& rnd) {
	Population<LocusXY>::Type *pop = curPopulationY;
	if (pop == NULL) 
		pop = sourcePopulationY;

	AlleleSource<LocusXY> *source = pop->at(rnd((long int)pop->size()));
	//we want to create a basic copy, with no allelic stuff to copy, so that nothing is wasted if we can't use it
	vector<size_t> empty;
	return source->Cross(source, empty);
}
bool CPoolXY::VerifyChromosomeIs(AlleleSource<LocusXY>* chromosome, bool isX) {
	stringstream ss1, ss2;
	int locCount = loci->LocusCount();
	bool allClear = true;

	for (int i=0; i<locCount; i++) {
		LocusXY *loc = loci->At(i);
		ss1<<chromosome->At(i)<<" ";
		if (isX) {
			if (loc->type == LocusXY::Y_Only || loc->type == LocusXY::Y_Homolog) {
				if (chromosome->At(i)) {
					ss2<<"X ";
					allClear = false;
				}
				else 
					ss2<<"  ";
			}
			else 
				ss2<<"  ";
		
		} else {
			if (loc->type == LocusXY::X_Only || loc->type == LocusXY::X_Homolog) {
				if (chromosome->At(i)) {
					ss2<<"X ";
					allClear = false;
				}
				else 
					ss2<<"  ";
			}
			else 
				ss2<<"  ";
		}

	}
	if (!allClear)
		cerr<<ss1.str()<<"\n"<<ss2.str()<<"\n";
	return allClear;
}
void CPoolXY::DrawX(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *dest, Utility::Random& rnd) {
	int originalSizeX = sourceX->size();
	int originalSizeY = sourceY->size();
//	int count = 0;
//	size_t eventCount = 0;
	while (ContinueGrowingX()) {		
//		count++;
		int matID = rnd(originalSizeX);
	
		vector<size_t> events;
		AlleleSource<LocusXY> *newChrom=NULL;
		//X Y -> X
		if (rnd.drand() < 0.5) {
			bool phase = par.BuildXOEventList(rnd, events);
//			eventCount+=events.size();
			int patID = rnd(originalSizeY);
			if (phase)
				newChrom = sourceY->at(patID)->Cross(sourceX->at(matID), events);
			else
				newChrom = sourceX->at(matID)->Cross(sourceY->at(patID), events);
			assert(VerifyChromosomeIs(newChrom, true));


		} else {
		//X X -> X
			BuildXOEvents(events, rnd, false);	
//			eventCount+=events.size();
			int patID = rnd(originalSizeX);
			if (rnd.drand() < 0.5)
				newChrom = sourceX->at(matID)->Cross(sourceX->at(patID), events);
			else
				newChrom = sourceX->at(patID)->Cross(sourceX->at(matID), events);
			assert(VerifyChromosomeIs(newChrom, true));

		}
		dest->push_back(newChrom);
	}
//cerr<<"X Complete. Avg Event Count: "<<(float)eventCount/(float)count<<"\n";
}

void CPoolXY::DrawY(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *dest, Utility::Random& rnd) {
	int originalSizeX = sourceX->size();
	int originalSizeY = sourceY->size();
	int count = 0;	
	int matCount  = 0;
	int totCount = 0;
//	size_t eventCount = 0;
	while (ContinueGrowingY()) {
		int matID = rnd(originalSizeX);
		int patID = rnd(originalSizeY);
		vector<size_t> events;
	
		AlleleSource<LocusXY> *newChrom=NULL;
		AlleleSource<LocusXY> *pat = sourceY->at(patID);
		AlleleSource<LocusXY> *mat = sourceX->at(matID);
		bool phase = par.BuildXOEventList(rnd, events);
//		eventCount+=events.size();
//		count++;
		//X Y -> Y	
		matCount+=phase;
		totCount+=events.size();
		count++;
		if (phase) {
			newChrom = mat->Cross(pat, events);
//cerr<<"F->M ("<<events.size()<<","<<events[0]<<")\t";
		}
		else {
//if (events.size() > 0)
//cerr<<"M->M ("<<events.size()<<")\t";
			newChrom = pat->Cross(mat, events);

		}
		assert(VerifyChromosomeIs(newChrom, false));
		dest->push_back(newChrom);
	}
//cerr<<"\nY Complete. Avg Mat phase: ("<<float (matCount)/(float)count<<"\n";
//cerr<<"Total XO: "<<totCount<<"\n";
//cerr<<"X Complete. Avg Event Count: "<<(float)eventCount/(float)count<<"\n";

}


int CPoolXY::GetPopulationSize() {
	return GetPopulationSizeX() + GetPopulationSizeY();
}
int CPoolXY::GetPopulationSizeX() {
	if (curPopulation)
		return curPopulation->size();
	else
		return sourcePopulation->size();
}

int CPoolXY::GetPopulationSizeY() {
	if (curPopulationY)
		return curPopulationY->size();
	else
		return sourcePopulationY->size();
}

struct AdvancementArgXY {
	Population<LocusXY>::Type *sourcePopX;
	Population<LocusXY>::Type *sourcePopY;
	Population<LocusXY>::Type newPopX;
	Population<LocusXY>::Type newPopY;
	Utility::Random rnd;
	uint countX;
	uint countY;
	CPoolXY *pool;

	AdvancementArgXY(Population<LocusXY>::Type *sourcePopX, Population<LocusXY>::Type *sourcePopY, uint countX, uint countY, CPoolXY *pool, Utility::Random& rnd) : 
					sourcePopX(sourcePopX), sourcePopY(sourcePopY), rnd(rnd), countX(countX), countY(countY), pool(pool) { 
		newPopX.reserve(countX);
		newPopY.reserve(countY);
	}
};
void *CPoolXY::PopulatePool(void *args) {
	AdvancementArgXY *thArg = (AdvancementArgXY*)args;
	thArg->pool->PopulatePool(thArg->sourcePopX, thArg->sourcePopY, &thArg->newPopX, &thArg->newPopY, thArg->rnd);

	return args;
}
void CPoolXY::PopulatePool(Population<LocusXY>::Type *sourceX, Population<LocusXY>::Type *sourceY, Population<LocusXY>::Type *destX, Population<LocusXY>::Type *destY, Utility::Random& rnd) {
	DrawX(sourceX, sourceY, destX, rnd);
	DrawY(sourceX, sourceY, destY, rnd);
}

void CPoolXY::BuildInitialPopulation(Random& rnd, int expressionCount) {
	int yCount = (uint)((float)expressionCount * 0.25);
	int xCount = expressionCount - yCount;

	int count = loci->LocusCount();
	boost::dynamic_bitset<> fixedAlleles(count, false);

	for (int i=0; i<count; i++) 
		fixedAlleles[i]=rnd.drand() < 0.5;

	//Build X chromosomes
	if (sourcePopulation) {
		Purge(sourcePopulation);
		sourceAlleles.Clear();
	}
	else
		sourcePopulation = new vector<AlleleSource<LocusXY> *>();

	for (int i=0; i < xCount; i++) {
		AlleleSource<LocusXY> *chrom = new AlleleSourceSingle<LocusXY>(loci, i);
		((AlleleSourceSingle<LocusXY>*)chrom)->InitLoci(rnd, false);
		assert(VerifyChromosomeIs(chrom, true));
		sourcePopulation->push_back(chrom);
		sourceAlleles.Add(chrom->GetSourceID(), chrom);

	}

	//Build the Y Chromosomes
	//Build X chromosomes
	if (sourcePopulationY) {
		Purge(sourcePopulationY);
	}
	else
		sourcePopulationY = new vector<AlleleSource<LocusXY> *>();

	for (int i=0; i<yCount; i++) {
		AlleleSource<LocusXY> *chrom = new AlleleSourceSingle<LocusXY>(loci, i+xCount);
		((AlleleSourceSingle<LocusXY>*)chrom)->InitLoci(rnd, true);
		assert(VerifyChromosomeIs(chrom, false));
		sourcePopulationY->push_back(chrom);
		sourceAlleles.Add(chrom->GetSourceID(), chrom);

	}

}



uint CPoolXY::AdvanceGenerations(uint generations, PopulationGrowth::GrowthRate* f, Utility::Random& rnd, uint tCount) {
	tCount--;
	
	bool continueRunning = true;
	uint goalPop = targetPop;
	if (goalPop == 0)
		goalPop = (uint)-1;

	Population<LocusXY>::Type *newPopulationX = new vector<AlleleSource<LocusXY> *>();
	Population<LocusXY>::Type *newPopulationY = new vector<AlleleSource<LocusXY> *>();

	Population<LocusXY>::Type *curX = curPopulation;
	Population<LocusXY>::Type *curY = curPopulationY;
	if (curX == NULL)
		curX = sourcePopulation;
	if (curY == NULL)
		curY = sourcePopulationY;
		//cur = new vector<AlleleSource<T>*>(*sourcePopulation);
	uint popSize = curY->size() + curX->size();
	Population<LocusXY>::Type *nextX = newPopulationX;
	Population<LocusXY>::Type *nextY = newPopulationY;
	
	for (uint currGen = 0; currGen < generations && popSize < goalPop && continueRunning; currGen++) {
		if (nextX == sourcePopulation) {
			nextX = new vector<AlleleSource<LocusXY>*>();
			nextY = new vector<AlleleSource<LocusXY>*>();
		}
		else {
			Purge(nextX);
			Purge(nextY);
		}

		uint endSize = (*f)(++generationCount);
		uint endY = (uint)((float)endSize * 0.25);
		uint endX = endSize - endY;
		if (endSize > goalPop)
			endSize = goalPop;

		remainingChroms = endX;
		remainingChromsY = endY;
		cout<<".";cout.flush();

		pthread_t threads[tCount];
		for (uint th=0; th<tCount; th++) {
			AdvancementArgXY *arg = new AdvancementArgXY(curX, curY, endX, endY, this, rnd);
			pthread_create(&threads[th], NULL, PopulatePool, (void*)arg);
		}
		nextX->reserve(endX);
		nextY->reserve(endY);
		PopulatePool(curX, curY, nextX, nextY, rnd);
	
		for (uint th=0; th<tCount; th++) {
			void *rtn;
			pthread_join(threads[th], &rtn);
			AdvancementArgXY *arg = (AdvancementArgXY *)rtn;
			nextX->insert(nextX->end(), arg->newPopX.begin(), arg->newPopX.end());
			nextY->insert(nextY->end(), arg->newPopY.begin(), arg->newPopY.end());
			delete arg;
		}
		if (nextX->size() > endX)
			nextX->resize(endX);
		if (nextY->size() > endY)
			nextY->resize(endY);
		popSize = endSize;
		Population<LocusXY>::Type *x = curX;
		Population<LocusXY>::Type *y = curY;
		curX = nextX;
		curY = nextY;
		nextX = x;
		nextY = y;
	}
	if (nextX != sourcePopulation) {
		Purge(nextX);
		delete nextX;
	}
	if (nextY != sourcePopulationY) {
		Purge(nextY);
		delete nextY;
	}
	this->curPopulation = curX;
	this->curPopulationY = curY;
cerr<<"Population X("<<curPopulation->size()<<")\tY("<<curPopulationY->size()<<")\n";
	CPool<LocusXY>::curGeneration+=generations;
	return generationCount;
}

void CPoolXY::DebugPrint(ostream& os, int start, int stop, Population<LocusXY>::Type* popX, Population<LocusXY>::Type* popY) {
	Population<LocusXY>::Type *pop = popX;
	if (popY==NULL) {
		pop = curPopulation;
		if (pop == NULL) 
			pop = sourcePopulation;
	}

	if (stop > pop->size())
		stop = pop->size();

	os<<"X: \n";
	for (int i=start; i<stop; i++) {	
		os<<"-- ";
		pop->at(i)->ShowGenotypes(os, 10, '\t');
		os<<"\n";
	}

	pop = popY;
	if (popY == NULL) {
		pop = curPopulationY;
		if (pop == NULL)
			pop = sourcePopulationY;
	}
	if (stop > pop->size())
		stop = pop->size();

	os<<"Y: \n";
	for (int i=start; i<stop; i++) {
		os<<"-- ";
		pop->at(i)->ShowGenotypes(os, 10, '\t');
		os<<"\n";
	}
		
}


void CPoolXY::WriteHaploview(ostream& os, Population<LocusXY>::Type *pop, int count) {
	set<int> usedIndices;
	int popCount = pop->size();

	if (count > popCount)
		count = popCount;

	if (count % 2 > 0)
		count--;

	int indID = 0;
	int i = 0;
	while (count > 0) {
		indID++;
		int c1 = i++;
		while (usedIndices.find(c1) != usedIndices.end())
			c1 = Random::globalGenerator(popCount);
		int c2 = i++;	//Random::globalGenerator(popCount);
		while (usedIndices.find(c2) != usedIndices.end())
			c2 = i++; 	//Random::globalGenerator(popCount);

		AlleleSource<LocusXY>* chrom1 = pop->at(c1)->Clone();
		AlleleSource<LocusXY>* chrom2 = pop->at(c2)->Clone();
		
		CPairXY pair(chrom1, chrom2, loci, false);
		os<<indID<<" "<<indID<<" 0 0 0 0 ";
		pair.ShowGenotypesPedigree(os, 0.0, -1, ' ', NULL);
		os<<"\n";
		count-=2;
	}
}

void CPoolXY::Save(int generation) {
	Population<LocusXY>::Type *popX=curPopulation;
	Population<LocusXY>::Type *popY=curPopulationY;
	if (generation == sourceGeneration) {
		popX=sourcePopulation;
		popY=sourcePopulationY;
	}
	
	ofstream file(GetFilename(generation).c_str(), ios::trunc|ios::binary);
	if (!file.good()) {
		stringstream ss;
		ss<<"Error when trying to save chromsoome pool, "<<label<<". The filename "<<GetFilename(generation)<<" doesn't seem to be valid.";
		throw Exception::General(ss.str().c_str());
	}
	
	file.write((char*)&id, 4);
	file.write((char*)&curGeneration, 4);
	assert(curGeneration == generation);
	size_t popSizeX=popX->size(), popSizeY=popY->size();
	file.write((char*)&popSizeX, 4);
	file.write((char*)&popSizeY, 4);
	int locCount = loci->LocusCount();
	file.write((char*)&locCount, 4);

	int popType = (generation != sourceGeneration);
	file.write((char*)&popType, 4);
	//Read the X population first


	Population<LocusXY>::Type::iterator itr = popX->begin();
	Population<LocusXY>::Type::iterator end = popX->end();
	while (itr != end) {
		(*itr)->WriteBinary(file);
		assert(VerifyChromosomeIs(*itr, true));
		itr++;
	}
	itr = popY->begin();
	end = popY->end();
	while (itr != end) {
		(*itr)->WriteBinary(file);
		assert(VerifyChromosomeIs(*itr, false));
		itr++;
	}

}


void CPoolXY::Suspend() {
	//For now, this is the same, but I can see it changing slightly, if we use suspension during advancement and need to use temporary files
	if (DoSuspend) {
		Save(0);
		if (curGeneration > 0)
			Save(curGeneration);

		CPool<LocusXY>::Suspend(sourcePopulation);
		CPool<LocusXY>::Suspend(curPopulation);
	
		CPool<LocusXY>::Suspend(sourcePopulationY);
		CPool<LocusXY>::Suspend(curPopulationY);

		isSuspended = true;
	}
}
void CPoolXY::Wake() {
	if (DoSuspend) {
		Refresh(0);
		if (curGeneration!=0)
			Refresh(curGeneration);
	}
cerr<<"CPoolXY::Wake\n";
}

void CPoolXY::Open(LocusManager<LocusXY>* loci, Random& rnd) {
	Init(loci, rnd);
	if (sourcePopulation) 
		Purge(sourcePopulation);
	else
		sourcePopulation = new vector<AlleleSource<LocusXY> *>();
	if (sourcePopulationY)
		Purge(sourcePopulationY);
	else
		sourcePopulationY = new vector<AlleleSource<LocusXY>*>();

//	Refresh(0);
}

size_t CPoolXY::CalculateAlleleFrequencies(uint maxIndividualCount) {
	
	Population<LocusXY>::Type *popX = curPopulation;
	Population<LocusXY>::Type *popY = curPopulationY;
	if (popX == NULL)
		popX = sourcePopulation;
	if (popY == NULL)
		popY = sourcePopulationY;

	uint indCount = popX->size();
	if (indCount == 0)
		return 0;
	size_t lociCount = loci->LocusCount();

	//If we want to perform "sampled" frequencies, we do that here
	Population<LocusXY>::Type::iterator startX = popX->begin();
	Population<LocusXY>::Type::iterator endX   = startX;
	Population<LocusXY>::Type::iterator startY = popY->begin();
	Population<LocusXY>::Type::iterator endY	= startY;

	if (maxIndividualCount == 0 || maxIndividualCount > indCount) {
		maxIndividualCount = indCount;	
		endX = popX->end();	
		endY = popY->end();
	}
	else {
		for (int i=0; i<maxIndividualCount; i++) {
			endX++;
			if (endY != popY->end())
				endY++;
		}
	}
	int fixedAlleles = CalculateAlleleFrequencies(startX, endX, false);
	fixedAlleles+=CalculateAlleleFrequencies(startY, endY, true);
	return fixedAlleles;
}
	
	
size_t CPoolXY::CalculateAlleleFrequencies(Population<LocusXY>::Type::iterator start, Population<LocusXY>::Type::iterator end, bool isY) {
	assert(start != end);
	uint fixedAlleles = 0;
	
	Population<LocusXY>::Type::iterator itr = start;
	int lociCount = loci->LocusCount();

	Array2D<int> freq(lociCount, 2);
	int indCount = 0;
	Population<LocusXY>::Type::iterator chrom = start;
	while (chrom != end) {
		for (uint i=0; i<lociCount; i++) 
			freq(i, (*chrom)->At(i))++;
		chrom++;
		indCount++;
	}
	for (int i=0; i<lociCount; i++) {
		float freq1 = freq(i,0);
		if (freq1 == 1.0 || freq1 == 0.0) 
			fixedAlleles++;
		double al1 = freq1/(double)indCount;
		double al2 = 0.0;
		LocusXY *locus = loci->At(i);
		if (isY) {
			al2 = al1;
			al1 = locus->Freq1X();
		} else {
			al2 = locus->Freq1Y();
		}
		locus->AssignFreq(al1, al2);
	}

	if (fixedAlleles == lociCount) {
		chrom = start;
		while (start != end){
			cerr<<"Locus Count ("<<lociCount<<")\t";
			(*chrom)->ShowGenotypes(cerr, 20, '>');
			cerr<<"\n";
			start++;
		}
		cerr<<"We have a problem. There are no variants in our pool at all. Every one of our "<<indCount<<" expressions is the same\n";
		abort();
	}

	return fixedAlleles;
}

int CPoolXY::Bake(Population<LocusXY>::Type *source, Population<LocusXY>::Type *cur, int &sourceID) {
	assert(source);
	
	int xoSum = 0;
	if (isSuspended) 
		Wake();

	if (cur) {
		int count = cur->size();	
		for (int i=0; i<count; i++) {
			AlleleSource<LocusXY> *originalSource = cur->at(i);
			xoSum+=originalSource->CountCrossOvers();
			AlleleSource<LocusXY> *realizedSource = originalSource->Realize(++sourceID);
			sourceAlleles.Set(sourceID, realizedSource);
			delete originalSource;
			(*cur)[i] = realizedSource;
		}
		xoSum/=count;
	}
	return xoSum;
}
void CPoolXY::Bake() {

	int sourceID = sourceAlleles.GetLast()->GetKey();
	sourceAlleles.Clear();
	//Do X first
	int xCount = Bake(sourcePopulation, curPopulation, sourceID);
	//Now y
	int yCount = Bake(sourcePopulationY, curPopulationY, sourceID);
	cerr<<"Average XO (X): "<<xCount<<"\nAverage XO (Y): "<<yCount<<"\n";
	//Purge the populations AFTER we've built them both, since there might be some common heritage
	Purge(sourcePopulation);
	Purge(sourcePopulationY);

	//Replace the original sources and update the generation
	sourcePopulation = curPopulation;
	curPopulation = NULL;
	sourcePopulationY = curPopulationY;
	curPopulationY = NULL;
	sourceGeneration = curGeneration;

}
void CPoolXY::Refresh(int generation) {
	Population<LocusXY>::Type *popX=curPopulation;
	Population<LocusXY>::Type *popY=curPopulationY;
	if (generation == sourceGeneration|| popX==NULL) {
		popX=sourcePopulation;
		popY=sourcePopulationY;
		assert(curPopulation == NULL);
		assert(curPopulationY == NULL);
	}
	Purge(popX);
	Purge(popY);
	
	
	ifstream file(GetFilename(generation).c_str(), ios::binary);
	if (!file.good()) {
		stringstream ss;
		ss<<"Error when trying to refresh chromsoome pool, "<<label<<". The file "<<GetFilename(generation)<<" doesn't seem to be valid.";
		throw Exception::General(ss.str().c_str());
	}
	
	file.read((char*)&id, 4);
	file.read((char*)&curGeneration, 4);
	assert(curGeneration == generation);
	size_t popSizeX=0, popSizeY=0;
	file.read((char*)&popSizeX, 4);
	file.read((char*)&popSizeY, 4);
	int locCount = 0;
	file.read((char*)&locCount, 4);
	assert(locCount == loci->LocusCount());
	int popType = 0;
	file.read((char*)&popType, 4);
	//Read the X population first

cerr<<"\nLoad Population X ("<<popSizeX<<"):\n";
	for (int i=0; i<popSizeX; i++) {
		AlleleSource<LocusXY> *chrom;
		if (popType == 0)
			chrom=new AlleleSourceSingle<LocusXY>(loci, i);
		else	
			chrom=new AlleleSourceInh<LocusXY>(loci, i);
		chrom->ReadBinary(file, &sourceAlleles);
		assert(VerifyChromosomeIs(chrom, true));
		popX->push_back(chrom);
	}
	
cerr<<"\nLoad Population Y ("<<popSizeY<<"):\n";
	for (int i=0; i<popSizeY; i++) {
		AlleleSource<LocusXY> *chrom;
		if (popType == 0)
			chrom=new AlleleSourceSingle<LocusXY>(loci, i);
		else 
			chrom=new AlleleSourceInh<LocusXY>(loci, i);
		chrom->ReadBinary(file, &sourceAlleles);
		assert(VerifyChromosomeIs(chrom, false));
		popY->push_back(chrom);
	}

	VerifyAlleleFrequencies();
}

/**
 * 
 */
bool CPoolXY::VerifyAlleleFrequencies() {
	//Start with X
	Population<LocusXY>::Type *pop=curPopulation;
	if (pop == NULL) {
		pop=sourcePopulation;
	}
	int locusCount = loci->LocusCount();
	Population<LocusXY>::Type::iterator itr = pop->begin();
	Population<LocusXY>::Type::iterator end = pop->end();
	Array2D<int> freq(locusCount, 2);
	bool frequenciesMatch = true;
	stringstream ss;
	while (itr != end) {
		for (int i=0; i<locusCount; i++) 
			freq(i, (*itr)->At(i))++;
		itr++;
	}
	for (int i=0; i<locusCount; i++) {
		float countA = freq(i,0);
		float freq1 = 0.0; 
		float counta = freq(i, 1);
		if (countA > 0.0) 
			freq1 = countA /(countA+counta);
		frequenciesMatch = frequenciesMatch && (loci->At(i)->Freq1X() - 0.0001 < freq1 && loci->At(i)->Freq1X() + 0.0001 > freq1);

		LocusXY *loc = loci->At(i);
		if (loc->type == LocusXY::Y_Only || loc->type == LocusXY::Y_Homolog) {
			if (freq1 != 1.0) {
				cerr<<"Invalid Pool contents in X at locus("<<i<<"), freq ("<<freq1<<") (should be fixed)\n";
			}
		}
	}
	if (!frequenciesMatch) {
		cerr<<"There are issues with allele frequencies not matching completely the locus file\n";
		freq.Print(cerr);
	}
	assert(frequenciesMatch);
	return frequenciesMatch;
}

void CPoolXY::WriteLdData(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh ) {
	WriteLdDataX(basefilename, summary, locReport, thresh);
	WriteLdDataY(basefilename, summary, locReport, thresh);
}

void CPoolXY::WriteLdDataX(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh ) {
	vector<Locus> locusData;
	//Get the X Locus array
	int fixedAlleles = loci->GetLocusArray(0, loci->LocusCount(), locusData, 0);
	string sampledDetails = "";
	if (ChromPool::fastLD)
		sampledDetails="-sampled";
	string filePrefix = basefilename + string("X") + sampledDetails;
	string htmlFilename = filePrefix + ".html";
	string locusFilename = filePrefix + ".loc";
	LdHtmlReport report(htmlFilename.c_str(), curGeneration);
	if (ChromPool::fastLD)
		filePrefix += "-sampled";

	Population<LocusXY>::Type *pop = curPopulation;
	if (pop == NULL) 
		pop = sourcePopulation;

	int locusCount = locusData.size();
	int threadCount = 1;				// We can correct this once everyhting else works
	//Do sampling, if we are told to do so
	Population<LocusXY>::Type::iterator popItr = pop->begin();
	Population<LocusXY>::Type::iterator popEnd = pop->end();
	int totalIndividuals = pop->size();

	if (ChromPool::fastLD) {
		int portion = pop->size() / threadCount;
		//We need to determine how many, and then find the correct offsets. 
		if (totalIndividuals > ChromPool::maxLDIndividuals) {
			totalIndividuals = ChromPool::maxLDIndividuals;
			//For now, we're just doing one thread, so we'll just set the end now
			popEnd = popItr;
			int count = 0;
			while (count++ < portion)
				popEnd++;
		}
	}

	BlockListHead<BlockListNodeFourGammetes> haplotypeHead;
	haplotypeHead.Reset(0, locusData.size() -1);
	LdCalculator<LocusXY> ld(popItr, popEnd, &locusData, &haplotypeHead);

	//Each SNP gets it's own vector which is all subsequent SNPs with LD values 
	vector<vector<LdResult> > results(locusCount - 1);
	int locusIdx = 0;
	
boost::timer progress;
	int lastLocus = locusCount-1;
	for (int i=0; i<lastLocus; i++) {
		ld.Init(i, lastLocus);
		ld.CalculateFrequencies();
		vector<LdResult>& result = results[locusIdx++];
		ld.GetLDResults(result);
//cerr<<"LD Results("<<i<<").size() == "<<result.size()<<"\n";
	}
cerr<<setw(6)<<"\t"<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<progress.elapsed()<<"\n";
	haplotypeHead.ValidateBlocks(locusData, NULL);
	BlockList *haplotypes = haplotypeHead.GetValidBlocks();
	sort(haplotypes->begin(), haplotypes->end(), BlockSizeEval());
	report.Init("Y Chromosome", 0.01, totalIndividuals, locusData.size(), fixedAlleles, haplotypes->size());
	LdPlotManager ldPlots(&locusData);;
	ldPlots.Init(results, &haplotypeHead);
	string textPlot = ldPlots.WriteLdReport(filePrefix.c_str(), NULL);

	uint imgHeight, imgWidth;
	string dprimeOver = ldPlots.WriteLdDPrime(filePrefix.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
	string rsOver     = ldPlots.WriteLdRSquared(filePrefix.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
	report.AddSummaryPlots(dprimeOver.c_str(), rsOver.c_str(), textPlot.c_str(), locusFilename.c_str());
	
	int blockCount = haplotypes->size();
	if (blockCount > numberOfBlocksToReport)
		blockCount = numberOfBlocksToReport;

	uint detailedReportCount = ChromPool::numberOfBlocksToReport;
	char plotFilename[4096];
	for (int i=0; i<blockCount; i++) {
		BlockListNode *block = haplotypes->at(i);
		sprintf(plotFilename, "%s.%d", filePrefix.c_str(), i+1);


		if (i<detailedReportCount) {
			string dpO = ldPlots.WriteLdDPrime(plotFilename, block, imgHeight, imgWidth, highBufferSize, "-ovr");
			string dpD = ldPlots.WriteLdDPrime(plotFilename, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			string rsO = ldPlots.WriteLdRSquared(plotFilename, block, imgHeight, imgWidth, highBufferSize, "-ovr");
			string rsD = ldPlots.WriteLdRSquared(plotFilename, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			report.AddBlockReport(block, ExtractFilename(dpO.c_str()).c_str(), ExtractFilename(rsO.c_str()).c_str(), ExtractFilename(dpD.c_str()).c_str(), ExtractFilename(rsD.c_str()).c_str(), imgHeight);
		}
		else
			report.AddBlockSummary(block);
	}
	/********* HAPLOVIEW X Chromosome***********/
	string markerFile = filePrefix + ".info";
	((LocusManagerFileBased<LocusXY>*)loci)->WriteMarkerInfo(markerFile.c_str());
	string hapFile = filePrefix + ".ped";
	ofstream hap(hapFile.c_str());
	WriteHaploview(hap, pop, 600);
	hap.close();
	report.Close();

	LOCKREPORT;
	summary<<"<TR><TD><TABLE><CAPTION>X Chromosome"<<sampledDetails<<"</CAPTION>\n";
	summary<<"\t<TR><TH>SNP Count</TH><TD>"<<locusData.size()<<"</TD>\n";
	summary<<"\t<TR><TH>Fixed Alleles</TH><TD>"<<fixedAlleles<<"</TD>\n";
	summary<<"\t<TR><TH>Block Count</TH><TD>"<<blockCount<<"</TD>\n";
	if (dprimeOver.length() > 0) 
		summary<<"</TABLE></TD><TD><img src='"<<ExtractFilename(dprimeOver.c_str())<<"' width=650 onClick=\"window.open('"<<ExtractFilename(dprimeOver.c_str())<<"');\"></IMG></TD></TR>\n";
	else if (rsOver.length() > 0) 
		summary<<"</TABLE></TD><TD><img src='"<<ExtractFilename(rsOver.c_str())<<"' width=650 onClick=\"window.open('"<<ExtractFilename(rsOver.c_str())<<"');\"></IMG></TD></TR>\n";
		
	UNLOCKREPORT;

	return;
}


void CPoolXY::WriteLdDataY(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh ) {
	vector<Locus> locusData;
	//Get the X Locus array
	int fixedAlleles = loci->GetLocusArray(0, loci->LocusCount(), locusData, 1);
	string sampledDetails = "";
	if (ChromPool::fastLD)
		sampledDetails="-sampled";
	string filePrefix = basefilename + string("Y") + sampledDetails;
	string htmlFilename = filePrefix + ".html";
	string locusFilename = filePrefix + ".loc";
	LdHtmlReport report(htmlFilename.c_str(), curGeneration);
	


	Population<LocusXY>::Type *pop = curPopulationY;
	if (pop == NULL) 
		pop = sourcePopulationY;

	int locusCount = locusData.size();
	if (locusCount < 1) {
		string message = "The Y Population has completely fixed. Unable to generate LD plots.";
		throw Utility::Exception::General(message.c_str());
	}



	int threadCount = 1;				// We can correct this once everyhting else works

	Population<LocusXY>::Type::iterator popItr = pop->begin();
	Population<LocusXY>::Type::iterator popEnd = pop->end();
	int totalIndividuals = pop->size();

	if (ChromPool::fastLD) {
		int portion = pop->size() / threadCount;
		//We need to determine how many, and then find the correct offsets. 
		if (totalIndividuals > ChromPool::maxLDIndividuals) {
			totalIndividuals = ChromPool::maxLDIndividuals;
			//For now, we're just doing one thread, so we'll just set the end now
			popEnd = popItr;
			int count = 0;
			while (count++ < portion)
				popEnd++;
		}
	}


	BlockListHead<BlockListNodeFourGammetes> haplotypeHead;
	haplotypeHead.Reset(0, locusData.size() -1);
	LdCalculator<LocusXY> ld(popItr, popEnd, &locusData, &haplotypeHead);

//	vector<Locus>::iterator locItr = locusData.begin();
//	vector<Locus>::iterator locEnd = locusData.end();

	//Each SNP gets it's own vector which is all subsequent SNPs with LD values 
	vector<vector<LdResult> > results(locusCount - 1);
	int locusIdx = 0;

	int lastLocus = locusCount-1;
	for (int i=0; i<lastLocus; i++) {
		ld.Init(i, lastLocus);
		ld.CalculateFrequencies();
		vector<LdResult>& result = results[locusIdx++];
		ld.GetLDResults(result);
	}
	haplotypeHead.ValidateBlocks(locusData, NULL);
	BlockList *haplotypes = haplotypeHead.GetValidBlocks();
	sort(haplotypes->begin(), haplotypes->end(), BlockSizeEval());

	LOCKREPORT;
	locReport.Init(locusData, haplotypes, ExtractFilename(htmlFilename.c_str()).c_str(), GetLabel().c_str());
	UNLOCKREPORT;


	LdPlotManager ldPlots(&locusData);
	ldPlots.Init(results, &haplotypeHead);
	string textPlot = ldPlots.WriteLdReport(filePrefix.c_str(), NULL);

	uint imgHeight, imgWidth;
	string dprimeOver = ldPlots.WriteLdDPrime(filePrefix.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
	string rsOver     = ldPlots.WriteLdRSquared(filePrefix.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
	
	report.Init("Y Chromosome", 0.01, totalIndividuals, locusData.size(), fixedAlleles, haplotypes->size());
	report.AddSummaryPlots(dprimeOver.c_str(), rsOver.c_str(), textPlot.c_str(), locusFilename.c_str());
	int blockCount = haplotypes->size();
	if (blockCount > numberOfBlocksToReport)
		blockCount = numberOfBlocksToReport;

	uint detailedReportCount = ChromPool::numberOfBlocksToReport;
	char plotFilename[4096];
	for (int i=0; i<blockCount; i++) {
		BlockListNode *block = haplotypes->at(i);
		sprintf(plotFilename, "%s.%d", filePrefix.c_str(), i+1);

		if (i<detailedReportCount) {
			string dpO = ldPlots.WriteLdDPrime(plotFilename, block, imgHeight, imgWidth, highBufferSize, "-ovr");
			string dpD = ldPlots.WriteLdDPrime(plotFilename, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			string rsO = ldPlots.WriteLdRSquared(plotFilename, block, imgHeight, imgWidth, highBufferSize, "-ovr");
			string rsD = ldPlots.WriteLdRSquared(plotFilename, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			report.AddBlockReport(block, ExtractFilename(dpO.c_str()).c_str(), ExtractFilename(rsO.c_str()).c_str(), ExtractFilename(dpD.c_str()).c_str(), ExtractFilename(rsD.c_str()).c_str(), imgHeight);
		}
		report.AddBlockSummary(block);
	}
	/********* HAPLOVIEW Y Chromosome***********/
	string markerFile = filePrefix + ".info";
	((LocusManagerFileBased<LocusXY>*)loci)->WriteMarkerInfo(markerFile.c_str());
	string hapFile = filePrefix + ".ped";
	ofstream hap(hapFile.c_str());
	WriteHaploview(hap, pop, 600);
	hap.close();

	LOCKREPORT;
	summary<<"<TR><TD><TABLE><CAPTION>Y Chromosome"<<sampledDetails<<"</CAPTION>\n";
	summary<<"\t<TR><TH>SNP Count</TH><TD>"<<locusData.size()<<"</TD>\n";
	summary<<"\t<TR><TH>Fixed Alleles</TH><TD>"<<fixedAlleles<<"</TD>\n";
	summary<<"\t<TR><TH>Block Count</TH><TD>"<<blockCount<<"</TD>\n";
	if (dprimeOver.length() > 0) 
		summary<<"</TABLE></TD><TD><img src='"<<ExtractFilename(dprimeOver.c_str())<<"' width=650 onClick=\"window.open('"<<ExtractFilename(dprimeOver.c_str())<<"');\"></IMG></TD></TR>\n";
	else if (rsOver.length() > 0) 
		summary<<"</TABLE></TD><TD><img src='"<<ExtractFilename(rsOver.c_str())<<"' width=650 onClick=\"window.open('"<<ExtractFilename(rsOver.c_str())<<"');\"></IMG></TD></TR>\n";
	UNLOCKREPORT;

	return;
}


#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(CPoolXYTest);


CPoolXYTest::CPoolXYTest() { }

CPoolXYTest::~CPoolXYTest() { }

void CPoolXYTest::setUp() { }

void CPoolXYTest::tearDown() { }

void CPoolXYTest::TestInitialization() {
	Utility::Random rnd;
	string filename("20XYloc.loc");
	WriteLocusXY(filename.c_str());
	LocusManagerFileBased<LocusXY> fb("test", 1);
	fb.Load(filename.c_str());

	//I have no idea how to test randomized stuff. So, I'm going to just exercise 
	//the functions and look at things I can test for. Verifying things like
	//the genotypes inside a given individual, though, is going to be hard to do
	CPoolXY pool(1, "test_pool", "Test Pool");
	pool.Init(&fb, rnd);
	pool.BuildInitialPopulation(rnd, 100);
	pool.Save();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpoolxy::open", 75, (int)pool.GetPopulationSizeX());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpoolxy::open", 25, (int)pool.GetPopulationSizeY());
	//OK, let's try opening the pool as another and make sure that they have the same
	//stuff inside
	CPoolXY pool2(2, "test_pool", "Test Pool");
	pool2.Open(&fb, rnd);
	
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpoolxy::open", 75, (int)pool2.GetPopulationSizeX());
	AlleleSource<LocusXY> *chrom1 = pool.At(0);
	AlleleSource<LocusXY> *chrom2 = pool2.At(0);

	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", 15, (int)chrom2->GetLocusCount());
	int countZeros=0;
	for (int i=0; i<20; i++) {
		CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", chrom1->At(i), chrom2->At(i));
		countZeros+=(chrom1->At(i)==1);
	}
	CPPUNIT_ASSERT_MESSAGE("num zeros", (countZeros > 0 && countZeros < 20));

	//These should be different once we advance into the future
	chrom1 = pool.At(0);
		int a = chrom1->At(0);
	PopulationGrowth::GrowthRate::minPoolSize = 100;
	//Grow the population
	PopulationGrowth::LinearGrowth growth(100, 5, 0.0);				///<Basic growth, 5/generation
	pool.AdvanceGenerations(10, &growth, rnd, 1);
	

	int popY = (int)(150.0 * 0.25);
	int popX = 150 - popY;
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", popY, (int)pool.GetPopulationSizeY());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", popX, (int)pool.GetPopulationSizeX());
	
	chrom2 = pool.At(0);
	
	countZeros = 0;
	int overlap =0;
	for (int i=0; i<15; i++) {
		int a = chrom1->At(i);
		int b = chrom2->At(i);
		overlap += chrom1->At(i) == chrom2->At(i);
		countZeros += chrom2->At(i) == 1;
	}
	CPPUNIT_ASSERT_MESSAGE("cpool::advance generation", (overlap > 0 && overlap < 20));
	CPPUNIT_ASSERT_MESSAGE("cpool::advance generation", (countZeros > 0 && countZeros < 20));
		
	pool.Suspend();
	pool.Wake();

	chrom2 = pool.At(0);

	int zeros=0, overlapping=0;
	for (int i=0; i<15; i++) {
		overlapping += chrom1->At(i) == chrom2->At(i);
		zeros += chrom2->At(i) == 1;
	}


	//Let's test drawing a Y CHromosome
	vector<size_t> events;
	
	AlleleSource<LocusXY> *newChrom=NULL;
	chrom1 = pool.At(0);
	chrom2 = pool.At(1);
	bool phase = false;

	//We want to make sure that we get an odd XO at the first place
	while (!phase) {
		events.clear();
		phase = pool.par.BuildXOEventList(rnd, events);
	}
	newChrom = chrom1->Cross(chrom2, events);
cerr<<"\n-----------------------------------------------\n";
cerr<<chrom1<<"\n";
cerr<<chrom2<<"\n";
cerr<<newChrom<<"\n";
newChrom->PrintCrossOvers(cerr);
cerr<<"\n-----------------------------------------------\n";
	char msg[1024];
	for (int i=3; i<11; i++) {
		sprintf(msg, "CPoolXY::XO (%d)", i);
		CPPUNIT_ASSERT_EQUAL_MESSAGE(msg, newChrom->At(i), chrom2->At(i));
	}


}
#endif //CPPUNIT



}
