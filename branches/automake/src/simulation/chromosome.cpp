// Chromosome.cpp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////

#include "chromosome.h"
#include "utility/types.h"
#include "minchrom.h"
#include <sstream>
#include "kinshipcalculator.h"
namespace Simulation {

using namespace Utility;
using namespace std;
//Random *Chromosome::generator = NULL;



Chromosome::Chromosome(LocusArray *loci, double pLambda, RecDistType *recombIndexLookup) : 
						loci(loci),	poissonLambda(pLambda), recombIndexLookup(recombIndexLookup), 
						xoEvents(NULL), hasCompleteGenotypes(false), phasedChrom(NULL), unphasedChrom(NULL), phase(false)	{ 	
	generator = &Utility::Random::globalGenerator;
	chrom.resize(loci->size()); 
}

Chromosome::Chromosome() : loci(NULL), poissonLambda(0.0), recombIndexLookup(NULL), 
						xoEvents(NULL), hasCompleteGenotypes(false), phasedChrom(NULL), unphasedChrom(NULL), phase(false) {
	generator = &Utility::Random::globalGenerator;
}

Chromosome::Chromosome(const Chromosome& other) : 
			generator(other.generator), chrom(other.chrom), loci(other.loci), 
			poissonLambda(other.poissonLambda), recombIndexLookup(other.recombIndexLookup), hasCompleteGenotypes(other.hasCompleteGenotypes),
			genotypePhase(other.genotypePhase), 
			phasedChrom(other.phasedChrom), 
			unphasedChrom(other.unphasedChrom),	
			phase(false), modelOnlyChrom(other.modelOnlyChrom)  {
	if (other.xoEvents) {
		xoEvents = new vector<size_t>;
		*xoEvents = *(other.xoEvents);
	}
	else 
		xoEvents = NULL;
	alleleSource = other.alleleSource;
	
}

void Chromosome::Distort(const double& frequency) {
	if (frequency > 0.0) {
	
		//How many snps to distort
		int snpCount = (int)((double)loci->size() * frequency);
		
		for (int i=0; i<snpCount; i++) {
			size_t idx = Random::globalGenerator.lrand(0, loci->size() - 1);
			chrom[idx]=!chrom[idx];
		}
	}	
}


Chromosome::Chromosome(const Chromosome &mother, const Chromosome &father) :  
			generator(mother.generator), chrom(mother.chrom), loci(mother.loci), 
			poissonLambda(mother.poissonLambda), 
			xoEvents(NULL), hasCompleteGenotypes(false), phase(false) {

	
	xoEvents = mother.BuildXOEvents();
	recombIndexLookup = mother.recombIndexLookup;
	size_t eventCount = xoEvents->size();
	size_t phaseCounts[] = {0,0};
	size_t curCount = 0, lastCount = 0;
	phase = false;
	//Count up how many cross overs we'll have for either In or Out of phase (F/T)
	for (size_t i=0; i<eventCount; i++) {
		curCount = (*xoEvents)[i];
		phaseCounts[phase]+=curCount - lastCount;
		lastCount = curCount;
		phase = !phase;
	}
	phaseCounts[phase]+=(loci->size() - lastCount);

	//Determine if we are supposed to start with M (false) or D (true)
	this->phase = generator->drand() < 0.5;

	//OK, if we get the wrong phase, we want to flip them mom/dad and the phase
	//This makes sure we start with the right one, but we still copy the least using
	//the [] operator
	if ((this->phase && phaseCounts[0] >= phaseCounts[1]) || (!this->phase && phaseCounts[0] < phaseCounts[1])) {
		chrom = father.chrom;
		phasedChrom = (Chromosome *)&mother;
		unphasedChrom = (Chromosome *)&father;
		this->phase = !this->phase;
	} else {
		chrom = mother.chrom;
		phasedChrom = (Chromosome *)&father;
		unphasedChrom = (Chromosome *)&mother;
	}
	hasCompleteGenotypes = phasedChrom->hasCompleteGenotypes;
	ExtractAlleleSource();
	bool phase = this->phase;
	genotypePhase.Insert(0, phase);
//cerr<<cout<<"Phase Count: "<<eventCount<<" ... -> ";
	//Set up the locations for events 
	for (size_t i=0; i<eventCount; i++ ){
		//Find the index 
		uint nextIdx = (*xoEvents)[i];
		//cout<<nextIdx<<" ";
		phase=!phase;
		genotypePhase.Insert(nextIdx+1, phase);
	}	
	genotypePhase.Insert(chrom.size(), !phase);
//cerr<<" | ";
//RenderPhase(cout, chrom.size());


	assert(eventCount != 0 || !this->phase);

}

void Chromosome::InitAlleleSource(const char *source) {
	assert(alleleSource.GetCount() == 0);
	alleleSource.Insert(0, source);
	alleleSource.Insert(LociCount(), "X");
}
int Chromosome::EvaluateKinship(Chromosome& other) {
	KinshipCalculator<string> k;
	return k.EvaluateKinship(alleleSource, other.alleleSource);
}
void Chromosome::ExtractAlleleSource() {
	//Phased chromosome gets the zero point
	bool phase = this->phase;
	int curPosition = 0;
	
	if (phasedChrom->alleleSource.GetCount() == 0 || unphasedChrom->alleleSource.GetCount()== 0)
		return;
	//To start with, we need to get the first node from the phased chromosome (after this, we alternate)
	AllSrcNodeType *node = NULL;
	int eventCount = xoEvents->size();

	if (eventCount > 0) {
		//Next, we iterate over the events. 
		for (int i=0; i<eventCount; i++ ){
			if (phase)
				node = phasedChrom->alleleSource.FindNearestMin(curPosition);
			else
				node = unphasedChrom->alleleSource.FindNearestMin(curPosition);
	
			int nextPosition = (*xoEvents)[i];
			//Record the node's string + the current XO point
			while (node && curPosition < nextPosition)  {
				//We don't want duplicates...that would be confusing
				AllSrcNodeType *oldNode = alleleSource.Find(curPosition);
				if (oldNode)
					alleleSource.Delete(oldNode);
					
				alleleSource.Add(curPosition, node->GetData());
				node = node->GetNext();
				if (node)
					curPosition = node->GetKey();
			}
			curPosition = nextPosition;
			//OK, now we switch to the other phase
			phase = !phase;
		}
		//We have to capture the tail of the other chromosome
		if (phase)
			node = phasedChrom->alleleSource.FindNearestMin(curPosition);
		else
			node = unphasedChrom->alleleSource.FindNearestMin(curPosition);
		
		while (node) {
			//We don't want duplicates...that would be confusing
			AllSrcNodeType *oldNode = alleleSource.Find(curPosition);
			if (oldNode)
				alleleSource.Delete(oldNode);

			alleleSource.Add(curPosition, node->GetData());
			node = node->GetNext();
			if (node)
				curPosition = node->GetKey();
			
		}
	}
	else
		alleleSource.Add(curPosition, phasedChrom->alleleSource.GetFirst()->GetData());
	if (alleleSource.Find(phasedChrom->LociCount())==NULL)
		alleleSource.Add(phasedChrom->LociCount(), "X");
}

void Chromosome::WriteXOPoints(std::ostream& os) {
	AllSrcNodeType *node = alleleSource.GetFirst();
	while (node) {
		os<<"\t"<<node->GetKey()<<"->"<<node->GetData()<<"\t";	
		node = node->GetNext();
	}
}


void Chromosome::ReferenceDistance(vector<double>& distances, int diseaseLocus) {
	LocusArray::iterator locItr = loci->begin();
	LocusArray::iterator locEnd = loci->end();

	Locus disease = loci->at(diseaseLocus);
	double diseaseDistance = disease.MapPosition();
	int idx = 0;
	while (locItr != locEnd) {
		distances[idx++]= fabs(locItr->MapPosition()-diseaseDistance);
		locItr++;
	}
}

void Chromosome::EvaluateRecombinants(vector<int>& r, int diseaseLocus) {
	RBTreeNode<size_t, bool> *xoBefore = genotypePhase.FindNearestMin(diseaseLocus);
	RBTreeNode<size_t, bool> *xoAfter = genotypePhase.FindNearestMax(diseaseLocus);
	xoAfter = xoBefore->GetNext();

	
	size_t xoPoint = xoBefore->GetKey();
//	if (xoPoint != 0 && xoPoint != diseaseLocus)
//		xoPoint++;

	for (size_t i=0; i<xoPoint; i++) {
//cerr<<"1 ";
		r[i]++;
	}

//for (int i=xoPoint; i<diseaseLocus; i++)
//cerr<<"0 ";
//cerr<<"[0] ";
	size_t locusCount = loci->size();
	xoPoint = xoAfter->GetKey();
	if (xoPoint != locusCount && xoPoint == (size_t)diseaseLocus)
		xoPoint++;
//for (int i=diseaseLocus+1; i<xoPoint; i++) 
//cerr<<"0 ";
	for (size_t i=xoPoint; i<locusCount; i++) {
//cerr<<"1 ";
		r[i]++;
	}

//cerr<<"\t\t";
}
void Chromosome::XOCount(vector<int>& counts, int diseaseLocus) {
	if (genotypePhase.GetCount() == 0)
		return;
	LocusArray::iterator locItr = loci->begin();
	LocusArray::iterator locEnd = loci->end();

	Locus disease = loci->at(diseaseLocus);
	RBTreeNode<size_t, bool> *diseaseNode = genotypePhase.FindNearestMax(diseaseLocus);
	int idx = 0;
	while (locItr != locEnd) {
		RBTreeNode<size_t, bool> *start = genotypePhase.FindNearestMin(idx);
		RBTreeNode<size_t, bool> *end = diseaseNode;

		int xoCount = 0;
		if (diseaseLocus == idx) {
			start = diseaseNode;
		} else if (diseaseLocus < idx) {
			end = start;
			start = diseaseNode;
		}
		
		while (start && start->GetKey() != end->GetKey()) {
			xoCount++;
			start = start->GetNext();
		}
		counts[idx++]+= xoCount;
		locItr++;
	}
}



vector<size_t> *Chromosome::BuildXOEvents() const {
	float start = (*loci)[0].MapPosition();
	float length = (*loci)[loci->size()-1].MapPosition();
	size_t eventCount = PoissonEventCount();

	vector<size_t> *events = new vector<size_t>;
	for (size_t i=0; i<eventCount; i++) {
		float loc = start + (*generator)(length);
		RecDistNodeType *node = recombIndexLookup->FindNearestMin(loc);
		size_t idx = 0;
		assert(node);
		idx = node->GetData();
		//If this index exists, let's get rid of it. two cross overs between two snps 
		//each other out
		vector<size_t>::iterator other = find(events->begin(), events->end(), idx);
		if (other != events->end()) 
			events->erase(other);
		else
			events->push_back(idx);
	}
	sort(events->begin(), events->end());	
	return events;
}

Chromosome &Chromosome::operator=(MinChrom &other ){
	MinChrom::Common *commonData = other.GetCommonData();

	loci 			= commonData->underlyingLoci;
	poissonLambda 	= commonData->poissonLambda;
	recombIndexLookup = commonData->recombIndexLookup;

	chrom.reset();
	chrom.resize(loci->size());

	if (xoEvents) {
		delete xoEvents;
		xoEvents = NULL;
	}
	xoEvents = NULL;
	hasCompleteGenotypes = false;
	phasedChrom = NULL;
	unphasedChrom = NULL;
	

	uint modelSize = other.GetModelSize();
	for (uint i=0; i<modelSize; i++) {
		uint idx = other.GetModelIndex(i);
		modelOnlyChrom[i] = other.At(idx);
		chrom[idx] = other.At(idx);
	}
//	alleleSource=other.alleleSource;

	return *this;
}

Chromosome &Chromosome::operator=(const Chromosome& other) {
	generator = other.generator;
	loci = other.loci;
	chrom = other.chrom;
	poissonLambda = other.poissonLambda;
	phase = other.phase;
	generator=other.generator;
	recombIndexLookup = other.recombIndexLookup;
	genotypePhase = other.genotypePhase;
	if (xoEvents) {
		delete xoEvents;
		xoEvents = NULL;
	}
	if (other.xoEvents) {
		xoEvents = new vector<size_t>;		
		*xoEvents = *(other.xoEvents);
	}
	hasCompleteGenotypes = other.hasCompleteGenotypes;
	phasedChrom = other.phasedChrom;
	unphasedChrom = other.unphasedChrom;

	modelOnlyChrom = other.modelOnlyChrom;
	if (alleleSource.GetCount() < other.alleleSource.GetCount()) 
		alleleSource=other.alleleSource;

	return *this;		
}    

/**
 * @brief Simply invert the alleles at each locus
 */
void Chromosome::Invert() {
	chrom.flip();
}

void Chromosome::ResolveTo(const Chromosome &other) {
	generator = other.generator;
	loci = other.loci;
	poissonLambda = other.poissonLambda;
	recombIndexLookup = other.recombIndexLookup;
	//Basically, we just want to get the bitset from the other one 
	chrom = other.chrom;
}

void Chromosome::ShowGenotypes( uint count, const char div ){
	if (count > chrom.size())
		count = chrom.size();
	cout<<" "<<div;
	if (xoEvents)
		cout<<"["<<xoEvents->size()<<"] ";
	else
		cout<<" "<<chrom.size()<<"  ";
	for (uint i=0; i<count; i++) 
		cout<<(*this)[i]<<" ";

	if (count < chrom.size()) {
		cout<<"...";
		for (uint i=chrom.size() - count; i<chrom.size(); i++) 
			cout<<(*this)[i]<<" ";
	}
}



void Chromosome::WriteBinary( std::ostream& file) {
	//We want 1 extra bit for the status
	//ShowGenotypes(64, '>');
	//cout<<"\n";

	assert(hasCompleteGenotypes);

	uint totalBlocks = (chrom.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
	dynamic_bitset<>::block_type raw[totalBlocks];

	//Convert the bitset to raw data
	to_block_range(chrom, raw);

	//cout<<"writing "<<totalBlocks<<" ints to the stream "<<(sizeof(dynamic_bitset<>::block_type)*totalBlocks)<<"\n";
	file.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	

}

void Chromosome::ReadBinary(std::istream &file, size_t lociCount, vector<uint> *modelLoci, bool retainChrom) {	
	uint totalBlocks = (lociCount / boost::dynamic_bitset<>::bits_per_block) + 1;
	dynamic_bitset<>::block_type raw[totalBlocks];


	boost::dynamic_bitset<> chrom;

	assert(!file.eof());
	file.read((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));

	chrom=dynamic_bitset<>(&raw[0], &raw[totalBlocks]);
	chrom.resize(lociCount);

	if (modelLoci) {
		size_t modelSize = modelLoci->size();

		modelOnlyChrom.clear();
		for (size_t i=0; i<modelSize; i++) 
			modelOnlyChrom[i] = chrom[(*modelLoci)[i]];
		
	}
	if (retainChrom)
		this->chrom = chrom;
	hasCompleteGenotypes = retainChrom;
	
	//ShowGenotypes( 64, '<');
}



// Use: Initializes specified locus according to allele
//      frequencies
// Arg: locusIndex -- index of locus to initialize
// ret: none
void Chromosome::InitLocus(int locusIndex){
	assert(generator);
	if(generator->drand()<=(*loci)[locusIndex].GetMinAlleleFreq())
    	chrom[locusIndex]=0;
	else
    	chrom[locusIndex]=1;
}


// Use: Initializes all loci in chromosome according to allele
//      frequencies
// Arg: none
// ret: none
void Chromosome::InitLoci(LocusAssociationGrid *grid){
	unsigned int totalLoci = Chromosome::loci->size();

	if (grid) {
		chrom = grid->GenerateChromosome(*generator);
	}
	else {
		for(unsigned int currLocus=0; currLocus < totalLoci; currLocus++){
			if(generator->drand()<=(*loci)[currLocus].Freq1())
				chrom[currLocus]=0;
			else
				chrom[currLocus]=1;
		}
	}

	hasCompleteGenotypes = true;
}

void Chromosome::InitLoci(vector<int>& alleles) {
	uint count = alleles.size();
	assert(count == loci->size());
	for (uint i=0; i<count; i++) {
		chrom[i] = alleles[i] == 1;
	}
	hasCompleteGenotypes = true;
}

/*
unsigned long random_poisson(double lamda)
{
       double p=exp(-lamda);
       double g=p;
       double u=Utility::Random::globalGenerator.drand();
       unsigned long k=0;
       while (u>g)             {
               p*=(lamda/(double)(++k));
               g+=p;
       }
       return k;
};
*/

size_t Chromosome::PoissonEventCount() const {
	double p=exp(-poissonLambda), 
		g=p, 
		u=generator->drand();

	size_t k=0;
	while (u>g)		{
		p*=(poissonLambda/(double)(++k));
		g+=p;
	}
	return k;
};


// Use: Produces new chromosome by crossing this chromosome with
//      a second one.  
// Arg: secondChrom - chromosome to cross
// Ret: Pointer to new chromosome
Chromosome Chromosome::Cross(Chromosome &secondChrom, size_t &crossoverEvents){
	assert(loci);
	Chromosome recombinant(loci, poissonLambda, recombIndexLookup);
assert(0);	
	Chromosome * currChrom = this;
	Chromosome * otherChrom = &secondChrom;
	Chromosome * tempChrom;
	
	// randomly choose a chromosome to start copying from
	if(generator->drand() < 0.5){
		currChrom = &secondChrom;
		otherChrom = this;
  	}
  


	crossoverEvents = 0;
	recombinant.hasCompleteGenotypes = true;

	// at each locus check for crossover event and
	// switch which chromosome is being used to copy
	// the alleles
	size_t totalLoci = LociCount();
	size_t lastLocus = totalLoci-1;
	assert(totalLoci > 0);
	for(size_t currLoc=0; currLoc<lastLocus; currLoc++){
	    recombinant[currLoc] = (*currChrom)[currLoc];
    	if(generator->drand() < (double)(*loci)[currLoc+1].MapDistance()){
			crossoverEvents++;
			tempChrom = currChrom;
			currChrom = otherChrom;
			otherChrom = tempChrom;
    	}  
  	}

	//cout<<"XOVR\tLoci Count\t"<<totalLoci<<"\tMBase Count\t"<<((*loci)[lastLocus].GetLocation())<<"\tX Count\t"<<crossoverEvents<<"\n";
	// set final locus, no need to check for crossover
	recombinant[lastLocus] = (*currChrom)[lastLocus];
	return recombinant;
}


void Chromosome::ResetToMinimal(vector<uint> &modelLoci) {
	uint modelSize = modelLoci.size();
	//modelOnlyChrom.resize(modelSize, false);
	modelOnlyChrom.clear();
	for (uint i=0; i<modelSize; i++) 
		modelOnlyChrom[i]=chrom[modelLoci[i]];
	chrom.resize(0);		
}
void Chromosome::ResolveGenotypes() {

	//We need to retrieve the data from the parents if we are working from an 
	//unloaded chromosome
//	if (modelOnlyChrom.size() > 0 && chrom.size() == 0)
//		chrom = unphasedChrom->chrom;
/*	cout<<"Resolve Genotypes: ";
	uint chromSize = chrom.size();
	if (phasedChrom) {
		//stringstream phasedPieces, unphasedPieces;
		for (uint i=0; i<chromSize; i++) 
			cout<<unphasedChrom->chrom[i]<<" ";
		cout<<"\t-\t";
		for (uint i=0; i<chromSize; i++) 
			cout<<phasedChrom->chrom[i]<<" ";
		cout<<"\t-\t";
	}
	for (uint i=0; i<chromSize; i++) 
		cout<<chrom[i]<<" ";
	cout<<"\t-\t"; 
*/


	if (xoEvents) {
		chrom = unphasedChrom->chrom;
//		cout<<"\tCh) "<<chrom<<"\tPh) "<<phasedChrom->chrom<<" -> ";
		int ooPhaseCounts = 0;
		int inPhaseCounts = 0;

		size_t crossoverEvents = xoEvents->size();
		uint lastLocus = loci->size() - 1;
		uint lastIndexConsidered = 0;
		hasCompleteGenotypes = true;
		bool phase = this->phase;
		//Set up the locations for events 
		for (size_t i=0; i<crossoverEvents; i++ ){
			//Find the index 
			uint nextIdx = (*xoEvents)[i] + 1;
			if (phase) {
//				cout<<"Processing partial: xoEvents="<<crossoverEvents<<". Remaining SNPS: "<<lastIndexConsidered<<" / "<<lastLocus<<"\n";
				//Nasty brute force approach. We need to do it in chunks
				for (; lastIndexConsidered<nextIdx; lastIndexConsidered++)  {
					ooPhaseCounts++;
					chrom[lastIndexConsidered] = phasedChrom->chrom[lastIndexConsidered];
	//				genotypePhase[lastIndexConsidered] =  true;
				}
			}else 	{
				inPhaseCounts += (nextIdx - lastIndexConsidered);
//				cout<<"Skipping partial: xoEvents="<<crossoverEvents<<". Remaining SNPS: "<<nextIdx<<" / "<<lastLocus<<"\n";

				lastIndexConsidered = nextIdx;
			}
			phase=!phase;
		}
	
		if (phase) {
//			cout<<"Finishing off resolution. xoEvents = "<<crossoverEvents<<". Remaining SNPS: "<<lastLocus - lastIndexConsidered<<" / "<<lastLocus<<"\n";
			while (lastIndexConsidered <= lastLocus) {
				ooPhaseCounts++;
				chrom[lastIndexConsidered] = phasedChrom->chrom[lastIndexConsidered];
				lastIndexConsidered++;			
			}
		}
		else {
			inPhaseCounts += (1 + lastLocus - lastIndexConsidered);
//			cout<<"Skipping the resulution: xoEvents = "<<crossoverEvents<<". Remaining SNPS: "<<lastLocus - lastIndexConsidered<<" / "<<lastLocus<<"\n";
		}
//		cout<<"\tCM) "<<chrom<<"\t";
//		RenderPhase(cout, 15);

/*		cout<<"\t:"<<crossoverEvents<<" (";
		RenderPhase(cout, 10);
		cout<<")\t";
	*/
		if (inPhaseCounts + 1 < ooPhaseCounts) {
/*
			cout<<"Well, we are about to go away!:\n";
			cout<<"Event Count: "<<xoEvents->size()<<"\tPhase: "<<this->phase<<"\n";
			cout<<"XO Events: ";
			for (size_t i=0; i<crossoverEvents; i++ )
				cout<<(*xoEvents)[i]<<" ";
			cout<<"\n";

			cout<<"In Phase: "<<inPhaseCounts<<" vs "<<ooPhaseCounts<<"\n";
*/
		}
//		assert(inPhaseCounts >= ooPhaseCounts);
	}



	delete xoEvents;
	xoEvents = NULL;
	phasedChrom = NULL;
	unphasedChrom = NULL;
/*	for (uint i=0; i<chrom.size(); i++) 
		cout<<chrom[i]<<" ";
	cout<<"\n";
*/
}


Chromosome Chromosome::DLCross(Chromosome &secondChrom, size_t &crossoverEvents) {
	assert(loci);
//	assert(secondChrom.hasCompleteGenotypes || modelOnlyChrom.size() > 0);
	Chromosome child(*this, secondChrom);
	crossoverEvents = child.xoEvents->size();
	return child;
}

Chromosome Chromosome::CrossImmediate(Chromosome &secondChrom, size_t &crossoverEvents){
	assert(loci);
	assert(secondChrom.hasCompleteGenotypes || modelOnlyChrom.size() > 0);
	Chromosome child(*this, secondChrom);
	crossoverEvents = child.xoEvents->size();
	child.ResolveGenotypes();
	return child;
}

// Use: Produces new chromosome by crossing this chromosome with
//      a second one.  
// Arg: secondChrom - chromosome to cross
// Ret: Pointer to new chromosome
/*
Chromosome Chromosome::CrossImmediate(Chromosome &secondChrom, size_t &crossoverEvents, RecDistType &recombIndexLookup){
	assert(loci);
	crossoverEvents = PoissonEventCount();
	Chromosome *first = this, *last = &secondChrom;

	if(generator->drand() < 0.5){
		first = &secondChrom;
		last = this;
  	}

	//If there are no cross over events, we can just return the first chrosome 
	//and save ourselves some unnecessary checking
	if (crossoverEvents == 0) 
		return *first;


	// randomly choose a chromosome to start copying from
  
	size_t totalLoci = LociCount();
	size_t lastLocus = totalLoci-1;
	size_t start = (*loci)[0].GetLocation();
	size_t stop = (*loci)[lastLocus].GetLocation();
	size_t lastIndexConsidered = 0;
	bool phase = false;

//	cout<<"\n\nXO Count: "<<crossoverEvents<<"\tSize: "<<totalLoci<<"\n";
//	cout<<"P: "<<*first<<"\nQ: "<<*last<<"\n";

	size_t *events = new size_t[crossoverEvents];
	for (int i=0; i<crossoverEvents; i++) 
		events[(size_t)i] = Utility::Random::globalGenerator.lrand(start, stop);
	sort(events, events+crossoverEvents);

	Chromosome recombinant = *first;
	
	
	//Set up the locations for events 
	for (int i=0; i<crossoverEvents; i++ ){
		//Find the index 
		RecDistNodeType * node = recombIndexLookup.FindNearestMax(events[i]);
		size_t nextIdx = node->GetData();

		if (phase) {
			//Nasty brute force approach. We need to do it in chunks
			for (lastIndexConsidered; lastIndexConsidered<nextIdx; lastIndexConsidered++)  {
				recombinant[lastIndexConsidered] = (*last)[lastIndexConsidered];
			}
		}else 	
			lastIndexConsidered = nextIdx;
		phase=!phase;
	}

	if (phase)
		while (lastIndexConsidered < lastLocus) {
			recombinant[lastIndexConsidered] = (*last)[lastIndexConsidered];
			lastIndexConsidered++;
		}
		
	hasCompleteGenotypes = true;

	return recombinant;
}*/

// Use: Output chromsome using overloaded operator
// Arg: os -- output stream
//      chrom -- chromosome to output
// Ret: output stream
std::ostream & operator << (std::ostream & os, Chromosome & chrom){
	unsigned int numLoci = chrom.LociCount();
	//assert(chrom.hasCompleteGenotypes);
	for(unsigned int i=0; i<numLoci; i++)
    	os << chrom[i]+1 << " ";
  	os << std::endl;
  	return os;
}




void Chromosome::WritePedFormat(std::ostream& os, float threshold, uint first, uint last) {
	//std::string vals[] = { "1 1 ", "1 2 ", "2 2 " };
	unsigned int numLoci = LociCount();
	if (last == 0 || last > numLoci)
		last = numLoci;

	for(unsigned int i=first; i<last; i++) {
		Locus &l = (*loci)[i];
		if (l.GetMinAlleleFreq() >= threshold) {
	    	os << chrom[i]+1<< " ";
		}
	}
  	os << std::endl;
}

void Chromosome::ExpandModelLoci(vector<uint> &modelLoci) {
	if (hasCompleteGenotypes) {
		ResolveGenotypes();
	} else {
		chrom.reset();
		chrom.resize(loci->size(), false);
		//cout<<"Expanding Model Loci: ";
		for (uint i=0; i<modelLoci.size(); i++) {
			if (phasedChrom) {
				uint locusIndex = modelLoci[i];
				//We need to figure out which phase we the genotype is at and ask our parent chromosome
				if (GetPhase(locusIndex))
					chrom[locusIndex] = (*phasedChrom)[locusIndex];
				else
					chrom[locusIndex] = (*unphasedChrom)[locusIndex];
			}
			else {
				chrom[modelLoci[i]] = modelOnlyChrom[i];
			}
		}
	}
		
}

void Chromosome::SetValue(uint locus, uint value) {
	assert(chrom.size() > locus);
	chrom[locus] = (value==2);
}


#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(ChromosomeTest);

ChromosomeTest::ChromosomeTest() :  ch1(NULL), ch2(NULL), ch3(NULL), ch4(NULL) { }


ChromosomeTest::~ChromosomeTest(){	
	if (ch1)
		delete ch1;
	if (ch2)
		delete ch2;
	if (ch3)
		delete ch3;
	if (ch4)
		delete ch4;
	if (ch5)
		delete ch5;
	if (ch6)
		delete ch6;
	if (ch7)
		delete ch7;
	if (ch8)
		delete ch8;

	if (ex1)
		delete ex1;
	if (ex2)
		delete ex2;
	if (ex3)
		delete ex3;
	if (ex4)
		delete ex4;

}

void ChromosomeTest::tearDown() {
}

void ChromosomeTest::setUp() { 
	RecDistType recombIndexLookup;
	int locusCount = 1000;
	for (int i=0; i<locusCount; i++) {
		loci.push_back(Locus(0.0001, 1, i, i*10));
		recombIndexLookup.Insert(i*10, i);
	}

	ch1 = new Chromosome(&loci, 5, &recombIndexLookup);
	ch2 = new Chromosome(&loci, 5, &recombIndexLookup);
	ex1	= new Chromosome(&loci, 5, &recombIndexLookup);
	ex2 = new Chromosome(&loci, 5, &recombIndexLookup);
	ex3 = new Chromosome(&loci, 5, &recombIndexLookup);
	ex4 = new Chromosome(&loci, 5, &recombIndexLookup);

	

	for (int i=0; i<locusCount; i++) {
		ch1->SetValue(i, 0);
		ch2->SetValue(i, 1);
		ex1->SetValue(i, 0);
		ex2->SetValue(i, 1);
		ex3->SetValue(i, 0);
		ex4->SetValue(i, 1);
	}
	ch1->InitAlleleSource("CH1");
	ch2->InitAlleleSource("CH2");
	ex1->InitAlleleSource("EX1");
	ex2->InitAlleleSource("EX2");
	ex3->InitAlleleSource("EX3");
	ex4->InitAlleleSource("EX4");
	ch3 = new Chromosome(*ch1, *ch2);
	ch4 = new Chromosome(*ch1, *ch2);

	ch5 = new Chromosome(*ch3, *ex1);
	ch6 = new Chromosome(*ch4, *ex2);
	
	ch7 = new Chromosome(*ch5, *ex3);
	ch8 = new Chromosome(*ch6, *ex4);


	cout<<"\nCH1 -";
	ch1->WriteXOPoints(cout);
	cout<<"\nCH2 -";
	ch2->WriteXOPoints(cout);
	cout<<"\nCH3 -";
	ch3->WriteXOPoints(cout);
	cout<<"\nCH4 -";
	ch4->WriteXOPoints(cout);
	cout<<"\nCH5 -";
	ch5->WriteXOPoints(cout);
	cout<<"\nCH6 -";
	ch6->WriteXOPoints(cout);
	cout<<"\nCH7 -";
	ch7->WriteXOPoints(cout);
	cout<<"\nCH8 -";
	ch8->WriteXOPoints(cout);
	cout<<"\n";
}

void ChromosomeTest::TestKinship() {
	//Test kinship with self
	float locusCount = 1000.0;
	float kinship = (float)ch1->EvaluateKinship(*ch1) / locusCount;
cout<<"Kinship (Self): "<<kinship<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (Self)", 1.0, kinship, 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (none)", 0.0, (float)ch1->EvaluateKinship(*ch2)/ locusCount, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (parent child)", 0.5, (float)ch1->EvaluateKinship(*ch4)/ locusCount, 0.15);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (sibling)", 0.5, (float)ch3->EvaluateKinship(*ch4)/ locusCount, 0.15);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (gran/child)", 0.25, (float)ch2->EvaluateKinship(*ch5)/ locusCount, 0.15);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (gg/child)", 0.125, (float)ch2->EvaluateKinship(*ch7)/ locusCount, 0.15);
}
#endif


}
