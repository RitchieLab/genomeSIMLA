//
// C++ Implementation: ldplotter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldplotter.h"
#include "ldprimepngwriter.h"
#include <iomanip>
#include <boost/timer.hpp>
#include <sstream>
#include "chrompool.h"
#include "poolmanager.h"
#include "ldlineplot.h"
#include "utility/array2d.h"

namespace Simulation {
namespace Visualization {

uint LdPlotter::maxSnpDistance = 500000;

struct SortByLocation {
	bool operator() (const Locus& A, const Locus& B) const {
		return A.GetLocation() < B.GetLocation();
	}
};



//We are assuming that the allele frequencies are already calculated prior to establishing the plotter object
LdPlotter::LdPlotter(float thresh, vector<Chromosome>& pool, vector<Locus>& loci) : pool(pool), threshold(thresh), ldStatistics(NULL), ldStart(NULL) {

	vector<Locus>::iterator i = loci.begin();
	vector<Locus>::iterator end = loci.end();
	for (; i != end; i++)  
		if (i->GetMinAlleleFreq() > thresh) 
			this->loci.push_back(*i);

	sort(this->loci.begin(), this->loci.end(), SortByLocation() );

}


LdPlotter::~LdPlotter() {
	//We don't want to delete pool nor loci since they are just copies of the real data

	if (ldStatistics)
		delete[] ldStatistics;
	if (ldStart)
		delete[] ldStart;

}

//vector<HaplotypeBlock *> *LdPlotter::Init() {
BlockList *LdPlotter::Init() {
	//First, we need to count the deepest row count
	haplotypes = haplotypeHead.GetValidBlocks();
	sort(haplotypes->begin(), haplotypes->end(), BlockSizeEval());
	return haplotypes;
}


void LdPlotter::CountRC(size_t first, size_t last, size_t &rows, size_t &cols) {
	rows = cols = 0;

	size_t validFirst = 0, validLast = 0;

	if (last >= loci.size()) {
		cout<<"Trying to count "<<last<<" rows even though there are only "<<loci.size()<<" available\n";
		cout<<"OK, we can't search past the end of the locus array\n";
		assert(0);
	}

	for (size_t p=first; p<=last; p++) {
		if (validFirst == 0)
			validFirst = p;
		validLast = p+1;
		cols++;
		size_t q=p;
		uint curPos = loci[p].GetLocation();
		
		while (++q<=last && loci[q].GetLocation() - curPos < maxSnpDistance);

		if (q - p> rows)
			rows=q - p;
	}

//	cout<<"\tFirst:    "<<setw(10)<<first<<"\tLast:   "<<setw(10)<<last<<"\n";
//	cout<<"\trows:     "<<setw(10)<<rows<<"\tCols:   "<<setw(10)<<cols<<"\n";
	if (rows > cols) {
		rows = cols;
		cout<<"Rows adjusted to match columns\n";
	}
}
void LdPlotter::CountRC(size_t &rows, size_t &cols) {
	size_t lociCount = loci.size() - 1;
	size_t start = 0;
	CountRC(start, lociCount, rows, cols);
}

/**
 * @brief This is used by the plotter to calculate multilocus genotype frequencies
 */
struct CountFreqArgs {
	uint A;				///<First locus of the pair
	uint B1; 			///<Starting range of the second locus
	uint B2;			///<The end of the range of the second locus
	Probabilities *p;	///<The probabilities associated with these ml genotypes
	LdPlotter *plotter;	
	uint firstInd;		///<where to start in the chromosome pool
	uint count;			///<How many chromosomes to look at

	//void Increment() { A++; B1++; p=NULL; }
	~CountFreqArgs() { delete[] p; }
	CountFreqArgs(uint A, uint B1, uint B2, LdPlotter *plotter, uint firstInd, uint count) : A(A), B1(B1), B2(B2), p(NULL), plotter(plotter), firstInd(firstInd), count(count) { }
};


void *LdPlotter::CountFrequencies(void *args) {
//	boost::timer progress;
	CountFreqArgs *thArg = (CountFreqArgs*)args;
	thArg->plotter->CountFrequencies(thArg->A, thArg->B1, thArg->B2, thArg->p, thArg->firstInd, thArg->count);
//	cout<<"Finished with thread: "<<progress.elapsed()<<"\n";
	return args;
}


Probabilities *LdPlotter::ThreadedFrequencyCount(uint A, uint B1, uint B2, vector<pthread_t>& threads) {
	uint count = B2-B1 + 1;
	if (count == 0)
		return NULL;
	uint threadCount = threads.size();
	uint totalIndividuals = pool.size();
	if (ChromPool::fastLD && totalIndividuals > ChromPool::maxLDIndividuals) {
		totalIndividuals = ChromPool::maxLDIndividuals;
	}
	uint portion = totalIndividuals / (threadCount+1);

	//We want to process some with the main thread.
	//threadCount--;
	
	//pthread_t threads[threadCount];
	//launch the individual threads
	
	uint firstIndividual = 0;
	for (uint i =0; i<threadCount; i++) {
		CountFreqArgs *args = new CountFreqArgs(A, B1, B2, this, firstIndividual, portion); 
		pthread_create(&threads[i], NULL, &CountFrequencies, (void*)args);
		firstIndividual += portion;
	}
	Probabilities *p = NULL; 
	CountFrequencies(A, B1, B2, p, firstIndividual, (totalIndividuals - firstIndividual));

 	
	for (uint th=0; th<threadCount; th++) {
		void *rtn;
		pthread_join(threads[th], &rtn);
		CountFreqArgs *args = (CountFreqArgs*)rtn;
			
		for (uint i=0; i<count; i++) {
			assert(p[i].P == args->p[i].P);
			assert(p[i].Q == args->p[i].Q);
			p[i].cAB+=args->p[i].cAB;
			p[i].cAb+=args->p[i].cAb;
			p[i].caB+=args->p[i].caB;
			p[i].cab+=args->p[i].cab;	
		}
		delete args;
	}

	for (uint i=0; i<count; i++) {
		p[i].CalcSum();
		if (((uint)p[i].sum) != totalIndividuals && (uint)p[i].sum != 0) {
			cout<<"We have a problem with "<<A<<", "<<B1<<"-"<<B2<<" ("<<i<<") "<<(int)(p[i].sum)<<" != "<<totalIndividuals<<"\n";
			p[i].Report(cout);
		}
		assert(((uint)(p[i].sum)) == totalIndividuals || (uint)p[i].sum == 0);
	}


	return p;
}

void LdPlotter::CountFrequencies(uint A, uint B1, uint B2, Probabilities* &p, uint firstIndividual, uint totalInd) {
	uint indCount = firstIndividual + totalInd;
	if (indCount > pool.size())
		indCount = pool.size();
	//locix2 array for each allele

	if (p) {
		cout<<"Deleting P\n";
		delete[] p;	
	}
	p = NULL;


	if (A == B2)
		return;


	uint count = B2-B1 + 1;
	p = new Probabilities[count];

	Utility::Array2D<uint> freq(count, 4);
	//uint freq[count][4];

	//memset((void*)freq, sizeof(uint)*(count*4), 0);
	/**for (uint i = 0; i<count; i++) 
		for (uint j=0; j<4; j++) 
			freq[i][j]=0;
	 */
	uint maxLoc = loci[A].GetLocation() + maxSnpDistance;
	size_t idA = loci[A].GetID();
	//os << setw(width) << right << "Loc  "<< setw(width) << " All1 "<< setw(width) << "All2 "<< setw(width) << "Recomb   " << std::endl;

//Just testing something
	if (ChromPool::fastLD && indCount > ChromPool::maxLDIndividuals) {
		indCount = ChromPool::maxLDIndividuals;
	}
	
//	cout<<"Counting multilocus genotypes for "<<firstIndividual<<" to "<<indCount<<"\n";

	for (uint ind = firstIndividual; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];
		bool doContinue = true;
		uint localP=0;
		for (uint q=B1; doContinue && q<=B2; q++) {
			Locus &loc = loci[q];
			if (doContinue = loc.GetLocation() < maxLoc)			
				//Basically, A is 0, a is 2   and B is +0 while b is +1
				freq(localP++,curChrom[idA]*2+curChrom[loc.GetID()])++;
				//freq[q-B1][curChrom.At(A)*2+curChrom.At(q)]++;
		}
	}


	for (uint i=0; i<count; i++) {
		p[i].cAB=freq(i,0);
		p[i].cAb=freq(i,1);
		p[i].caB=freq(i,2);
		p[i].cab=freq(i,3);
		p[i].CalcSum();

		p[i].P = A;
		p[i].Q = i+B1;

/*		cout<<A<<"x"<<i+B1<<"\t "<<setw(8)<<p[i].cAB<<setprecision(4)<<setw(8)<<p[i].pAB()
			<<setw(8)<<p[i].cAb<<setprecision(4)<<setw(8)<<p[i].pAb()
			<<setw(8)<<p[i].caB<<setprecision(4)<<setw(8)<<p[i].paB()
			<<setw(8)<<p[i].cab<<setprecision(4)<<setw(8)<<p[i].pab()
			<<"\t"<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<setw(8)<<p[i].mafA<<setw(8)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(2)<<p[i].mafB<<setw(8)<<p[i].sum<<"\n";		
	*/	
		uint sumGenotypes = p[i].cAB+p[i].cAb+p[i].caB+p[i].cab;
		if (sumGenotypes != totalInd && sumGenotypes != 0) {
			cerr<<" An inconsistency was encountered when counting frequencies for "<<A<<"x"<<i + B1<<" ("<<i<<")\n";
			cerr<<"\t"<<p[i].cAB<<" + "<<+p[i].cAb<<" + "<<p[i].caB<<" + "<<p[i].cab<<" != "<<totalInd<<"\n";
		}
		assert(sumGenotypes == totalInd || sumGenotypes == 0);

	}
	

}

void LdPlotter::CountFrequencies(uint A, uint B1, uint B2, Probabilities* &p) {
	uint indCount = pool.size();
	//locix2 array for each allele

	if (p) {
		cout<<"Deleting P\n";
		delete[] p;	
	}
	p = NULL;


	if (A == B2)
		return;


	uint count = B2-B1 + 1;
	p = new Probabilities[count];

	uint freq[count][4];

	//memset((void*)freq, sizeof(uint)*(count*4), 0);
	for (uint i = 0; i<count; i++) 
		for (uint j=0; j<4; j++) 
			freq[i][j]=0;

	uint maxLoc = loci[A].GetLocation() + maxSnpDistance;
	size_t idA = loci[A].GetID();
	//os << setw(width) << right << "Loc  "<< setw(width) << " All1 "<< setw(width) << "All2 "<< setw(width) << "Recomb   " << std::endl;

//Just testing something
	if (ChromPool::fastLD && indCount > ChromPool::maxLDIndividuals) {
		indCount = ChromPool::maxLDIndividuals;
	}
	


	for (uint ind = 0; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];
		bool doContinue = true;
		uint localP=0;
		for (uint q=B1; doContinue && q<=B2; q++) {
			if (doContinue = loci[q].GetLocation() < maxLoc)			
				//Basically, A is 0, a is 2   and B is +0 while b is +1
				freq[localP++][curChrom[idA]*2+curChrom[loci[q].GetID()]]++;
				//freq[q-B1][curChrom.At(A)*2+curChrom.At(q)]++;
		}
	}


	for (uint i=0; i<count; i++) {
		p[i].cAB=freq[i][0];
		p[i].cAb=freq[i][1];
		p[i].caB=freq[i][2];
		p[i].cab=freq[i][3];
		p[i].CalcSum();

		p[i].P = A;
		p[i].Q = i+B1;

/*		cout<<A<<"x"<<i+B1<<"\t "<<setw(8)<<p[i].cAB<<setprecision(4)<<setw(8)<<p[i].pAB()
			<<setw(8)<<p[i].cAb<<setprecision(4)<<setw(8)<<p[i].pAb()
			<<setw(8)<<p[i].caB<<setprecision(4)<<setw(8)<<p[i].paB()
			<<setw(8)<<p[i].cab<<setprecision(4)<<setw(8)<<p[i].pab()
			<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<setw(8)<<p[i].mafA<<setw(8)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(2)<<p[i].mafB<<setw(8)<<p[i].sum<<"\n";		
		if (p[i].cAB+p[i].cAb+p[i].caB+p[i].cab != indCount) {
			cerr<<" An inconsistency was encountered when counting frequencies for "<<A<<"x"<<i + B1<<" ("<<i<<")\n";
			cerr<<"\t"<<p[i].cAB<<" + "<<+p[i].cAb<<" + "<<p[i].caB<<" + "<<p[i].cab<<" != "<<indCount<<"\n";
		}
		assert(p[i].cAB+p[i].cAb+p[i].caB+p[i].cab == indCount);
*/
	}
	

}
void LdPlotter::CountFrequencies(uint A, uint B, Probabilities &p) {
	uint indCount = pool.size();
	//locix2 array for each allele

	uint freq[4];

	for (uint i=0; i<4; i++) 
		freq[i]=0;

	//os << setw(width) << right << "Loc  "<< setw(width) << " All1 "<< setw(width) << "All2 "<< setw(width) << "Recomb   " << std::endl;
	for (uint ind = 0; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];
		//cout<<curChrom<<"\n";
		//Basically, A is 0, a is 2   and B is +0 while b is +1
		freq[curChrom.At(A)*2+curChrom.At(B)]++;
	}

	p.cAB=freq[0];
	p.cAb=freq[1];
	p.caB=freq[2];
	p.cab=freq[3];
	p.CalcSum();
//	cout<<setw(8)<<p.cAB<<setw(8)<<p.cAb<<setw(8)<<p.caB<<setw(8)<<p.cab<<setw(8)<<p.mafA<<setw(8)<<p.mafB<<setw(8)<<p.sum<<"\n";
	assert(p.cAB+p.cAb+p.caB+p.cab == indCount);
	
}
/*
uint LdPlotter::BuildHaplotypeBlocks() 	{
	uint chromCount = loci.size();
	uint P	 		= 0;
	uint Q			= 0;
	
	haplotypeHead.Reset(Q, chromCount);
	//HaplotypeBlock *currBlock = NULL;

	//Using the 4 gamete rule, we will inch through each neighboring snp to look for missing cells and call the contiguous sections blocks
	minBlockDensity = (size_t)-1;
	maxBlockDensity = 0;
	while (Q<chromCount) {
		P=Q;

		Locus &A=loci[P];

		Probabilities *probs = NULL;
		CountFrequencies(P, P+1, chromCount, probs);
		uint pidx = 0;
		if (++Q < chromCount) {
			Locus &B=loci[Q];

			//cout<<"MAF(A)="<<A.GetMinAlleleFreq()<<"  MAF(B)="<<B.GetMinAlleleFreq()<<"\n";

			Probabilities &prob = probs[pidx++];

			//Get the basic counts and determine if they should be part of a block
			//CountFrequencies(A.GetID(), B.GetID(), prob);

			prob.mafA=A.GetMinAlleleFreq();
			prob.mafB=B.GetMinAlleleFreq();


		}
		delete[] probs;
	}	
	

	//For reporting purposes, we'll sort them according to the number of snps found in each block
	//sort(haplotypes.begin(), haplotypes.end(), BlockSizeEval());

	return haplotypes->size();
}
*/
/**
 * @note This is leftover from an original approach and is no longer used
 */
void LdPlotter::CalculateAlleleFrequencies(vector<Locus> &loci) {
	uint indCount = pool.size();
	uint lociCount = loci.size();
	//locix2 array for each allele
	uint alleleCount = 2 * lociCount;
	//MMmmmm, MS won't bother to comply with C99, so I gotta go back to mid 90s 
	uint *freq = new uint[alleleCount];


	for (uint i=0; i<alleleCount; i++) {
		freq[i]=0;
	}

	//os << setw(width) << right << "Loc  "<< setw(width) << " All1 "<< setw(width) << "All2 "<< setw(width) << "Recomb   " << std::endl;
	for (uint ind = 0; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];
		for(uint i=0; i<lociCount; i++) {
//			cout<<"Setting "<<i*2+curChrom.At(loci[i].GetID())<<" Up one\n";
			freq[i*2+curChrom.At(loci[i].GetID())]++;
		}
	}

	uint idx = 0;
	for (uint i=0;i<alleleCount; ) {
		double al1 = (double)(freq[i++])/(double)(indCount);
		double al2 = (double)(freq[i++])/(double)(indCount);
//		cout<<"Assigning "<<al1<<" "<<al2<<" to locus: "<<idx<<"\n";
		loci[idx++].AssignFreq(al1, al2);
	}

	delete[] freq;
}

void LdPlotter::SetLabel(const char *lbl) {
	label = lbl;
}

//string LdPlotter::WriteLdDPrime(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl) {
string LdPlotter::WriteLdDPrime(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl) {
	char adjustedFilename[1024];
	sprintf(adjustedFilename, "%s%s-dp", filename, lbl);

	uint lociCount = loci.size();
	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	size_t first = 0;
	if (block)
		first = block->GetFirstIndex();
	

	if (first > ldBufferSize)
		first-=ldBufferSize;
	else
		first = 0;

	size_t last = lociCount - 1;
	if (block)
		last = block->GetLastIndex() + ldBufferSize;

	if (last >= lociCount)
		last = lociCount - 1;

	assert(last<lociCount);

	//Since we want to limit the size of the plot to heights that will actually be used, we need to 
	//scan over all possible couplings to determine who is within range for LD consideration
	size_t cols, rows;
	CountRC(first, last, rows, cols);
	

	string pngFilename = "";

	if (rows && cols) {
		//Debugging verification to make sure we have correctly sorted our loci
		if (loci[first].GetLocation() > loci[last-1].GetLocation()) {
			cout<<"We have a problem\n";
			vector<Locus>::iterator i = loci.begin();
			vector<Locus>::iterator e = loci.end();
			uint idx =0;
			for (; i!=e; i++) {
				cout<<"\t"<<idx++<<"\t"<<setw(20)<<i->GetLabel()<<setw(10)<<i->GetLocation()<<"\n";
			}
			abort();
		}

		//We'll set up the image parameters for the current region.
		ImageParameters *params = new ImageParameters;
		params->Init(cols, rows);

		LdWriter *png = new LdPrimePngWriter(params, loci[first].GetLocation(), loci[last].GetLocation());

		stringstream lbl;
		lbl<<label;
		if (block)
			lbl<<"("<<block->GetLabel()<<")";
		lbl<<": D'";
		png->SetLabel(lbl.str().c_str());

		//((LdDPrimePng *)png)->SetLabel(label.c_str());
		imgHeight = params->dimensions.y;
		imgWidth  = params->dimensions.x;

		//We'll write the LD values to file
		pngFilename = WriteLdData( adjustedFilename, first, last, png, params->snpDepth);


		//Then add each of the blocks- making the central snp more boldly marked
		haplotypes = haplotypeHead.GetValidBlocks();
		haplotypeHead.GetBlockDensities(minBlockDensity, maxBlockDensity);

		((LdPrimePngWriter *)png)->SetBlockDensity(minBlockDensity, maxBlockDensity);
		BlockList::iterator itr=haplotypes->begin();
		BlockList::iterator end=haplotypes->end();

		while (itr != end) {
			png->AddBlock(*(itr++), 1);
		}
		if (block)
			png->AddBlock( block, 3 );
		png->Close();
		delete png;
	}
	else {
		cout<<"Unable to write an empty png!\n";
		abort();
	}
	return pngFilename;
}
//string LdPlotter::WriteLdRSquared(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl) {
string LdPlotter::WriteLdRSquared(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl) {
	char adjustedFilename[1024];
	sprintf(adjustedFilename, "%s%s-rs", filename, lbl);

	uint lociCount = loci.size();

	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	size_t first = 0;
	if (block)
		first = block->GetFirstIndex();

	if (first > ldBufferSize)
		first-=ldBufferSize;
	else
		first = 0;

	size_t last = lociCount - 1;
	if (block)
		last = block->GetLastIndex() + ldBufferSize;

	if (last >= lociCount)
		last = lociCount - 1;

	assert(last<lociCount);

	//Calculate the number of snps to be considered as well as the largest pairing based on the max distance setting
	size_t cols, rows;
	CountRC(first, last, rows, cols);

	string pngFilename = "";

	if (rows && cols) {
		//Verify that we have a sorted array of loci
		if (loci[first].GetLocation() > loci[last-1].GetLocation()) {
			cout<<"We have a problem\n";
			vector<Locus>::iterator i = loci.begin();
			vector<Locus>::iterator e = loci.end();
			uint idx =0;
			for (; i!=e; i++) {
				cout<<"\t"<<idx++<<"\t"<<setw(20)<<i->GetLabel()<<setw(10)<<i->GetLocation()<<"\n";
			}
			abort();
		}

		//Create the image parameters and set up the renderer 
		ImageParameters *params = new ImageParameters;
		params->Init(cols, rows);
		LdWriter *png = new LdRSquaredPngWriter(params, loci[first].GetLocation(), loci[last].GetLocation());
		((LdRSquaredPngWriter *)png)->SetBlockDensity(minBlockDensity, maxBlockDensity);
		//((LdDPrimePng *)png)->SetLabel(label.c_str());
		imgHeight = params->dimensions.y;
		imgWidth  = params->dimensions.x;

		stringstream lbl;
		lbl<<label;
		if (block)
			lbl<<"("<<block->GetLabel()<<")";
		lbl<<": R^2";
		png->SetLabel(lbl.str().c_str());

		pngFilename = WriteLdData( adjustedFilename, first, last, png, params->snpDepth);

		//Iterate over all blocks and add them to the plot. If it's the central block, we'll make it darker
		BlockList::iterator itr=haplotypes->begin();
		BlockList::iterator end=haplotypes->end();
		while (itr != end) {
			png->AddBlock(*(itr++), 1);
		}
		if (block)
			png->AddBlock( block, 3 );
		png->Close();
		delete png;
	}
	else {
		cout<<"Unable to write an empty png!\n";
		abort();
	}
	return pngFilename;
}

string LdPlotter::WriteLdReport(const char *filename, BlockListNode*block) {
	LdWriter *txtReport = new LdTextReport();
	string r = WriteLdData(filename, block, txtReport, 0);
	delete txtReport;
	return r;
}

/**
 * @brief this function was left over from a previous approach, but it's unclear whether we'll need it
 */
string LdPlotter::WriteLdData(const char *filename, BlockListNode *block, LdWriter* writer, uint maxSnpDepth) {
	char adjustedFilename[1024];
	sprintf(adjustedFilename, "%s-LD-Report.txt", filename );

	uint lociCount = loci.size();

	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	size_t first = 0;
	if (block)
		first = block->GetFirstIndex();

/*	if (first > ldBufferSize)
		first-=ldBufferSize;
	else
		first = 0;
*/
	size_t last = lociCount - 1;
	if (block)
		last = block->GetLastIndex();// + ldBufferSize;

	if (last >= lociCount)
		last = lociCount - 1;

	assert(last<lociCount);

	//Calculate the number of snps to be considered as well as the largest pairing based on the max distance setting
	size_t cols, rows;
	CountRC(first, last, rows, cols);

	string reportFilename = "";
	if (rows && cols) {
		//Verify that we have a sorted array of loci
		if (loci[first].GetLocation() > loci[last-1].GetLocation()) {
			cout<<"We have a problem\n";
			vector<Locus>::iterator i = loci.begin();
			vector<Locus>::iterator e = loci.end();
			uint idx =0;
			for (; i!=e; i++) {
				cout<<"\t"<<idx++<<"\t"<<setw(20)<<i->GetLabel()<<setw(10)<<i->GetLocation()<<"\n";
			}
			abort();
		}

		reportFilename = WriteLdData(filename, first, last, writer, maxSnpDepth);		
					  
	}
	return reportFilename;
}


void LdPlotter::CalculateLD2(const char *filename) {
	boost::timer progress;

	uint first = 0;
	assert(loci.size() > 0);
	uint last = loci.size() - 1;
	size_t cols, rows;
	CountRC(first, last, rows, cols);


	LdLinePlot dpLP(NULL, maxSnpDistance);
	dpLP.SetYAxisLabel( (char *)"DPrime" );

	LdLinePlot rsLP(NULL, maxSnpDistance);
	rsLP.SetYAxisLabel( (char *)"RSquared" );

	uint snpWidth = last-first;

	haplotypeHead.Reset(first, last);
	//Initialize the LD cache
	ldStatistics = new LDStatistics[snpWidth*rows + 1];
	ldStart      = new size_t[snpWidth + 1];
	size_t idx = 0;
	bool withinDist = true;

	float log10 = 1/log(10.0);
	
	//pthread_t freqCalcThread;
	//CountFreqArgs args(first, first+1, last, this);

	//pthread_create(&freqCalcThread, NULL, CountFrequencies, (void*)&args);

	Probabilities *probabilities = NULL;

	double sumLocalProcessing = 0.0;
	vector<pthread_t> threads(PoolManager::threadsPerChrom - 1);
	//Walk over each snp in the selected region
	for (uint P=first; P<=last && ChromPool::continueRunning; P++) {
		probabilities = ThreadedFrequencyCount(P, P+1, last, threads);

		boost::timer localProcessing;

		withinDist = true;
		Locus &A=loci[P];


		float PdotA=A.Freq1();
		float Pdota=A.Freq2();

		//Add the Snp to the chart so it's label and MAF will be seen	
//		writer->WriteHeader( P, A );
		ldStart[P] = idx;
		//uint q1 = loci[P+1].GetID();
		//uint qLast = loci[last].GetID();		
		//From "optimization"
		

		
		//CountFrequencies(P, P+1, last, probabilities);
		uint qidx = 0;
		//Since we are looking at haplotypes, let's inch forward from P one step at a time
		for (uint Q=P+1; Q<=last && withinDist; Q++) {
			Locus &B=loci[Q];
			//Probabilities prob;
			Probabilities &prob = probabilities[qidx++];
			

			uint dist = abs((int)B.GetLocation() - (int)A.GetLocation());

			withinDist = dist < maxSnpDistance;

			//Make sure we don't waste time considering SNPs that are too far away from one another
			if (withinDist) {
				//CountFrequencies(A.GetID(), B.GetID(), prob);
				//float Pdotb=B.GetMinAlleleFreq();
				//float PdotB=1.0 - Pdotb;
				float PdotB=B.Freq1();
				float Pdotb=B.Freq2();

//Sanity check
				float marginal = prob.pAB() + prob.pAb();
				if (marginal < PdotA - 0.0001 || marginal > PdotA + 0.0001) {
					cout<<"Marginals aren't adding up for "<<A.GetID()<<"x"<<B.GetID()<<" "<<P<<"x"<<Q<<":\n";
					cout<<"P(A): "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<PdotA<<" != P(AB) "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<prob.pAB()<<" + P(Ab) "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<prob.pAb()<<"\n";
					cout<<"\nP(A): "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<PdotA<<" P(a): "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<Pdota<<", P(B) "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<PdotB<<" P(b) "<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
				marginal = prob.paB() + prob.pAB();
				if (marginal < PdotB - 0.0001 || marginal > PdotB + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(B): "<<PdotB<<" != P(aB) "<<prob.paB()<<" + P(AB) "<<prob.pAB()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
				marginal = prob.paB() + prob.pab();
				if (marginal < Pdota - 0.0001 || marginal > Pdota + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(a): "<<Pdota<<" != P(aB) "<<prob.paB()<<" + P(ab) "<<prob.pab()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}	
				marginal = prob.pAb() + prob.pab();
				if (marginal < Pdotb - 0.0001 || marginal > Pdotb + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(b): "<<Pdotb<<" != P(Ab) "<<prob.pAb()<<" + P(ab) "<<prob.pab()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
//End Sanity Check

				float dprime=0.0;
				float delta=((float)prob.pAB()*prob.pab())-((float)prob.pAb()*prob.paB());

				if (delta>0) {
					float m = fmin((float)Pdotb*PdotA, (float)Pdota*PdotB);
					if (m != 0.0)
						dprime=delta/m;
					else
						dprime=0.0;
				} else {
					float m = fmin((float)PdotA*PdotB, (float)Pdota*Pdotb);
					if (m!=0.0)
						dprime=(delta/m)*-1.0;
					else 
						dprime=0.0;
				}

//More sanity checking
				if (dprime > 1.01) {
					cout<<"\t"<<A.GetLabel()
						<<"\t"<<B.GetLabel()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<dprime
						//<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<lod
						//<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<rsquared
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<PdotA
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<Pdota
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<PdotB
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<Pdotb
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pAB()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pAb()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.paB()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pab()
						<<"\t"<<prob.cAB
						<<"\t"<<prob.cAb
						<<"\t"<<prob.caB
						<<"\t"<<prob.cab			
						<<"\n";
				}
//End sanity checking

				float lod=0.0;
				if (prob.cAB > 0)
					lod+=(prob.cAB*log((float)prob.pAB()/((float)PdotA*PdotB)));
				if (prob.cAb > 0)
					lod+=(prob.cAb*log((float)prob.pAb()/((float)PdotA*Pdotb)));
				if (prob.caB > 0)
					lod+=(prob.caB*log((float)prob.paB()/((float)Pdota*PdotB)));
				if (prob.cab > 0)
					lod+=(prob.cab*log((float)prob.pab()/((float)Pdota*Pdotb)));

				if (lod > 0.0)
					lod*=log10;


				float r=(prob.pAB()*prob.pab()-prob.pAb()*prob.paB());
				float rsquared = r*r/(PdotA*Pdota*PdotB*Pdotb);

				ldStatistics[idx++] = LDStatistics(&A, &B, dprime, lod, rsquared);
				haplotypeHead.Append(P, Q, loci, dprime);
	
				dpLP.AddPoint(dist, dprime, A.GetLabel().c_str(), B.GetLabel().c_str());
				rsLP.AddPoint(dist, rsquared, A.GetLabel().c_str(), B.GetLabel().c_str());
			}
			else {
				haplotypeHead.TruncateBlock( P, Q);
				if (P==Q-1)
					ldStatistics[idx++] = LDStatistics(true);
					
			}

			
		}

		if (ChromPool::continueRunning) {
			ldStatistics[idx] = LDStatistics();
			ldStart[snpWidth]=idx;
		}

		delete[] probabilities;	
		probabilities = NULL;
		sumLocalProcessing+=localProcessing.elapsed();
//		cout<<setw(6)<<P<<" "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(8)<<progress.elapsed()<<"\n";;
	}
	//cout<<"\nCalculating LD Took: "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(8)<<progress.elapsed()<<" ("<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(8)<<sumLocalProcessing<<")\n";
	
	dpLP.Close();
	rsLP.Close();

	haplotypeHead.ValidateBlocks(loci, NULL);

}

void LdPlotter::CalculateLD() {
	assert(0);
	uint first = 0;
	uint last = loci.size() - 1;
	size_t cols, rows;
	CountRC(first, last, rows, cols);


	uint snpWidth = last-first;

	haplotypeHead.Reset(first, last);
	//Initialize the LD cache
	ldStatistics = new LDStatistics[snpWidth*rows + 1];
	ldStart      = new size_t[snpWidth + 1];
	size_t idx = 0;
	bool withinDist = true;

	float log10 = 1/log(10.0);
	
	cout<<"Attempting to write LD data for SNPs "<<loci[first].GetLabel()<<"\t"<<loci[last].GetLabel()<<"\n";
	boost::timer progress;

	//Walk over each snp in the selected region
	for (uint P=first; P<=last; P++) {
		progress.restart();

		withinDist = true;
		Locus &A=loci[P];


		float PdotA=A.Freq1();
		float Pdota=A.Freq2();

		//Add the Snp to the chart so it's label and MAF will be seen	
//		writer->WriteHeader( P, A );
		ldStart[P] = idx;
		//uint q1 = loci[P+1].GetID();
		//uint qLast = loci[last].GetID();		
		//From "optimization"
		//Probabilities *probabilities = NULL;
		//CountFrequencies(P, P+1, last, probabilities);

		//uint qidx = 0;
		//Since we are looking at haplotypes, let's inch forward from P one step at a time
		for (uint Q=P+1; Q<=last && withinDist; Q++) {
			Locus &B=loci[Q];
			Probabilities prob;
			//Probabilities &prob = probabilities[qidx++];
			

			uint dist = abs((int)B.GetLocation() - (int)A.GetLocation());

			withinDist = dist < maxSnpDistance;

			//Make sure we don't waste time considering SNPs that are too far away from one another
			if (withinDist) {
				CountFrequencies(A.GetID(), B.GetID(), prob);
				//float Pdotb=B.GetMinAlleleFreq();
				//float PdotB=1.0 - Pdotb;
				float PdotB=B.Freq1();
				float Pdotb=B.Freq2();

//Sanity check
				float marginal = prob.pAB() + prob.pAb();
				if (marginal < PdotA - 0.0001 || marginal > PdotA + 0.0001) {
					cout<<"Marginals aren't adding up for "<<A.GetID()<<"x"<<B.GetID()<<" "<<P<<"x"<<Q<<":\n";
					cout<<"P(A): "<<PdotA<<" != P(AB) "<<prob.pAB()<<" + P(Ab) "<<prob.pAb()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
				marginal = prob.paB() + prob.pAB();
				if (marginal < PdotB - 0.0001 || marginal > PdotB + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(B): "<<PdotB<<" != P(aB) "<<prob.paB()<<" + P(AB) "<<prob.pAB()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
				marginal = prob.paB() + prob.pab();
				if (marginal < Pdota - 0.0001 || marginal > Pdota + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(a): "<<Pdota<<" != P(aB) "<<prob.paB()<<" + P(ab) "<<prob.pab()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}	
				marginal = prob.pAb() + prob.pab();
				if (marginal < Pdotb - 0.0001 || marginal > Pdotb + 0.0001) {
					cout<<"Marginals aren't adding up for "<<P<<"x"<<Q<<":\n";
					cout<<"P(b): "<<Pdotb<<" != P(Ab) "<<prob.pAb()<<" + P(ab) "<<prob.pab()<<"\n";
					cout<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
					A.WriteMarkerInfo(cout, 8);
					B.WriteMarkerInfo(cout, 8);
					prob.Report(cout);
					abort();
				}
//End Sanity Check

				float dprime=0.0;
				float delta=((float)prob.pAB()*prob.pab())-((float)prob.pAb()*prob.paB());

				if (delta>0) {
					float m = fmin((float)Pdotb*PdotA, (float)Pdota*PdotB);
					if (m != 0.0)
						dprime=delta/m;
					else
						dprime=0.0;
				} else {
					float m = fmin((float)PdotA*PdotB, (float)Pdota*Pdotb);
					if (m!=0.0)
						dprime=(delta/m)*-1.0;
					else 
						dprime=0.0;
				}

//More sanity checking
				if (dprime > 1.01) {
					cout<<"\t"<<A.GetLabel()
						<<"\t"<<B.GetLabel()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<dprime
						//<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<lod
						//<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<rsquared
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<PdotA
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<Pdota
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<PdotB
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<Pdotb
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pAB()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pAb()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.paB()
						<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(3)<<prob.pab()
						<<"\t"<<prob.cAB
						<<"\t"<<prob.cAb
						<<"\t"<<prob.caB
						<<"\t"<<prob.cab			
						<<"\n";
				}
//End sanity checking

				float lod=0.0;
				if (prob.cAB > 0)
					lod+=(prob.cAB*log((float)prob.pAB()/((float)PdotA*PdotB)));
				if (prob.cAb > 0)
					lod+=(prob.cAb*log((float)prob.pAb()/((float)PdotA*Pdotb)));
				if (prob.caB > 0)
					lod+=(prob.caB*log((float)prob.paB()/((float)Pdota*PdotB)));
				if (prob.cab > 0)
					lod+=(prob.cab*log((float)prob.pab()/((float)Pdota*Pdotb)));

				if (lod > 0.0)
					lod*=log10;


				float r=(prob.pAB()*prob.pab()-prob.pAb()*prob.paB());
				float rsquared = r*r/(PdotA*Pdota*PdotB*Pdotb);

				ldStatistics[idx++] = LDStatistics(&A, &B, dprime, lod, rsquared);
				haplotypeHead.Append(P, Q, loci, dprime);

			}
			
		}
		cout<<setw(6)<<P<<" "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<progress.elapsed()<<"\n";;
		ldStatistics[idx] = LDStatistics();
		ldStart[snpWidth]=idx;

		//delete[] probabilities;	
	}
	//cout<<"\nCalculating LD Took: "<<progress.elapsed()<<"\n"
	
	haplotypeHead.ValidateBlocks(loci, NULL);

}
/**
 * @NOTE it is very important to realize that the cache is not smart at all. If you initialize the LD Cache to a 
 *       subset of the Loci, and later attempt to render the whole chromosome, you will likely get a crash. Definitely, you
 *       will get questionable results. 
 */
string LdPlotter::WriteLdData(const char *filename, uint first, uint last, LdWriter* writer, uint maxSnpDepth) {
	//This will open the "writer" object and write up any header information
	string pngFilename = writer->Open(filename, loci, first, last);

	if (ldStatistics == NULL) 
		CalculateLD();

	Locus *currentA = NULL;
	size_t lastLdValud = ldStart[last];
	size_t P = first;
	size_t Q = first +1;
	for (size_t idx=ldStart[first]; idx < lastLdValud; idx++ ){

		LDStatistics &s = ldStatistics[idx];
//			cout<<"--"<<s.A->GetLabel()<<" x "<<s.B->GetLabel()<<" "<<s.dprime<<" "<<s.lod<<" "<<s.rsquared<<"\n";
		
		if (s.A || s.B) {
			if (s.A != currentA) {
				currentA = s.A;
//					cout<<"-- "<<setw(10)<<currIdx<<" "<<s.A->GetLabel()<<"\n";
				writer->WriteHeader(P++, *(s.A));
				Q = P;
			}
//				cout<<"\t\t"<<s.A->GetLabel()<<" : "<<s.B->GetLabel()<<"\n";
			if (Q++ <= last)
				writer->Write(*(s.B), s.dprime, s.lod, s.rsquared);
		}
		//Basically, if we have a spot where there is no LD
		else if (s.deadSpace)
			writer->WriteHeader(P++, *(s.A));
	}

	return pngFilename;
}

}
}
