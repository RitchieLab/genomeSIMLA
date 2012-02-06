//
// C++ Implementation: snppool
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snppool.h"

namespace Genetics {

uint		SnpPool::snpID = 0;							///<As we create new snps, we give it a unique ID
SnpPool 	*SnpPool::_instance = NULL;					///<Singleton instance
uint 		SnpPool::_instanceCount = 0;				///<The number of uses for cleanup
uint 		SnpPool::cachesize = 250;					///<The size of the original pool
uint 		SnpPool::growby = 25;						///<The number of new items to add when the available pool is empty
uint 		SnpPool::defGenotypeCount = 4;				///<This is the number of genotypes most of the SNPs are expected to have



SnpAligned *SnpPool::GetSnp(uint genotypeCount, uint individualCount /*= 0 */) {

	//This should block anything else from getting in until after the lock goes out of scope
	//boost::mutex::scoped_lock lock(poolMutex);
	//Let's check to see if we can return the last used. This will speed up certain patterns of usage
	SnpAligned *snp= NULL;	
	if (genotypeCount <MAX_EXPECTED_GENOTYPECOUNT && lastSnp[genotypeCount]) {
		snp=lastSnp[genotypeCount];
		snp->Init(genotypeCount, individualCount);
		snp->IncrementInstanceCount();
		lastSnp[genotypeCount]=NULL;	
		return snp;
	}
	//If we can't use that one, let's grab one from the list of available snps
	
	//If the available list is empty, let's get a new chunk
	if (snpsAvailable.size() == 0)
		ResizeAvailable(growby);

	//Let's work through the grap the top snp and mark it as in use
	snp = snpsAvailable.back();
	snp->Init(genotypeCount, individualCount);
	snpsInUse[snp->GetID()] = snp;
	snpsAvailable.pop_back();
	snp->IncrementInstanceCount();
	return snp;

}


void SnpPool::ReleaseSnp(SnpAligned *snp) {
	//This should block anything else from getting in until after the lock goes out of scope
	//boost::mutex::scoped_lock lock(poolMutex);
	if (snp && snp->ReduceInstanceCount() < 1) {
		uint gtCount = snp->GetBitsetCount();
		//If we already have a last snp, we have to prepare it for use
		if (gtCount< MAX_EXPECTED_GENOTYPECOUNT && lastSnp[gtCount] == NULL) 
			lastSnp[gtCount] = snp;

		else {
			SnpMap::iterator itr = snpsInUse.find(snp->id);
			assert(itr != snpsInUse.end());
			snpsAvailable.push_back(itr->second);
			snpsInUse.erase(itr);
		}
	}
}

void SnpPool::ResizeAvailable(uint count) {
	int diff = count - snpsAvailable.size();
	SnpAligned *p = NULL;
	if (diff > 0){
		for (int i=diff; i>0; i--) {
			p = new SnpAligned(defGenotypeCount, snpID++);
			snpsAvailable.push_back(p);
		}
	}
	else	{
		for (int i=diff; i<0; i++)	{
			p = snpsAvailable.back();
			snpsAvailable.pop_back();
			delete p;
			
		}
	}
}

void SnpPool::Reduce() {
	//This should block anything else from getting in until after the lock goes out of scope
//	boost::mutex::scoped_lock lock(poolMutex);
	SnpPool *p = Instance();
	uint size = uint(cachesize * 1.5f);
	if (p->snpsAvailable.size() > size)
		p->ResizeAvailable(size);
	p->Release();
}

void SnpPool::Purge() {
	//This should block anything else from getting in until after the lock goes out of scope
	//boost::mutex::scoped_lock lock(poolMutex);
	ResizeAvailable(0);
	//int count = snpsInUse.size();
	SnpAligned *p = NULL;
	SnpMap::iterator e=snpsInUse.end();
	
	for (uint i =0; i<MAX_EXPECTED_GENOTYPECOUNT; i++) 
		if (lastSnp[i]) { 
			//lastSnp[i]->ReleaseInstance();
			lastSnp[i]=NULL;
		}

	for (SnpMap::iterator i=snpsInUse.begin(); i!=e; i++){
		p=i->second;
		cout<<"Snp in use when pool was destroyed: "<<p->GetID()<<"\n";
		delete p;
	}
	snpsInUse.clear();
	snpID=0;
}

SnpPool::~SnpPool()
{
	Purge();
}


}
