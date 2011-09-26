//
// C++ Implementation: regionspline
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "regionspline.h"

namespace Biofilter {

using namespace std;

map<float, uint> RegionSpline::dprime;							///< dprime->popID
map<float, uint> RegionSpline::rsquared;						///< rsquared->popID


RegionSpline::RegionSpline(int geneID, int chrom, uint start, uint stop) : start(start), end(stop),chrom(chrom),geneID(geneID) {
	InitBoundaries();
}

RegionSpline::~RegionSpline() { }


void RegionSpline::AddRS(float ldValue, uint popID) {
	rsquared[ldValue] = popID;
}

void RegionSpline::AddDP(float ldValue, uint popID) {
	dprime[ldValue] = popID;
}

bool RegionSpline::AddSnps(uint first, uint last, int chromosome, float dprime, float rsquared) {
	bool withinBounds = false;
	//First, we have to make sure they are both A) on the same chromosome as we are
	if (chromosome == chrom && first < end && last > start) {
		vector<RegionBoundary>::iterator itr = dprimeBounds.begin();
		vector<RegionBoundary>::iterator end = dprimeBounds.end();
		
		bool success = withinBounds =  true;				///< Cheap short circuit
		if ((first < start && abs((int)first - (int)start) > 600000) || (last > this->end && abs((int)last - (int)this->end) > 600000)){
			cerr<<"OK, we've encountered an attempt to extend a region by over 600,000 bases\n";
			cerr<<"- Real Boundaries: "<<start<<" "<<this->end<<"\n";
			cerr<<"- Extended Bounds: "<<first<<" "<<last<<"\n";
		}
		//DPrime
		while (success && itr != end) {
			//First SNP falls within the real bounds
			if (first >= start && first <= this->end)
				success = itr->Evaluate(first, last, dprime);
			else if (last >= start && last <= this->end)
				success = itr->Evaluate(last, first, dprime);
			itr++;
		}

		success = true;
		itr = rsquaredBounds.begin();
		end = rsquaredBounds.end();
		//RSquared
		while (success && itr != end) {
			//First SNP falls within the real bounds
			if (first >= start && first <= this->end)
				success = itr->Evaluate(first, last, rsquared);
			else if (last >= start && last <= this->end)
				success = itr->Evaluate(last, first, rsquared);
			itr++;
		}
		
	}
	return withinBounds;
}

bool RegionSpline::AddSnps(SNP_Details* first, SNP_Details* last, float dprime, float rsquared) {
	bool withinBounds = false;
	//First, we have to make sure they are both A) on the same chromosome as we are
	if (first->chromosome == chrom && last->chromosome == chrom && first->position < end && last->position > start) {
		bool success = withinBounds =  true;				///< Cheap short circuit

		if (abs((int)first->position - (int)start) > 600000 || abs((int)last->position - (int)this->end) > 600000) {
			cerr<<"OK, we've encountered an attempt to extend a region by over 600,000 bases\n";
			cerr<<"- Real Boundaries: "<<start<<" "<<this->end<<"\n";
			cerr<<"- Extended Bounds: "<<first->position<<" "<<last->position<<"  --- rs"<<first->rsID<<" - rs"<<last->rsID<<"\n";
		}
		//if (first->position < start || last->position>end) {
		{
			vector<RegionBoundary>::iterator itr = dprimeBounds.begin();
			vector<RegionBoundary>::iterator end = dprimeBounds.end();

			//DPrime
			while (success && itr != end) {
				//First SNP falls within the real bounds
				if (first->position >= start && first->position <= this->end)
					success = itr->Evaluate(first->position, last->position, dprime);
				else if (last->position >= start && last->position <= this->end)
					success = itr->Evaluate(last->position, first->position, dprime);
				itr++;
			}

			success = true;
			itr = rsquaredBounds.begin();
			end = rsquaredBounds.end();
			//RSquared
			while (success && itr != end) {
				//First SNP falls within the real bounds
				if (first->position >= start && first->position <= this->end)
					success = itr->Evaluate(first->position, last->position, rsquared);
				else if (last->position >= start && last->position <= this->end)
					success = itr->Evaluate(last->position, first->position, rsquared);
				itr++;
			}
		}
	}
	return withinBounds;
}

void RegionSpline::Commit(std::ostream& os) {
	vector<RegionBoundary>::reverse_iterator itr = dprimeBounds.rbegin();
	vector<RegionBoundary>::reverse_iterator end = dprimeBounds.rend();

cout<<"-> "<<geneID<<" ( "<<itr->PointCount()<<" ) ["<<start<<" "<<this->end<<"] ";
	while (itr != end ) {
		itr++->Commit(os, geneID, start, this->end);
	}

	itr = rsquaredBounds.rbegin();
	end = rsquaredBounds.rend();
cout<<"\t";
	while (itr != end) {
		itr++->Commit(os, geneID, start, this->end);
	}

cout<<"\n";
}


void RegionSpline::Commit(soci::session& sociDB) {
	vector<RegionBoundary>::reverse_iterator itr = dprimeBounds.rbegin();
	vector<RegionBoundary>::reverse_iterator end = dprimeBounds.rend();

//cout<<"-> "<<geneID<<" ( "<<itr->PointCount()<<" ) ["<<start<<" "<<this->end<<"] ";
	while (itr != end ) {
		//itr++->Commit(ss, geneID, start, this->end);
		itr++->Commit(sociDB, geneID, start, this->end);
	}

	itr = rsquaredBounds.rbegin();
	end = rsquaredBounds.rend();
//cout<<"\t";
	while (itr != end) {
//		itr++->Commit(ss, geneID, start, this->end);

		itr++->Commit(sociDB, geneID, start, this->end);
	}

//cout<<"\n";
}

//We could probably do some sort of copy to make this faster, but it would require a bit more intelligence in the boundary class
void RegionSpline::InitBoundaries() {
	map<float, uint>::iterator itr=dprime.begin();
	map<float, uint>::iterator end=dprime.end();

	while (itr!=end) {
		cerr<<"DPrime: "<<itr->second<<"\t"<<itr->first<<"\n";
		dprimeBounds.push_back(RegionBoundary(start, this->end, itr->second, itr->first));
		itr++;
	}

	itr = rsquared.begin();
	end = rsquared.end();
	while (itr!=end) {
		cerr<<"DPrime: "<<itr->second<<"\t"<<itr->first<<"\n";
		rsquaredBounds.push_back(RegionBoundary(start, this->end, itr->second, itr->first));
		itr++;
	}
}
	



}
