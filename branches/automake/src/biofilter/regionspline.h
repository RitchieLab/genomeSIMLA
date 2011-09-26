//
// C++ Interface: regionspline
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERREGIONSPLINE_H
#define BIOFILTERREGIONSPLINE_H

#include <iomanip>
#include <iostream>
#include <soci.h>
#include <soci-sqlite3.h>

#include <vector>
#include <map>

#include "snpmanager.h"

namespace Biofilter {

using namespace soci;

struct SplinePoint {
	SplinePoint(uint point) : upper(point), lower(point), boundaryUpper((uint)-1), boundaryLower(0) { }
	SplinePoint() : upper(0), lower(0), boundaryUpper(-1), boundaryLower(0) { }
	SplinePoint(const SplinePoint& other)
		: upper(other.upper), lower(other.lower), boundaryUpper(other.boundaryUpper), boundaryLower(other.boundaryLower) { }
	uint upper;				///< Highest known good location (this can be increased as long as it is lower than boundary Upper)
	uint lower;				///< Lowest known good location (this can be decreased as long as it is higher than the lower boundary)
	uint boundaryUpper;		///< The lowest (greater than the current) point that fails to meet the minimum LD value
	uint boundaryLower;		///< The highest point (lower than the current) that fails to meet the minimum LD value
};
class RegionBoundary {
public:
	RegionBoundary(uint start, uint stop, uint popID, float value) : start(start), end(stop), popID(popID), ldValue(value) { }


	~RegionBoundary() {}
	
	/**
	 * @Brief Very basic evaluation. This assumes the SNP is on chromosome and that the other SNP does lie within the real boundaries
	 * @Return False if the LD Value doesn't meet the requirements (doesn't indicate that we actually modified anything)
	 */
	bool Evaluate(uint cur, uint other, float ld) {
		if (points.find(cur) == points.end()) {
			points[cur] = SplinePoint(cur);
		}
		SplinePoint& point = points[cur];
		bool success = false;
		//Expanding
		if (ld >= ldValue) {
			//We are Expanding to the right
			if (cur < other) {
				if (other < point.boundaryUpper && other > point.upper) {
					point.upper = other;
					success = true;
				}
			}
			//Expanding to the left
			else if (cur > other) {
				if (other > point.boundaryLower && other < point.lower) {
					success = true;
					point.lower = other;
				}
			}
		}
		//Contraction
		else {
			//Upper boundary
			if (cur < other) {
				if (other < point.boundaryUpper) {
					point.boundaryUpper = other;
					if (other <= point.upper)
						point.upper = cur;
					success = true;
				}
			}
			else {
				if (other > point.boundaryLower) {
					point.boundaryLower = other;
					if (other > point.lower)
						point.lower = cur;
					success = true;
				}
			}
		}
		return true;
	}

	void Commit(std::ostream& os, uint geneID, uint origStart, uint origEnd) {
		if (points.size() == 0)
			return;

		std::map<uint, SplinePoint>::iterator itr = points.begin();
		std::map<uint, SplinePoint>::iterator end = points.end();

		while (itr != end) {
			if (itr->second.lower < start) {
				start = itr->second.lower;
				if (origStart - itr->second.lower > 200000 )
					std::cout<<" { "<<origStart<<" : "<<itr->second.lower<<" : "<<itr->second.boundaryLower<<" } ";
			}
			if (itr->second.upper > this->end) {
				this->end = itr->second.upper;
				if (this->end - origEnd > 200000)
					std::cout<<" { "<<origEnd<<" : "<<itr->second.upper<<" : "<<itr->second.boundaryUpper<<" } ";
			}
			itr++;
		}

		assert(this->end >= origEnd);
		assert(this->start <= origStart);
		assert(this->start < this->end);

		if (this->end - origEnd > 200000)
			std::cout<<"**";
		if (origStart - this->start > 200000)
			std::cout<<"**";


		if (origStart > start || origEnd < this->end) {
			std::cout<<"<(";
			if (origStart > start)
				std::cout<<std::setw(10)<<(origStart - start);
			else
				std::cout<<std::setw(10)<<0;
			std::cout<<") (";
			if (this->end > origEnd)
				std::cout<<std::setw(10)<<(this->end - origEnd);
			else
				std::cout<<std::setw(10)<<0;
			std::cout<<")>\t";
			os<<"UPDATE region_bounds SET start="<<start<<", end="<<this->end<<" WHERE gene_id="<<geneID<<" AND population_id="<<popID<<";\n";
			//sociDB << "UPDATE region_bounds SET start=:start, end=:end WHERE gene_id=:geneID AND population_id=:popID", use((int)start), use((int)this->end), use((int)geneID), use((int)popID);
		}
		else { std::cout<<"<("<<std::setw(10)<<0<<") ("<<std::setw(10)<<0<<")>\t"; }

	}
	/**
	 * @brief Only performs the commit if start and stop differe from the local versions
	 */
	uint Commit(soci::session& sociDB, uint geneID, uint origStart, uint origEnd) {
		if (points.size() == 0)
			return 0;
		std::map<uint, SplinePoint>::iterator itr = points.begin();
		std::map<uint, SplinePoint>::iterator end = points.end();

		while (itr != end) {
			if (itr->second.lower < start) {
				start = itr->second.lower;
//				if (origStart - itr->second.lower > 200000 )
//					std::cout<<" { "<<origStart<<" : "<<itr->second.lower<<" : "<<itr->second.boundaryLower<<" } ";
			}
			if (itr->second.upper > this->end) {
				this->end = itr->second.upper;
//				if (this->end - origEnd > 200000)
//					std::cout<<" { "<<origEnd<<" : "<<itr->second.upper<<" : "<<itr->second.boundaryUpper<<" } ";
			}
			itr++;
		}

		assert(this->end >= origEnd);
		assert(this->start <= origStart);
		assert(this->start < this->end);

		if (origStart > start || origEnd < this->end) {
			if (abs(origStart - start) > 600000 | abs(origEnd - this->end) > 600000) {
				std::cerr<< "Oversized LD region identified:\n";
				std::cerr<< " - "<<start<<" - "<<origStart<<" = "<<abs(origStart-start)<<"\n";
				std::cerr<< " - "<<this->end<<" - "<<origEnd<<" = "<<abs(origEnd-this->end)<<"\n";
			}
			//assert(abs(origStart - start) < 600000);
			//assert(abs(origEnd - this->end) < 600000);
/*			std::cout<<"<(";
			if (origStart > start)
				std::cout<<std::setw(10)<<(origStart - start);
			else
				std::cout<<std::setw(10)<<0;
			std::cout<<") (";
			if (this->end > origEnd)
				std::cout<<std::setw(10)<<(this->end - origEnd);
			else
				std::cout<<std::setw(10)<<0;
			std::cout<<")>\t";
*/
			sociDB << "UPDATE region_bounds SET start=:start, end=:end WHERE gene_id=:geneID AND population_id=:popID", use((int)start), use((int)this->end), use((int)geneID), use((int)popID);
		}
//		else { std::cout<<"<("<<std::setw(10)<<0<<") ("<<std::setw(10)<<0<<")>\t"; }
		return abs(origStart - start) + abs(origEnd-this->end);
	}

	uint PointCount() { return points.size(); }
protected:
	std::map<uint, SplinePoint> points;	///< points associated with known boundaries
	uint start;
	uint end;
	uint popID;
	float ldValue;
};	

/**
@Brief Store all regions boundaries and their bounds for each different LD makeup

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class RegionSpline{
public:
    RegionSpline(int geneID, int chrom, uint start, uint stop);
    ~RegionSpline();

	/**
	 * @brief Add pairwise LD and potentially expand the boundaries
	 */
	bool AddSnps(SNP_Details* first, SNP_Details* last, float dprime, float rsquared);
	bool AddSnps(uint first, uint last, int chromosome, float dprime, float rsquared);
	void Print(std::ostream& os) {
		os<<"["<<start<<","<<end<<"] ";
	}

	/**
	 * @Brief Write the regions which were updated by LD
	 */
	void Commit(soci::session& sociDB);
	void Commit(std::ostream& os);
	
	static void AddRS(float ldValue, uint popID);
	static void AddDP(float ldValue, uint popID);

	static std::map<float, uint> dprime;			///< dprime->popID
	static std::map<float, uint> rsquared;		///< rsquared->popID
	uint start;								///< real gene start
	uint end;								///< real gene end
protected:
	void InitBoundaries();

	int chrom;								///< Used to eliminate invalid SNPs 
	int geneID;								///< Used for tracking
	std::vector<RegionBoundary> dprimeBounds;	///< Boundaries for each dprime threshold
	std::vector<RegionBoundary> rsquaredBounds;	///< Boundaries for each rsquared threshold
};




}

#endif
