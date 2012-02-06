//
// C++ Interface: templatedpedigree
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PEDIGREETEMPLATESTEMPLATEDPEDIGREE_H
#define PEDIGREETEMPLATESTEMPLATEDPEDIGREE_H

#include <vector>
#include <map>
#include <queue>
#include "pedigreesample.h"
namespace Simulation {
namespace PedigreeTemplates {

struct IndID {
	uint pedID;
	uint indID;

	IndID(const IndID& other) : pedID(other.pedID), indID(other.indID) {}
	IndID(uint pID, uint iID) : pedID(pID), indID(iID) { }
	~IndID() { }

	void Print() const {
		cout<<pedID<<"\t"<<indID<<"\n";
	}
	
	bool operator<(const IndID& other) const {
		if (pedID == other.pedID) 
			return indID < other.indID;
		else
			return pedID < other.pedID;
	}
};

class TemplatedIndividual {
public:
	uint pedID;				///< Pedigree ID
	uint id;				///< Individual ID
	uint patID;				///< Paternal ID
	uint matID;				///< Maternal ID
	int status;				///< Status   Unknown =0, Unaffected=1, Affected=2, Ambiguous=3
	int gender;				///< Gender
	bool validStatus;		///< Person of Interest....was the individual's status evaluated by a clinician?
	bool hasGenotypes;		///< Indicates whether or not this person was genotyped
	
	TemplatedIndividual() : pedID(0), id(0), patID(0), matID(0), status(0), gender(0), validStatus(false), hasGenotypes(false) { 
	}
	~TemplatedIndividual() { 
	}

	/**
	 * @brief Load from pedigree data, return true if there is a valid pedigree ID (> 0)
	 */
	bool Load(const char *line) {
		if (line[0] == '#' || strlen(line) < 1)
			return false;
		string gt1, gt2;
		stringstream ss(line);
		ss>>pedID>>id>>patID>>matID>>gender>>status>>gt1>>gt2;
		validStatus		= status > 0;
		hasGenotypes 	= (gt1 != "0" && gt2 != "0");
		return pedID > 0;
	}
	bool IsAffected() { return status == 2; }

	bool IsAmbiguous() { return status == 0; }	

	bool IsFounder() { return matID == 0 && patID == 0; }

	IndID GetPatID() { return IndID(pedID, patID); }
	IndID GetMatID() { return IndID(pedID, matID); }
	IndID GetID() { return IndID(pedID, id); }
};


/**
@Classes associated with a templated pedigree system

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class TemplatedPedigree {
public:
	struct StatusCounts {
		vector<int> status;				///< Unk, Unaff, Aff, Amb
		int totals;						///< Total number of entries
		int cycles;						///< Cycles, used for getting average for missed reporting
		int attempts;					///< Used to pass the total number of attempts to the calling routine
		void Append(int status) {
			assert(status >= 0 && status < 4);
			this->status[status]++;
			totals++;
		}
		StatusCounts() : status(4,0), totals(0), cycles(1), attempts(1) { }
		StatusCounts(int a, int u, int amb, int unk) : status(4, 0){
			status[0]=amb;
			status[1]=u;
			status[2]=a;
			status[3]=unk;
			totals=(unk+u+a+amb);
		}
		int CountPOI() {
			return totals-status[3];
		}
		int AffectedCount() {
			return status[2];
		}
		int AmbiguousCount() {
			return status[0];
		}
		float AmbiguousRatio() {
			return (float)status[0]/(float)(totals-status[3]);
		}
		void Reset() {
			status = vector<int>(4, 0);
			cycles = 1;
			totals = 0;
		}

	
		StatusCounts& operator+(const StatusCounts& other) {
			if (totals>0) {
				attempts++;
				cycles++;
			}
			status[0]+=other.status[0];
			status[1]+=other.status[1];
			status[2]+=other.status[2];
			status[3]+=other.status[3];
			totals+=other.totals;
			return *this;
		}
		void Header(ostream& os) {
			os<<setw(40)<<" "
				<<setw(15)<<"Unknown"
				<<setw(15)<<"Unaffected:"
				<<setw(15)<<"Affected:"
				<<setw(15)<<"Ambiguous"
				<<setw(15)<<"Total"
				<<setw(15)<<"Attempts"<<"\n";
		}
		void Report(ostream& os, const char *comment) {
			os<<setw(40)<<comment
				<<setw(15)<<status[3]/cycles
				<<setw(15)<<status[1]/cycles
				<<setw(15)<<status[2]/cycles
				<<setw(15)<<status[0]/cycles
				<<setw(15)<<totals/cycles
				<<setw(15)<<attempts;
			if (cycles > 1)
				os<<" ("<<cycles<<") tries";
			os<<"\n";
		}

		bool Evaluate(StatusCounts& other, float tolerance) {
			// For now, we'll just look at affected count
			int variation = abs(status[2] - other.status[2]) + abs(status[1] - other.status[1]);
			return (float)variation/(float)(status[2]+status[1]) <= tolerance;
		}

	};
	typedef std::map<IndID, TemplatedIndividual*> IndividualLookup;
    TemplatedPedigree();
    ~TemplatedPedigree();

	/**
	 * @Brief Loads template from file, returning number of individuals present
	 */
	int Load(const char *filename);

	/**
	 * @brief Grabs an individual if both parents are 0 or are present in members, otherwise returns NULL
	 */
	Individual *DrawIndividual(PoolManager &pools, TemplatedIndividual& tpl);

	/**
	 * @brief Intialize the founders
	 */
	void InitFounders(PoolManager& pools);

	/**
	 * @brief Initialize a new pedigree (using founders)
	 */	
	StatusCounts BuildSample(PoolManager& pools, vector<Individual*>& individuals, PenetranceModel* model, bool failOnError = true);

	/**
	 * @brief return the individual associated with the ped,ind IDs. 
	 * @return If none exists, NULL is returned
	 */
	TemplatedIndividual *GetMember(uint pedID, uint indID);
	
	/**
	 * @Brief Return the number of individuals in the reference set
	 */
	size_t GetMemberCount();

	void PurgeMembers();



	static float tolerance;							///< Tolerance betwee
	static int maxAttempts;							///< # tries used to match the template status layout
	uint ResetAttemptCount();						///< Resets the totalAttempts and returns the number of attempts made so far
protected:
	void ApplyPresentGenotypes(PoolManager &pools, Individual &person);
	/**
	 * @brief Empties the vector, deleting any memory present
	 */
	void PurgeIndividuals(std::vector<Individual*>& individuals);

	IndividualLookup individuals;					///< The templated individuals
	std::vector<Individual*> founders;				///< The founders (individuals without parents)
	std::queue<IndID> origOrder;					///< Used to order the individuals in the pedigree
	StatusCounts statusCounts;						///< A/U counts for persons of interest
	map<IndID, Individual*> members;				///< The current dataset, indexed by IndIDs
	uint totalAttempts;								///< # of attempts associated with building a set of replications
};





}
}

#endif
