/* 
 * File:   analyzer.h
 * Author: torstees
 *
 * The analyzer is a simple structure that allows the user to Analyze
 * a set genes for significance using permutation tests. 
 *
 * Created on February 8, 2010, 10:30 AM
 */

#ifndef _ANALYZER_H
#define	_ANALYZER_H
#include <set>
#include <map>
#include "utility/strings.h"
#include "gene.h"
#include <math.h>
#include "utility/random.h"
namespace Paris {

class Analyzer {
public:
	Analyzer();
	Analyzer(const Analyzer& orig);
	virtual ~Analyzer();

	/**
	 * @brief Basic set of counts useful for reporting
	 */
	struct Result {
		uint complexFeatures;								///< Number of features > 1 SNP
		uint simpleFeatures;									///< Number of 1 SNP features
		uint sigComplex;										///< Number of significant features > 1 SNP
		uint sigSimple;										///< Number of significant 1 SNP Features
		uint geneCount;										///< Number of genes
		uint sigGenes;											///< Number of significan genes
		uint groupID;											///< Tag back to the group
		uint sigPerms;											///< Number of significant permutations
		uint totPerms;											///< Total number of significant
		bool isRefinement;									///< Indicates that this is a refinement
		float stddev;											///< Std dev for N runs

		/**
		 * @brief operator required for STL sorting
		 */
		bool operator<(const Result& other) const {
			assert(isRefinement || other.isRefinement || totPerms == other.totPerms);
			return sigPerms < other.sigPerms;
		}

		/**
		 * @brief Calculates pvalue
		 */
		double PValue() const{
			return double(sigPerms) / double(totPerms);
		}

		bool IsBorderline() const  {
			float pv = PValue();
			return Analyzer::refinementThreshold.first != refinementThreshold.second
					  && pv > Analyzer::refinementThreshold.first
					  && pv < Analyzer::refinementThreshold.second;
		}

		/**
		 * @brief Returns the pvalue in string form, appending the < in cases where there were no permutations to outperform the local one
		 */
		std::string GetPValue() const{
			assert(totPerms > 0);
			if (!isRefinement && sigPerms == 0)
				return "< " + Utility::ToString(1.0/float(totPerms), (int)log10(totPerms));
			else
				return Utility::ToString(float(sigPerms)/float(totPerms), (int)log10(totPerms));
		}

		void SetAsRefinement(float mean, float stdev) {
			sigPerms = (uint)mean;
			this->stddev = stdev/(float)totPerms;
			isRefinement = true;
		}

		Result() : complexFeatures(0), simpleFeatures(0), sigComplex(0), sigSimple(0),
										geneCount(0), sigGenes(0), groupID(0), sigPerms(0), totPerms(0), isRefinement(false), stddev(0.0) {}
		Result(uint groupID) : complexFeatures(0), simpleFeatures(0), sigComplex(0), sigSimple(0),
										geneCount(0), sigGenes(0), groupID(groupID), sigPerms(0), totPerms(0), isRefinement(false), stddev(0.0) {}
	};


	/**
	 * @brief Perform analysis.
	 * @param reportID will be written to the result object that is produced
	 * @param genes one or more genes to be permuted
	 * @param bins The bins used for permutation
	 * @param permutationCount How many permutations will be run
	 * @param verbose When true, a summary for the analysis will be printed to std::out
	 */
	Analyzer::Result Analyze(uint reportID, std::set<Gene*>& genes, std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose=true);

	/**
	 * @brief Perform analysis.
	 * @param reportID will be written to the result object that is produced
	 * @param genes one or more genes to be permuted
	 * @param bins The bins used for permutation
	 * @param ncCount The number of negative control permutations to perform. If this is zero, then the regular analysis should be performed
	 * @param permutationCount How many permutations will be run
	 * @param verbose When true, a summary for the analysis will be printed to std::out
	 */
	Analyzer::Result AnalyzeForNegativeControls(uint reportID, std::set<Gene*>& genes, std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose=false);

	Feature *GetNextBin(std::set<uint>& visited, std::vector<Feature*>& features);
	/**
	 * @brief Required for permutations
	 * The seed must be explicitly assigned on OSX, since there seems to be some differences in how static references
	 * are handled.
	 */
	static Utility::Random &rnd;
	static std::pair<float, float> refinementThreshold;		///< inclusive boundaries for means which are considered borderline. Set these to be the same to never return borderline status
private:

};


} //Paris
#endif	/* _ANALYZER_H */

