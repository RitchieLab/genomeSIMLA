//
// C++ Interface: penetranceeval
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPENETRANCEEVAL_H
#define SIMULATIONPENETRANCEEVAL_H

#include <map>
#include <vector>
#include <string>
#include <math.h>

namespace Simulation {

using namespace std;

/**
@brief Handles stores and evaluates a penetrance based on a number of basic statistics

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PenetranceEval{
public:

class PenCell;

struct Allele {
	double freq;
	char id;
	
	Allele() : freq(0.0), id('\0') {}
	Allele(const Allele& other) : freq(other.freq), id(other.id) {} 
	Allele(double freq, char id) : freq(freq), id(id) { }

};

struct Locus {
	Allele al1, al2;
	
	Locus() { }
	Locus(const Locus& other) : al1(other.al1), al2(other.al2) { } 

	Locus(double f1, double f2, char id1, char id2) : al1(f1, id1), al2(f2, id2) { }

	/**
	 * @brief less than for sorting purposes
	 */
	bool operator<(const Locus& other) const { 
		return (al1.id < other.al1.id) || (al2.id < other.al2.id); 
	}

	/**
	 * @brief returns the genotype label for given type (0, 1, 2) for 3 possible genotypes
	 */
	string GetGenotypeLabel(int gt);
};

struct PenCell {
	double penetrance;
	double frequency;
	string id;

	PenCell() : penetrance(0.0), frequency(0.0), id("")  {} 
	PenCell(const char *id, double freq) : penetrance(0.0), frequency(freq), id(id) {}
	PenCell(const PenCell& other) : penetrance(other.penetrance), frequency(other.frequency), id(other.id) { }
	
	/**
	 * @brief Creates a new node that represents the cell associated with the local cell and the new allele
	 */
	PenCell Append(const Allele& al1, const Allele& al2) const;

	/**
	 * @brief returns the prevalence for the cell (case)
	 */
	double GetPrevalence();

	/**
	 * @brief returns the prev. of the cell if, and only if, the genotypeLabel is a member of the cell
	 * @return prev or 0.0 if genotype isn't found in the id
	 */
	double GetPrevalence(string genotypeLabel);
	/**
	 * @brief returns the prevalence for the cell (control)	
	 */
	double GetAntiPrevalence();
	
	/**
	 * @brief returns the control prev. of the cell if, and only if, the genotypeLabel is a member of the cell
	 * @return prev(control) or 0.0 if genotype isn't found in the id
	 */
	double GetAntiPrevalence(string genotypeLabel);

	//void Normalize(map<string, double>& prev, const double& overall);

};


struct NmiCalculator {
	struct NmiCell {
		double cases; 
		double controls;
		string label;

		NmiCell() : cases(0.0), controls(0.0), label("") { }
		NmiCell(double& cases, double& controls, const char *label) : cases(cases), controls(controls), label(label) { }
		
		double GetSum() { return cases+controls;} 
		double GetEntropy() { double r1=cases+controls; return r1*log(r1); }
		double GetT() { double r1=cases/(cases+controls); return r1*log(r1); }
		double GetU() { double r1=controls/(cases+controls); return r1*log(r1);}

		double GetCofAB() { 
			double sum=cases+controls;
			double percCases = cases/sum;
			double percControls = controls/sum;
			return ( (percCases*log(percCases) + percControls*log(percControls)) * sum);
		}
		
	};

	NmiCalculator() : sumCases(0.0), sumControls(0.0) { }
	void Append(double& cases, double& controls, const char *string) { cells.push_back(NmiCell(cases, controls, string)); sumCases+=cases; sumControls+=controls; }	
	double Sum();	
	double NMI();

	void Report(ostream& os);
private:
	double sumCases;
	double sumControls;
	vector<NmiCell> cells;
};




    PenetranceEval() {}
    ~PenetranceEval() {}

	/**
	 * @brief Appends a new locus to the locus array
	 */
	void AddLocus(const double& freq1, const double& freq2, char id1, char id2);

	/**
	 * @brief Clears the loci from the table and removes all penetrance values
	 */
	void Clear();

	/**
	 * @brief Initializes the penetrance cells based on the loci stored locally
	 */
	int InitializeCells();
	
	/**
	 * @brief recursively build the penetrance table
	 */
	void AddPenCell(const PenCell &curCell, vector<Locus>::iterator itr);

	/**
	 * @brief Returns the number of loci associated with the table
	 */
	int GetLocusCount() { return loci.size(); }
	
	/**
	 * @brief returns a reference to the cell at label, label. 
	 */
	PenCell& GetCell(const char *label);

	void Normalize(map<string, double>& prev, const double& overall);

	/**
	 * @brief performs analysis on the table using various statistical tests
	 * @return score associated with the evaluation
	 */
	int Evaluate(ostream& os);

	/**
	 * @brief Evaluates risk of 1 or more genotypes from a prevalence table
	 */
	double EvaluateRisk(map<string, double> &table, vector<string>& genotypes);

	void EvaluateSingleLocus(ostream& os, const double& prev, const double& anti);
	void EvaluateNMI(ostream& os);
	void DumpPrevalences(ostream& os);

protected:
	vector<Locus> loci;
	map<string, PenCell> penCells;
	double EvaluateMarginals(ostream &os);
};






}

#endif
