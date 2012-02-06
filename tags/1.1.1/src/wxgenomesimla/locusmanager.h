//
// C++ Interface: locusmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUILOCUSMANAGER_H
#define GENOMESIM_GUILOCUSMANAGER_H

#include <string>
#include <vector>
#include <map>
#include "simulation/locus.h"
#include "simulation/chromosome.h"

namespace GenomeSIM {

namespace GUI {

using namespace std;
using namespace Simulation;


struct FileBasedChromosome {

typedef enum FileStatus {
	Unknown, 
	Closed,
	Open,
	Error
} FileStatus;

	FileBasedChromosome(uint chrId) : filename(""), chrID(chrId), status(Unknown){}
	FileBasedChromosome(uint chrId, const char *label, const char *filename) 
							: label(label), filename(filename), chrID(chrId), status(Unknown) { 
		cout<<"FileBasedChromosome("<<chrId<<", "<<label<<", "<<filename<<")\n";
	}
	FileBasedChromosome(uint chrID, const char *label, const char *filename, 
							LocusArray *loci) : label(label), 
							filename(filename), chrID(chrID), status(Unknown) {
		cout<<"FileBasedChromosome("<<chrID<<", "<<label<<", "<<filename<<",...)\n";
		size_t count = loci->size();
		for (size_t i=0; i<count; i++) {
			Locus *l = new Locus((*loci)[i]);
			this->loci.push_back(l);
		}
	 }
	~FileBasedChromosome();
	/**
	 * @brief Load the contents of the file into memory
	 */
	bool Load();

	/**
	 * @brief Save the contents to file 
	 */
	void Save();


	bool Load(const char *filename);

	/**
	 * @brief Saves the locus report under the specified
	 */
	void Save(const char *filename);
	
	/**
	 * @brief Remove all loci from the vector
	 */
	void Purge();
	
	
	string label;				///<The short named associated with the chromosome
	string filename;			///<The filename
	uint chrID;

	size_t GetLocusCount() { return loci.size(); }
	Locus *GetLocus(const char *lbl);
	Locus *GetLocus(int idx) { if (loci.size() > 0) return loci[idx]; else return NULL; }
	int GetLocusIdx(const char *lbl) { return locusLookup[lbl]; }
	vector<Locus *> &GetLoci() { return loci; }
	vector<Locus> ExtractLoci();
	map<string, int> &GetLocusLookup() { return locusLookup; }

	bool HasValidFile();
	
	/**
	 * @brief Allow the chromosome to request it's filename.
	 * @return t/f indicating a valid file was selected
	 * @TODO This should be moved into a subclass that is wx-specific
	 */
	bool RequestFile();


protected:
	vector<Locus *> loci;				///<The loci associated with it
	map<string, int> locusLookup;	///<For looking up loci by the label
	Locus *locBegin;
	Locus *locEnd;
	FileStatus status;
};




}

}

#endif
