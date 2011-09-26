//
// C++ Implementation: genegenemodelarchive.h
//
// Description: Serializes gene gene models from binary file
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>
#include <fstream>
#include "genegenemodel.h"
#include "utility/filebuffer.h"
#ifndef __BIO_GENE_GENE_MODEL_ARCHIVE_H
#define	__BIO_GENE_GENE_MODEL_ARCHIVE_H

namespace Biofilter {

/**
	@Brief Manages the gene models
	@Note This only considers 2 genes at a time. If we want 3 gene models, we need to replace this with a  more dynamic structure
	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class GeneGeneModelArchive {
public:
	GeneGeneModelArchive(const char *tempFilename = "models.gg", uint modelInit=1000000, uint modelMax=500000000, bool binaryStorage=true);
	~GeneGeneModelArchive();

	std::map<float, uint> Archive(const char *geneArchive, const char *modelArchive);
	void DisplayContents(std::ostream& os);
	void Insert(GeneGeneModel& model);

	uint SummarizeModelCounts(std::map<uint, uint>& scores, std::map<uint, Region*>& regions);


	class iterator {
		friend class GeneGeneModelArchive;
	public:
		iterator(const iterator& other);
		bool operator==(const iterator& other) const;
		iterator& operator++();
		GeneGeneModel& operator*();
	protected:
		iterator(std::ifstream& source, bool isEnd = false);
		GeneGeneModel curModel;
		bool isEnd;
		std::ifstream& source;
	};
	iterator begin();
	iterator end();



	void WriteGeneArchive(std::ostream& os, bool useBinary);
	void WriteGeneArchive(const char *filename);
	//void WriteGeneArchiveBinary(const char *filename);
protected:

	Utility::FileBuffer<GeneGeneModel> models;
	std::ifstream source;
	bool binaryStorage;											///< Binary or ascii storage
};

inline
GeneGeneModelArchive::iterator::iterator(const GeneGeneModelArchive::iterator& other)
	: isEnd(other.isEnd), source(other.source) { }

inline
GeneGeneModelArchive::iterator::iterator(std::ifstream& source, bool isEnd) : isEnd(isEnd), source(source) { }

inline
bool GeneGeneModelArchive::iterator::operator==(const GeneGeneModelArchive::iterator& other) const {
	if (isEnd || other.isEnd)
		return isEnd == other.isEnd;
	return source == other.source;
}

inline
GeneGeneModelArchive::iterator& GeneGeneModelArchive::iterator::operator++() {
	if (source.good() && !source.eof()) {
		curModel.LoadBinary(source);
	}
	return *this;
}


inline
GeneGeneModel& GeneGeneModelArchive::iterator::operator*() {
	return curModel;
}
}

#endif //__BIO_GENE_GENE_MODEL_ARCHIVE_H
