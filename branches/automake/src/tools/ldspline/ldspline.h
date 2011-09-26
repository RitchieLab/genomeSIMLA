/* 
 * File:   ldspline.h
 * Author: torstees
 *
 * Created on August 30, 2010, 12:16 PM
 */

#ifndef LDSPLINE_H
#define	LDSPLINE_H

#include <map>
#include <vector>
#include <fstream>
#include "utility/strings.h"
#include "locuslookup.h"



class LdSplineApp {
public:
	LdSplineApp();
	LdSplineApp(const LdSplineApp& orig);
	virtual ~LdSplineApp();


	void LoadHeadersFromHapmap(const char *chrom, const char* filename);
	void LoadLdFromHapmap(const char *filename);
	void LoadFromBinary(const char *filename);
	void SaveToBinary(const char *filename);
	void SaveToCopyBinary(const char *newFilename);
	void RunReport(std::ostream& os);
	void RunReport(const char *chrom, int position, float value, const char *ldType, std::ostream& os);
	void ReportBounds(const char *chrom, int position, float value, const char *type, std::ostream& os);
	std::string ConvertHaploLD(const char *filename);
	void ExportForLiftOver(const char *bimOrig);
	void ImportLiftOver(const char *loFilename, const char *loUnmapped);
	void OpenBinary(const char *filename, bool loadFullHeaders = false);
	void Summarize(std::ostream& os, const char *chrom);
	std::vector<SnpSpline> GetLocusRange(const char *chrom, int start, int stop);
private:
	std::map<std::string, LocusLookup> loci;										///< This is how we record our loci. The spline will return only indexes
	std::map<std::string, std::string> filenames;								///< Used to record the hapmap files so that we can parse them for LD after we've set up the headerspace
	std::fstream file;
};

/**
class SnpSpline {
public:
	SnpSpline(int idx, int rsid, int pos) :pos(pos), rsid(rsid), idx(idx) {}
	~SnpSpline() {}

	void AddStatsUpstream(int idx, float dp, float rs) {
		int lIdx = idx - this->idx;
		if (lIdx >= upstream.size())
			upstream.resize(lIdx);
		upstream[lIdx] = LdStat(dp, rs);
	}

	void AddStatsDownstream(int idx, float dp, float rs) {
		int lIdx = this->idx - idx;
		if (lIdx >= downstream.size())
			downstream.resize(lIdx);
		downstream[lIdx] = LdStat(dp, rs);
	}

	std::vector<int> GetSplineDP(float minDP) {
		std::vector<int> indexes;
		std::vector<Locus>::iterator itr = downstream.begin();
		std::vector<Locus>::iterator end = downstream.end();

		int i=0;
		while (itr!=end && (itr->dp == -1 || itr->dp>=minDP)) {
			i++;
			if (itr->dp > 0)
				indexes.push_back(this->idx - i);
			itr++;
		}

		std::reverse(indexes.begin(), indexes.end());

		itr = upstream.begin();
		end = upstream.end();
		i = 0;
		while (itr!=end && (itr->dp == -1 || itr->dp>=minDP)) {
			i++;
			if (itr->dp > 0)
				indexes.push_back(this->idx + i);
			itr++;
		}

		return indexes;
	}

	int pos;
	int rs;
	int idx;
protected:
	std::vector<LdStat> upstream;
	std::vector<LdStat> downstream;
};

class LocusLookup {
	LocusLookup() { }
	~LocusLookup() { }

	void AddLocus(int rs, int pos) {
		
		std::map<int, int>::iterator itr = posToIdx.find(pos);
		
		if (itr != loci.end()) {
			int idx = loci.size();
			loci.push_back(SnpSpline(idx, pos, rs));
			posToIdx[pos] = idx;
		}
	}

	int GetIndex(int pos) {
		if (posToIdx.find(pos) != posToIdx.end())
			return posToIdx[pos];
		return -1;
	}
	int GetRS(int idx) {
		return loci[idx].rs;
	}
	int GetPos(int idx) {
		return loci[idx].pos;
	}

	void AddLdValue(int pos1, int pos2, int dp, int rs) {
		int id1 = GetIndex(pos1);
		int id2 = GetIndex(pos2);

		//For now, we are under the assumption that pos1<pos2
		loci[id1].AddStatsUpstream(id2, dp, rs);
		loci[id2].AddStatsDownstream(id1, dp, rs);
	}

protected:
	std::vector<SnpSpline> loci;
	std::map<int, int> posToIdx;
};


*/



#endif	/* LDSPLINE_H */

