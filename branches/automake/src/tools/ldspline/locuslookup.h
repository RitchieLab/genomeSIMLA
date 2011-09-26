#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include "snpspline.h"
#include <set>



class LocusLookup {
public:
	typedef std::map<int, float> SplineLookup;

	typedef std::map<int, int> RsToPos;
	typedef std::map<int, int> PosToRS;

	LocusLookup(std::fstream* file = NULL, const char *chromosome = "");
	~LocusLookup();

	void AddLocus(int rs, int pos, uint64_t offset);

	int GetIndex(int pos);
	int GetRS(int idx);
	int GetPos(int idx);

	int GetSnpCount();
	void IncrementSplineCounts(int pos1, int pos2);
	void AddLdValue(int pos1, int pos2, float dp, float rs);
	int64_t DumpBinaryHeader(int64_t offset);
	int64_t DumpBinaryHeader(std::fstream* file, int64_t offset);
	void DumpBinary();
	void DumpBinary(std::fstream* file);

	std::pair<int, int> GetLdSplineBoundsRS(int pos, float value);
	std::pair<int, int> GetLdSplineBoundsDP(int pos, float value);

	void LoadBinaryHeader();
	void LoadBinary();

	void LoadHeadersFromHapmap(const char *filename);
	void LoadLdFromHapmap(const char *filename);
	void WriteMapForLiftOver(std::ostream& os);
	void LoadMapFromLiftOver(std::istream& infile, std::set<int>& droppedPositions);
	void Release();
	void SkimBinaryHeader();
	void LoadLocusDetails();
	std::string Chromosome();

	void Summarize(std::ostream& os);
	std::vector<SnpSpline> GetLocusRange(int start, int stop);
	/**
	 * Returns the spline range in positions for dprime, dp
    * @param pos - the position of the item you are
    * @return
    */
	SplineLookup GetLdSplineDP(int pos, float dp);
	SplineLookup GetLdSplineRS(int pos, float rs);

	RsToPos GetRsToPos();			///< Returns rs->pos
	int GetPosToRS(int Pos);		///< Returns the RS number or -1
	PosToRS GetPosToRS();			///< Returns pos->rs

	int LocusCount();

protected:
	std::string chromosome;
	std::vector<SnpSpline> loci;
	std::map<int, int> posToIdx;
	std::fstream* file;
	int locusCount;
	std::fstream::off_type headerLoadPosition;
};



