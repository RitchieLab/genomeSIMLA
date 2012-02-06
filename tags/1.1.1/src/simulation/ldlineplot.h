//
// C++ Interface: ldlineplot
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LDLINEPLOT_H
#define LDLINEPLOT_H

#include <map>
#include <string>
#include <fstream>
#include "ldpngcomponent.h"

namespace Simulation {
namespace Visualization {

using namespace std;


/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdLinePlot {
public:

struct BinType;

struct BinType {
	float total;
	uint count;
	vector<string> labels;
	BinType() : total(0), count(0) { }

	BinType(const BinType& other) : total(other.total), count(other.count) { }

	void Append(float value, const char *lbl1, const char *lbl2) { 
		char newLabel[2048];
		sprintf(newLabel, "%sx%s(%.4f)", lbl1, lbl2, value);
		total+=value; 
		count++;
		labels.push_back(newLabel);
	}
/*	BinType &operator=(const BinType& other) { 
		total=other.total; 
		count=other.count; 
		return *this;
	}*/
	float operator+=(float value) { total+=value; count++; return total; }
	float GetMean() { return total / (float)count; }
	void Report(ostream& os) {
		os<<GetMean()<<"\t[";
		for (uint i=0; i<labels.size(); i++) 
			os<<labels[i]<<"\t";
		os<<"\t]\n";
	}

};

	/**
	 * @brief Builds a line plot for visualizing R2 and D'
	 * @param parameters This contains the details relating to the image itself
	 * @param maxDistance the maximum distance between two snps to be considered
	 */
    LdLinePlot(ImageParameters *parameters, uint maxX);
    ~LdLinePlot();

	void Close();
	string Open(const char *filename);

	void AddPoint(uint y, float value, const char *lbl1, const char *lbl2);


	void SetLabel(char *lbl);
	void SetXAxisLabel(char *lbl);
	void SetYAxisLabel(char *lbl);
protected:
	void Prep();
	void PlotPoint(uint x, uint y);
	uint maxX;
	std::map<uint, BinType> points;
	string label;
	string xAxisLabel;
	string yAxisLabel;
	ImageParameters *imgParams;
	uint binWidth;
	uint granularity;	
	
	uint minY;
	uint maxY;

	float yWidth;
	ofstream txtFile;
};

inline
LdLinePlot::LdLinePlot(ImageParameters *imgParams, uint maxX) : maxX(maxX), imgParams(imgParams) {
	if (imgParams) {
		granularity = imgParams->dimensions.x - 2 * imgParams->margin.x;
		binWidth = maxX / granularity;
	}
	else {
		granularity = maxX;
		binWidth 	= 1000;
	}
}

inline
void LdLinePlot::AddPoint(uint x, float value, const char *lbl1, const char *lbl2) {
	if (x <= maxX)  {
		uint bin = x / binWidth;
		points[bin].Append(value, lbl1, lbl2);
	}
}




inline
void LdLinePlot::SetLabel(char *lbl) {
	label = lbl;
}

inline
void LdLinePlot::SetXAxisLabel(char *lbl) {
	xAxisLabel = lbl;
}

inline
void LdLinePlot::SetYAxisLabel(char *lbl) {
	yAxisLabel = lbl;
}

}
}

#endif
