//
// C++ Implementation: ldlineplot
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldlineplot.h"

namespace Simulation {
namespace Visualization {

string LdLinePlot::Open(const char *filename) {
	string pngFilename = string(filename) + string("-") + yAxisLabel + string(".stat");
	if (imgParams)
		string pngFilename = imgParams->Open(filename);
	else 
		txtFile.open(pngFilename.c_str(), ios::out);
	
	//We need to frame the graph and write the labels
	return pngFilename;
}


void LdLinePlot::Prep() {
	std::map<uint, BinType>::iterator itr = points.begin();
	std::map<uint, BinType>::iterator end = points.end();

	minY=1000; 
	maxY = (uint)-1;
	while (itr != end) {
		uint y = (uint)itr->second.GetMean();
		if (y < minY)
			minY = y;
		if (y > maxY)
			maxY = y;
		itr++;
	}
	if (imgParams)
		yWidth = (float)maxY / (imgParams->dimensions.y - imgParams->margin.y);
	else
		yWidth = maxY;

	//Graph and add lines
}

void LdLinePlot::PlotPoint(uint x, uint y) {
	pngwriter *png = imgParams->png;

	png->plot(x + imgParams->margin.x, (uint)(y * yWidth) + imgParams->margin.y, 0, 0, 0);
}

void LdLinePlot::Close() {
	Prep();

	uint count = 0;

	std::map<uint, BinType>::iterator itr = points.begin();
	std::map<uint, BinType>::iterator end = points.end();

	while (itr != end) {
		uint y = (uint)itr->second.GetMean();
		uint x = itr->first;

		if (imgParams)
			PlotPoint(x, y);
		else {
			txtFile<<"1\t"<<x<<"\t";
			itr->second.Report(txtFile);//<<y<<"\t"<<itr->second.count<<"\n";
			count++;
		}
			
		itr++;
	}

	if (imgParams)	
		imgParams->png->close();
}


LdLinePlot::~LdLinePlot()
{
}

}
}


