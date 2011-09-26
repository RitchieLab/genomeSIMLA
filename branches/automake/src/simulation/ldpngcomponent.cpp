//
// C++ Implementation: lddprimepngwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldpngcomponent.h"
#include "utility/strings.h"
namespace Simulation {

namespace Visualization {

uint ImageParameters::standardBlockSize 	= 32;
uint ImageParameters::tinyBlockSize     	= 2;
uint ImageParameters::medBlockSize	   		= 8;
string ImageParameters::font 		   		= "FreeMonoBold.ttf";
uint ImageParameters::maxSnpsPerRow	   		= 3500;
uint ImageParameters::imageCompressionLevel	= 3;

ImageParameters::ImageParameters() : margin(xy(65, 30)), snpDepth(0), snpCount(0), dimensions(xy(0, 0)), 
						spacerSize(2), fontSize(10), halfFontSize(5), headerFontSize(14), tinyFontSize(6),
						blockSize(standardBlockSize), halfWidth((uint)(standardBlockSize*0.5)+spacerSize), png(NULL),
						showMAF(true), showLabels(true), showPositions(true), showLdValues(true), linkPositions(true), blockIdThickness(3) {}

void ImageParameters::Init(int snpCount, int snpDepth) {
	//Here we work out the necessary details for each of the various sizes
	spacerSize 			= 2;
	
	//This represents the total room for all header components except chart label
	headerHeight 		= 300;

	blockSize 			= standardBlockSize;

	labelFontSize		= 30;

	fontSize			= 12;
	//Gotta scale things down if we get huge things to display
	if (snpCount > 500) {
		blockSize		= tinyBlockSize;
		headerHeight 	= 75;
		spacerSize 		= 0;
		fontSize		= 1;
		halfFontSize	= 1;
		showLdValues = showPositions = showMAF = showLabels = false;
		linkPositions = false;
		blockIdThickness = 1;
		labelFontSize 	= 12;
		margin.x = 10;
	}
	else if (snpCount > 100)  {
		blockSize 		= medBlockSize;
		headerHeight 	= 250;
		spacerSize 		= 1;
		fontSize		= 6;
		halfFontSize	= 3;
		blockIdThickness = 2;
		showLdValues 	= false;
		labelFontSize	= 20;
		margin.x = 35;
	}

	halfWidth = (int)(0.5 * (float)blockSize + spacerSize);

	snpWidth = snpCount;
	this->snpDepth = snpDepth;
	this->snpCount = snpCount;
	wrapCount = 1;

	if (snpWidth > maxSnpsPerRow) {
		snpWidth = maxSnpsPerRow;
		wrapCount = snpCount / snpWidth + 1;

		if (snpCount % maxSnpsPerRow < (float)maxSnpsPerRow * 0.5) {
			snpWidth = snpCount / (wrapCount - 1);
			wrapCount--;
		}
		else {
			snpWidth = snpCount / wrapCount;
		}
	}
	dimensions.x = (int)(margin.x * 2 + ((snpWidth + 1) * halfWidth * 2));
	maxChartWidth = snpWidth * halfWidth * 2;
	maxChartHeight = snpDepth * halfWidth + headerHeight + halfWidth;
	//Each sub-chart is snpDepth * halfwidth + some header and spacing details
	dimensions.y = (int)(margin.y * 2 + maxChartHeight * wrapCount);	

}

/**
 * @brief Returns the xy coordinates for the top left corner of a given block
 * @note  This will take wrap around in consideration
 */
xy ImageParameters::GetBlockPosition(size_t idx1, size_t idx2) {
	//First, determine which "segment" we start with
	xy pos(margin.x + (idx1 + 1) * halfWidth * 2 + halfWidth * (idx2-idx1), dimensions.y - margin.y - headerHeight - (idx2 - idx1) * halfWidth);

	//If this is to be wrapped....
	if (wrapCount > 1) {
		size_t column = idx1 + (size_t)((idx2 - idx1) *0.5);
		size_t seg = column / snpWidth;

		if (seg > wrapCount)
			seg=wrapCount;

		pos.x-= maxChartWidth * seg;
		pos.y-=maxChartHeight*seg;
	}
	
	return pos;	
}

void ImageParameters::PlotText(int fontSize, const char *text, int x, int y, float rotation, int r, int g, int b) {
	if (Utility::FileExists((char *)font.c_str()))
		png->plot_text((char *)font.c_str(), fontSize, x, y , rotation, (char *)text, r, g, b);
}	

void ImageParameters::PlotText(const char *text, int x, int y, float rotation, int r, int g, int b) {
	PlotText(fontSize, text, x, y, rotation, r, g, b);

}

/**
 * @brief Returns the xy coordinates for the top center of a snps header information
 * @note This will take wrap around in consideration
 */
xy ImageParameters::GetHeaderPosition(size_t idx) {
	//First, find the theoretical X coordinate 
	xy pos(margin.x + (idx + 1) * halfWidth * 2, dimensions.y - margin.y - headerHeight);

	if (wrapCount > 1) {
		size_t seg = idx / snpWidth;

		if (seg >= wrapCount)
			seg=wrapCount-1;
		pos.x-= maxChartWidth * seg;
		pos.y-=maxChartHeight*seg;
	}
	return pos;	
}

string ImageParameters::Open(const char *filename) {
	char pngFilename[2048];
	sprintf(pngFilename, "%s.png", filename);

	png = new pngwriter(dimensions.x, dimensions.y, 0.65, pngFilename);
	png->setcompressionlevel(imageCompressionLevel);
	return pngFilename;
}
	
ImageParameters::~ImageParameters() {
	if (png) {
		png->close();
		delete png;
	}
}



}

}
