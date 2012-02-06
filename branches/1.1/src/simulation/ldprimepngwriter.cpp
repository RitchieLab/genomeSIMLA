//
// C++ Implementation: ldprimepngwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldprimepngwriter.h"
#include <sstream>
#include <iomanip>

namespace Simulation {

namespace Visualization {

using namespace std;

LdPrimePngWriter::LdPrimePngWriter(ImageParameters *params, uint start, uint stop) : 
	imgParams(params), labels(params), locations(params, start, stop, true), maf(params), 
	footer(params), bRed(0), bGreen(0), bBlue(0)
{
}


LdPrimePngWriter::~LdPrimePngWriter()
{
	delete imgParams;
}


void LdPrimePngWriter::AddBlock(BlockListNode *block, uint classification) {
//AddBlock(HaplotypeBlock *block, uint classification) {

	int firstSnp = block->GetFirstIndex();
	int firstIdx = firstSnp - firstSnpIndex;
	int lastSnp = 1 + block->GetLastIndex();
	int lastIdx = lastSnp - firstSnpIndex;

	if (firstIdx < 0 || lastSnp > (int)(lastSnpIndex+1)) {
		return;
	}

//	cout<<"Block "<<block->GetLabel()<<": "<<firstSnp<<"-"<<lastSnp<<" into "<<firstSnpIndex<<"-"<<lastSnpIndex<<"\n";
	
	int blockThickness = classification * imgParams->blockIdThickness;

	int nextSnp = firstIdx;
	size_t blockWidth = imgParams->halfWidth;
	size_t rightEdge = imgParams->GetHeaderPosition(0).x;


	if (firstSnp >= (int)lastSnpIndex)
		return;

	xy nextXY = imgParams->GetHeaderPosition( firstIdx );
	xy firstXY = nextXY;

	labels.AddBlock( block, classification );

	uint xOffset = 0;
	uint yOffset = 0;

	//Draw left side 
	while (nextSnp < lastIdx && nextSnp < (int)(lastSnpIndex - firstSnpIndex)) {
		firstXY = nextXY;

		size_t x1 = firstXY.x-blockWidth-xOffset;
		size_t y1 = firstXY.y - yOffset;
		size_t x2 = x1;
		size_t y2 = y1;

 		do {
			x2+=blockWidth;
			y2-=blockWidth;
			yOffset+=blockWidth;
			xOffset+=blockWidth;
			nextXY = imgParams->GetHeaderPosition( ++nextSnp);

		} while (nextSnp < lastIdx && nextSnp <= (int)(lastSnpIndex - firstSnpIndex) && x2 <= imgParams->margin.x + imgParams->maxChartWidth); 

		//Draw the line	
		for (int i=-1; i<blockThickness; i++) {
			imgParams->png->line(x1, y1-i, x2, y2-i, bRed, bGreen, bBlue);
		}
	}	
	
	//Draw right side 
	nextSnp = lastIdx;
	nextXY = imgParams->GetHeaderPosition( nextSnp-1);
	nextXY.x+=blockWidth*2;
	yOffset = 0;
	xOffset = 0;
	while (nextSnp < (int)lastSnpIndex && nextSnp > firstIdx && nextSnp > 0) {
		firstXY = nextXY;

		size_t x1 = firstXY.x-blockWidth+xOffset;
		size_t y1 = firstXY.y + yOffset;
		size_t x2 = x1;
		size_t y2 = y1;
		//cout<<"Searching for the end of the block on this 'segment'\n";
		//Find the various pieces
		do {
			x2-=blockWidth;
			y2-=blockWidth;
			if (nextSnp > lastSnp) {
			//if (x1 > imgParams->dimensions.x - imgParams->margin.x) {
				x1 = x2;
				y1 = y2;
			}
				
			xOffset+=blockWidth+imgParams->halfWidth;
			yOffset-=blockWidth-imgParams->halfWidth;
			nextXY = imgParams->GetHeaderPosition( --nextSnp );
		} while (nextSnp > firstIdx && x2 > rightEdge && nextSnp > 0);		

		if (x1 > x2) {
			//Draw the line	
			for (int i=-1; i<blockThickness; i++) 
				imgParams->png->line(x1, y1-i, x2, y2-i, bRed, bGreen, bBlue);
		}
	}	
	
}

string LdPrimePngWriter::Open(const char *filename, vector<Locus> &loci, uint first, uint last) {
	string pngFilename = imgParams->Open(filename);
	Init(loci, first, last);
	labels.Open( filename );
	locations.Open( filename );
	maf.Open(filename);
	return pngFilename;
}

void LdPrimePngWriter::Init(vector<Locus> &loc, uint first, uint last) {
	vector<Locus>::iterator cur = loc.begin();
	vector<Locus>::iterator end = loc.end();

	int initialY = imgParams->headerHeight;			// - imgParams->margin.y;
	locations.SetTracerToLabel(imgParams->linkPositions);
	initialY = locations.Init(initialY, first, last);
	maf.Init(initialY, first, last);
	labels.Init(initialY, first, last);
	footer.Init(initialY, first, last);
	uint i = 0;

	for (; cur != end && i<= last; cur++) {
		if (i++>=first) {
			Locus &l = *cur;
			locations.AddLocus(l);
			labels.AddLocus(l);
			maf.AddLocus(l);
		}
	}

	idx1 = idx2 = (size_t)-1;

	firstSnpIndex = first;
	lastSnpIndex = last;
}
	
void LdPrimePngWriter::Close() {
	labels.Close();
	locations.Close();
	maf.Close();
	footer.Close();
	imgParams->png->close();
}

void LdPrimePngWriter::SetBlockDensity(uint min, uint max) {
	labels.SetBlockDensityBounds(min, max);
}

void LdPrimePngWriter::WriteHeader( uint idx, Locus &A) {
	idx1++;
	idx2 = idx1+1;
}

void LdPrimePngWriter::Write(Locus &B, float &dprime, float &lod, float &rsquared) {
	xy pos = imgParams->GetBlockPosition( idx1, idx2++);
	float r=1.0, g=.090, b=0.90;
	//Shades of red/pink
	if (lod > 2) {
		if (dprime < 0.5) {
			g=0.878431373;
			b=0.956862745;
		}
		else {	
			g=b=((float)((255-32)*2)*(1.0-dprime)*0.003921569);
		}
	}
	//Light blue
	else if (dprime > 0.99) {
		r=g=0.75;
		b=0.878431373;

	}
	//Whatever drive down the b/g as dprime gets higher.
	else {
		r=g=b=0.95;
		g=b=((float)((255-32)*2)*(1.0-dprime)*0.003921569);
	}
	stringstream dpValue;
	dpValue<<(int)(100*dprime);
	imgParams->png->filleddiamond(pos.x, pos.y, imgParams->blockSize, imgParams->blockSize, r,g,b);
	
	if (imgParams->showLdValues && dprime < 1.0) {
		int textSize = (int)(imgParams->png->get_text_width((char *)imgParams->font.c_str(), imgParams->halfFontSize, (char *)dpValue.str().c_str()) * 0.5);
		int textHeight = imgParams->tinyFontSize;

//		cout<<"plot_text("<<(char *)imgParams->font.c_str()<<", "<<imgParams->halfFontSize<<", 0, "<<
//				(uint)pos.x - textSize<<", "<<pos.y - textHeight<<", 0, "<<(char *)dpValue.str().c_str()<<", "<<
//				bRed<<", "<<bGreen<<", "<<bBlue<<"\n";
/*		imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->halfFontSize, 
					(uint)pos.x - textSize, pos.y - textHeight, 
					0, (char *)dpValue.str().c_str(), 
					bRed, bGreen, bBlue);	
*/
		imgParams->PlotText(imgParams->halfFontSize, dpValue.str().c_str(), (int)pos.x-textSize, pos.y-textHeight);
	}


}



LdRSquaredPngWriter::LdRSquaredPngWriter(ImageParameters *params, uint start, uint stop) : 
	LdPrimePngWriter(params, start, stop)
{
	bRed 	= 1.0;
	bGreen 	= 0.2;
	bBlue 	= 0.2;
}

LdRSquaredPngWriter::~LdRSquaredPngWriter() { }

void LdRSquaredPngWriter::Write(Locus &B, float &dprime, float &lod, float &rsquared) {
	xy pos = imgParams->GetBlockPosition( idx1, idx2++);
	float r, g, b;
	r=g=b = 1.0 - rsquared;
	imgParams->png->filleddiamond(pos.x, pos.y, imgParams->blockSize, imgParams->blockSize, r,g,b);

	stringstream dpValue;
	dpValue<<(int)(100*rsquared);
	imgParams->png->filleddiamond(pos.x, pos.y, imgParams->blockSize, imgParams->blockSize, r,g,b);


	if (imgParams->showLdValues && rsquared < 1.0) {
		int textSize = (int)(imgParams->png->get_text_width((char *)imgParams->font.c_str(), imgParams->halfFontSize, (char *)dpValue.str().c_str()) * 0.5);
		int textHeight = imgParams->tinyFontSize;

/*//Just to figure out which are uninitialized
if (imgParams->halfFontSize != 0)
	cerr<<imgParams->halfFontSize<<"\n";
if (dpValue.str().length() > 0)
	cerr<<"THat works. That's good\n";
if (pos.x > 0)
	cerr<<"Pos x\n";
if (textSize > 0)
	cerr<<"Text Size\n";
if (pos.y > 0)
	cerr<<"pos y\n";
if (textHeight > 0)
	cerr<<"Text Height\n";
*/
		imgParams->PlotText(imgParams->halfFontSize, dpValue.str().c_str(), (int)pos.x-textSize, pos.y-textHeight);
/*		imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->halfFontSize, 
					(uint)pos.x - textSize, pos.y - textHeight, 
					0, (char *)dpValue.str().c_str(), 
					bRed+rsquared, bGreen+rsquared, bBlue+rsquared);	
*/
	}
}


}

}
