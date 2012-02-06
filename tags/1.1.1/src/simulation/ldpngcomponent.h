//
// C++ Interface: lddprimepngwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDDPRIMEPNGWRITER_H
#define SIMULATION_VISUALIZATIONLDDPRIMEPNGWRITER_H
#include "utility/types.h"
#include "pngwriter/pngwriter.h"
#include <string>

namespace Simulation {

namespace Visualization {

using namespace std;

struct xy {
	size_t x;
	size_t y;
	xy(size_t x, size_t y) : x(x), y(y) {}
};

//Most of the properties will be configurable, but there are a few
//that will be static for general settings per configuration
struct ImageParameters {
	xy margin;							///<Number of pixels are associated with the margins

	int snpDepth;						///<This is the real depth for a given string of ld blocks
	int snpCount;						///<Total number of snps. 

	xy dimensions;						///<x- width  y- height

	int spacerSize;						///<number of pixels between 
	int fontSize;						///<The size for normal text
	int halfFontSize;					///<
	int headerFontSize;					///<The size for header information
	int tinyFontSize;					///<For labeling graph axes

	int blockSize;						///<Actual size of the diamonds
	int halfWidth;						///<1/2 the number of pixels between each block

	size_t headerHeight;				///<This is the total height allocated per all header components that appear with each segment	
	size_t snpWidth;					///<Actual number of snps per chart segment
	size_t wrapCount;					///<How many chart segments there are
	size_t maxChartWidth;				///<Just the max width of any segment (width - 2*margin)
	size_t maxChartHeight;				///<Height of a given chart segment

	pngwriter *png;

	bool showMAF;						///<Allow MAF to be turned on/off according to size of graph
	bool showLabels;					///<Allow labels to be turned on/off according to size of graph
	bool showPositions;					///<Allow positions to be turned on/off according to size of graph
	bool showLdValues;					///<Write LD values onto the grid? (only on highest resolution graphs)
	bool linkPositions;					///<Do we want to render the line between the position and the snps?
	int  labelFontSize;

	static uint standardBlockSize;		///<Configurable size for normal blocks
	static uint tinyBlockSize;			///<Configurable size for tiny blocks
	static uint medBlockSize;			///<Configurable size for med. sized blocks

	static uint maxSnpsPerRow;			///<Maximum number of snps allowed on a single row
	static uint imageCompressionLevel;	///<0 no compression. 9- highest compression
	
	size_t blockIdThickness;			///<How wide a line is for designating a block
	
	string filename;					///<The name of the file associated with the image object

	ImageParameters();

	void Init(int snpCount, int snpDepth);

	/**
	 * @brief Returns the xy coordinates for the top left corner of a given block
	 * @note  This will take wrap around in consideration
	 */
	xy GetBlockPosition(size_t idx1, size_t idx2);

	/**
	 * @brief Returns the xy coordinates for the top center of a snps header information
	 * @note This will take wrap around in consideration
	 */
	xy GetHeaderPosition(size_t idx);

	string Open(const char* filename);


	void PlotText(const char *text, int x, int y, float rotation = 0.0, int r=0, int g=0, int b=0);
	void PlotText(int fontSize, const char *text, int x, int y, float rotation=0.0, int r=0, int g=0, int b=0);
	
	~ImageParameters();


	static string font;					///<The font used to write render text to the screen

protected:

};

class LdPngComponent {
public:
	LdPngComponent(ImageParameters *param) : imgParams(param) {}
	virtual ~LdPngComponent() {}
	virtual void Open(const char *text) = 0;
	virtual void Close() = 0;
	/**
	 * @brief Give the component it's starting offset, and it will return the lower offset
	 * @note This is strictly for stacked components
	 */
	virtual int Init(int yOffsetStart, uint first, uint last) = 0;

protected:
	ImageParameters *imgParams;
};



}

}

#endif
