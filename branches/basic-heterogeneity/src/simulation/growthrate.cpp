//
// C++ Implementation: growthmodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, 
//				(C) Marylyn Ritchie 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <cmath>
#include "growthrate.h"
#include "exponentialgrowth.h"
#include "logisticgrowth.h"
#include "lineargrowth.h"
#include "richardslogistic.h"
#include "ldpngcomponent.h"
#include "chrompool.h"
#include "utility/strings.h"


namespace Simulation {
namespace PopulationGrowth {


using namespace Visualization; 

uint  GrowthRate::maxPoolSize 	   = 100000;	///<App wide max size. This is valid for ALL models
uint  GrowthRate::minPoolSize 	   = 500;		///<App wide min
//float GrowthRate::variation   	   = 0.0005;	///<App wide variation to be applied to size calculation
uint  GrowthRate::randomSeed 	   = 1371;
uint  GrowthRate::adjustmentPeriod = 2048;



void GrowthRate::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Fixed Population Size: " << maxPoolSize<<endl;
}

/*** No Growth *************************************************/
GrowthRate::GrowthRate(float variation) : rnd(randomSeed), modelType(GrowthRateUnknown) { 
	draws = new float[adjustmentPeriod];

	SetVariation(variation);
}

GrowthRate::~GrowthRate() { 
	if (draws)
		delete[] draws;
}

uint GrowthRate::GetPopulationSize(uint t) { 
	return maxPoolSize;
}

uint GrowthRate::operator()(uint t) {
	return maxPoolSize;
}

float GrowthRate::GetVariation() {
	return variation;
}

void GrowthRate::SetVariation(float var) { 
	variation = var;

	double min = 1.0 - (variation * 0.5);


	for (uint i=0; i<adjustmentPeriod; i++) {
		draws[i] = rnd(variation) + min;
	}
}

double GrowthRate::AdjustedGrowth(double nsize) {

	if (!isnormal(nsize))
		return maxPoolSize;
	//double adjustment = ((float)GetInitialPopulationSize() * rnd(variation)) - ((float)GetInitialPopulationSize()* variation/2.0);
	return ((float)nsize * draws[(uint)nsize%adjustmentPeriod]);
}

uint GrowthRate::GetInitialPopulationSize() {
	return maxPoolSize;
}

string GrowthRate::GetType() {
	return "Fixed";
}

string GrowthRate::GetHtmlReport() {
	stringstream ss;
	ss<<"<CENTER>Fixed population: "<<maxPoolSize<<"</CENTER>\n";
	return ss.str();
}

GrowthRate *GrowthRate::LoadModel(istream& ss) {
	GrowthRate *model = NULL;
	string type;
	ss>>type;
	
	if (strcmp(type.c_str(), "EXPONENTIAL") == 0) {
		float rate, var;
		uint initialPop;
		ss>>initialPop>>var>>rate;
		if (rate > 0.0) 
			model = new ExponentialGrowth(initialPop, rate, var);
		else {
			cout<<"Unable to create an exponential growth model with a growth rate of 0.0. Please use linear if you want no growth\n";
			model = NULL;		
		}

	}
	else if (strcmp(type.c_str(), "LOGISTIC") == 0) {
		uint carry;
		float rate, var; 
		uint initialPop;
		ss>>initialPop>>var>>rate>>carry;
	
		if (rate > 0.0) {
			model = new LogisticGrowth(initialPop, rate, carry, var);
//			((LogisticGrowth *)model)->FindInflection(0, 2000);
		}
		else {
			cout<<"Unable to create an exponential growth model with a growth rate of 0.0. Please use linear if you want no growth\n";
			model = NULL;
		}
	}
	else if (strcmp(type.c_str(), "RICHARDS") == 0) {
		float rate, polarity, var;
		uint carry, initialPop, timeOfMaxGrowth;
		uint resetGen=0;
		ss>>initialPop>>var>>rate>>carry>>timeOfMaxGrowth>>polarity>>resetGen;
		
		if (rate > 0.0) {
			model = new RichardsLogistic(initialPop, carry, timeOfMaxGrowth, rate, polarity, var, resetGen);
//			((RichardsLogistic *)model)->FindInflection(0, 2000);
		}
		else		{
			cout<<"Unable to create an exponential growth model with a growth rate of 0.0. Please use linear if you want no growth\n";
			model = NULL;
		}
	}
	else if (strcmp(type.c_str(), "LINEAR") == 0) {
		float rate, var;
		uint initialPop;
		ss>>initialPop>>var>>rate;

		model = new LinearGrowth(initialPop, rate, var);
	}
	return model;
}
string GrowthRate::DrawGrowthChart(const char *project, int start, int stop, int interval, uint width, uint height, bool returnPage) {
	if (interval < 1)
		interval = 1;

	uint fontSize1 = 16, fontSize2 = 12, fontSize3 = 10;
	bool showLabels = true;
	if (width < 600 || height < 350) {
		showLabels = false;
		fontSize1 = 8; fontSize2 = 8; fontSize3=6; }

	uint popBreadth = (uint)((GetPopulationSize(stop) - GetPopulationSize(start)) * 1.05);
	if (popBreadth == 0)
		popBreadth = (uint)(GetPopulationSize(start) * 1.5);
	uint genBreadth = stop - start;

	string pngFilename = string(project) + ".growth.png";
	string htmlFilename = string(project) + ".growth.html";
	xy dimension(width, height);									//Image dimensions

	float margin = 0.1, length = 0.8;

	if (!showLabels) {
		margin = 0.05; 
		length = 0.9;
	}
	xy axis((uint)((float)dimension.x * margin), (uint)((float)dimension.y * margin));		//The x/y positions of the axes
	xy lengths((uint)((float)dimension.x * length),(uint)( (float)dimension.y * length));		//The lengths of the lines
	float xScale = (float)(dimension.x * length) / genBreadth;
	float yScale = (float)(dimension.y * length) / popBreadth;

	string font = Simulation::Visualization::ImageParameters::font;
	
	pngwriter png(dimension.x, dimension.y, 0.95, pngFilename.c_str());
	ofstream report(htmlFilename.c_str(), ios::out);
	report<<"<HTML><HEADER><TITLE>Population Growth Rate</TITLE></HEADER><BODY>\n";
	report<<"<LINK REL='stylesheet' TYPE='text/css' HREF='"<<ChromPool::cssFilename<<"'>\n";
//	report<<"<DIV CLASS='spacer'></DIV><DIV CLASS='image-frame'>\n";
	report<<"<CENTER><H2>Growth Rate Report</H2>\n";
	report<<GetHtmlReport()<<"\n";
	report<<"<DIV CLASS='spacer'></DIV>\n";

	report<<"<IMG SRC='"<<Utility::ExtractFilename(pngFilename.c_str())<<"'></IMG>\n<DIV CLASS='spacer'></DIV>\n";
	report<<"<TABLE><TR><TH>Generation</TH><TH>PopuationSize</TH></TR>\n";
	//Let's draw the axes
	png.line(axis.x, axis.y, dimension.x - axis.x, axis.y, 0, 0, 128);
	png.line(axis.x, axis.y, axis.x, dimension.y - axis.y , 0, 0, 128);
	
	char label[1024];
	uint startY = GetPopulationSize(start);
	xy lastPosition(axis.x + (int)(start * xScale), axis.y);
	bool invert = false;
	bool fontExists = false;

	if (showLabels) {
		sprintf(label, "%s", Utility::ExtractFilename(project).c_str());
		if (Utility::FileExists(font.c_str())) {
			fontExists = true;
			png.plot_text((char *)font.c_str(), fontSize1, (uint)(dimension.x * 0.5 - (0.5 * png.get_text_width((char *)font.c_str(), fontSize1, label))), dimension.y - 20, 0.0, label, 0,0,0);
		
			sprintf(label, "Growthrate - %s", GetType().c_str());
			png.plot_text((char *)font.c_str(), fontSize2, (uint)(dimension.x * 0.5 - (0.5 * png.get_text_width((char *)font.c_str(), fontSize2, label))), dimension.y - 35, 0.0, label, 0,0,0);
		
			sprintf(label, "Generation");
			png.plot_text((char *)font.c_str(), fontSize3, (uint)(dimension.x * 0.5 - (0.5 * png.get_text_width((char *)font.c_str(), fontSize3, label))), 5, 0.0, label, 0,0,0);
		
			sprintf(label, "Population (K)");
			png.plot_text((char *)font.c_str(), fontSize3, 15, (uint)(dimension.y * 0.5 - (0.5 * png.get_text_width((char *)font.c_str(), fontSize3, label))), 1.570795, label, 0,0,0);
		}
		else
			cout<<"Unable to open file, '"<<font<<"'. Charts will not have labels\n";
	
	}
	for (int i=start; i<=stop; i+=interval)  {
		//Draw the line from the old position to the new
		xy newPosition ((int)(axis.x + ((i - start) * xScale)), (int)(axis.y + (GetPopulationSize(i) - startY) * yScale));

		if (invert)
			report<<"<TR CLASS='invert'>";
		else
			report<<"<TR               >";
	
		report<<"<TD>"<<i<<"</TD><TD>"<<GetPopulationSize(i)<<"</TD></TR>\n";
		sprintf(label, "G(%u)=%u", (uint)newPosition.x, (uint)newPosition.y);

		png.line(axis.x+1, newPosition.y, newPosition.x, newPosition.y, 0.75, 0.75, 0.75);
		png.line(newPosition.x, axis.y+1, newPosition.x, newPosition.y, 0.75, 0.75, 0.75);

//		png.line(lastPosition.x, lastPosition.y, newPosition.x, newPosition.y, 0, 0, 0);
//		png.filledsquare(lastPosition.x - 2, lastPosition.y +2, lastPosition.x + 2, lastPosition.y-2, 0, 0, 0);
		lastPosition = newPosition;		
	}		
	report<<"\n</TABLE>\n";
	png.filledsquare(lastPosition.x - 2, lastPosition.y +2, lastPosition.x + 2, lastPosition.y-2, 0, 0, 0);


	lastPosition = xy(axis.x + (int)(start * xScale), axis.y);
	for (int i=start; i<=stop; i+=interval)  {
		//Draw the line from the old position to the new
		xy newPosition (int(axis.x + ((i - start) * xScale)), (int)(axis.y + (GetPopulationSize(i) - startY) * yScale));

		sprintf(label, "G(%u)=%u", (uint)newPosition.x, (uint)newPosition.y);

		png.line(lastPosition.x, lastPosition.y, newPosition.x, newPosition.y, 0, 0, 0);
		png.filledsquare(lastPosition.x - 2, lastPosition.y +2, lastPosition.x + 2, lastPosition.y-2, 0, 0, 0);
		lastPosition = newPosition;		
	}		

	//Let's draw the axes
	png.line(dimension.x - axis.x, dimension.y - axis.y, dimension.x - axis.x, axis.y, 0, 0, 128);
	png.line(dimension.x - axis.x, dimension.y - axis.y, axis.x, dimension.y - axis.y, 0, 0, 128);


	float gridPoints = 0.1;
	if (!showLabels)
		gridPoints = 0.25;

	if (showLabels) {
		//Let's draw labels and the grid!
		for (float i=0.0; i<=1.0; i+=gridPoints) {
			uint curY = (uint)(axis.y + (i * popBreadth * yScale));
			png.line(axis.x, curY, dimension.x - axis.x, curY, .35, .35, 0.35);
			sprintf(label, "%d", (int)(i*popBreadth)+GetPopulationSize(start));
	
			if (showLabels && fontExists) {
				png.plot_text((char *)font.c_str(), 7, axis.x - png.get_text_width((char*)font.c_str(), 7, label) - 5, curY, 0.0, label, 0, 0, 0);
			}
	
			uint curX = (uint)(axis.x + (i*genBreadth * xScale));
			png.line(curX, axis.y, curX, dimension.y - axis.y, .35, .35, 0.35);
			
			if (showLabels && fontExists) {
				sprintf(label, "%d", (int)((double)(i*genBreadth))+start);
				png.plot_text((char *)font.c_str(), 7, curX, axis.y - png.get_text_width((char*)font.c_str(), 7, label) - 5, 1.570795, label, 0, 0, 0);
			}
		}
	}

	

	png.close();
	if (returnPage)
		return htmlFilename;
	else 
		return pngFilename;
}


void GrowthRate::DiagramGrowth(std::ostream& os, int start, int stop, int interval) {
	os<<"Generation,Population Size\n";
	for (int i=start; i<=stop; i+=interval) 
		os<<i<<","<<GetPopulationSize(i)<<endl;
	//rnd.Seed( GrowthRate::randomSeed );
}


}

}
