//
// C++ Implementation: lineparsing
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineparsing.h"
#include "utility/lineparser.h"
namespace Utility {

namespace Test {

LineParsing::LineParsing()
{
}


LineParsing::~LineParsing()
{
}

void LineParsing::TestExclusionList() {
	LineParser lp;
	BitSetType values(10, 0ul);
	values=~values;				//Let's set it back to all ones 
	
	ExclusionList exclusions(&values, false);
	uint lines=lp.Parse("ExcludeSNPs.lst", &exclusions);
	CPPUNIT_ASSERT(lines == 3);
	CPPUNIT_ASSERT(values.count() == 7);
	CPPUNIT_ASSERT(!values[3]);			//4
	CPPUNIT_ASSERT(!values[4]);			//5
	CPPUNIT_ASSERT(!values[9]);			//10
	
}


}

}
