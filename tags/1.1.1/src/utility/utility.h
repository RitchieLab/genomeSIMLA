//
// C++ Interface: Utility
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
/** \page utility Utility Library
 * \section utility_intro Introduction
 * This Library offers a number of basic classes usable for many of the classes and libraries
 * expected to be found in the various programs.
 * - High Level Classes \ref util_design
 * 
 */

/**
 * \page util_design High Level Classes
 * Most of the classes within should be relatively simple. However, below is a 
 * list of the more interesting groups:
 * <P><UL>
 *  	<LI>\ref Utility::StringLookup - This class offers the ability to transform
 * 		textual fields into distrete integer fields. It also offers the ability to 
 * 		dump the mapping to a stream for reporting purposes. </LI>
 *    </UL>
 */
#ifndef UTILITY_H
#define UTILITY_H
/**
 * Basic library containing various utilities used throughout the applications and libraries
 */
#include "utility/types.h"
#include "utility/stringlookup.h"
#include "utility/basiclog.h"
#include "utility/genolookup.h"
#include "utility/genobplookup.h"
#include "utility/generatereport.h"
#include "utility/genotypeparser.h"
#include "utility/binarrayparser.h"
#include "utility/lineparser.h"
#include "utility/configurationparser.h"
#include "utility/application.h"
#include "utility/casecontrolstatus.h"



#endif //UTILITY_H
