//
// C++ Interface: snprecipient
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESESNPRECIPIENT_H
#define ESESNPRECIPIENT_H

#include "snpaligned.h"

namespace Genetics {

/**
This is the interface for any potential objects used to recieve products of SNP manipulation (such as searches). What happens to the snps once "Appended" is totally dependant on the implimentation (some might get delivered via the network, others written to file, others might be stored for use in more advanced searches)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpRecipient{
public:
    SnpRecipient();

    virtual ~SnpRecipient();

	/**
	 * @brief Sticks a new model into the array. 
	 * The array will grow by a predefined amount in the event of adding more snps than had been initially
	 * allowed for
	 * @param newSnp The snp being added
	 * @param foldID The fold where this is being added
	 */
	virtual void Append(SnpAligned *newSnp) = 0;

	
};



inline
SnpRecipient::SnpRecipient()
{
}

inline
SnpRecipient::~SnpRecipient()
{
}

}

#endif
