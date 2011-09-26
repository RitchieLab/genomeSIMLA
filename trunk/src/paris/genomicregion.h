/* 
 * File:   genomicregion.h
 * Author: torstees
 *
 * Created on January 5, 2010, 5:11 PM
 */

#ifndef _GENOMICREGION_H
#define	_GENOMICREGION_H

#include "utility/types.h"
#include <string>


namespace Paris {
class GenomicRegion {
public:
	GenomicRegion();
	GenomicRegion(uint ID, const char *chr, uint beg, uint end);
	GenomicRegion(const GenomicRegion& other);
	virtual ~GenomicRegion();
	uint id;
	std::string _chromosome;					///< Which chromosome this belongs to
	uint _begin;									///< First base pair location
	uint _end;										///< last base pair location



	static float pvThreshold;					///< Threshold to be considered "statistically interesting"

protected:
};
}
#endif	/* _GENOMICREGION_H */

