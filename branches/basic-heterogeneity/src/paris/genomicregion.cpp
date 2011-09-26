/* 
 * File:   genomicregion.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 5:11 PM
 */

#include "genomicregion.h"
#include "chromosome.h"

namespace Paris {
float GenomicRegion::pvThreshold = 0.05;

GenomicRegion::GenomicRegion() : id(0), _chromosome(0), _begin(0), _end(0) {}
GenomicRegion::GenomicRegion(uint id, const char *chr, uint beg, uint end) :id(id),  _chromosome(chr), _begin(beg), _end(end) { }
GenomicRegion::GenomicRegion(const GenomicRegion& other) : id(other.id), _chromosome(other._chromosome), _begin(other._begin), _end(other._end) { }

GenomicRegion::~GenomicRegion() { }


}
