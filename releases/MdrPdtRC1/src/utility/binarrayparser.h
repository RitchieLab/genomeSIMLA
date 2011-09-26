//
// C++ Interface: genobinparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYBINARRAYPARSER_H
#define UTILITYBINARRAYPARSER_H
#include <assert.h>
namespace Utility {

/**
@brief Parses an integer array for compressed integers of a constant width size.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

template <class T>
class BinArrayParser
{
public:
	/**	
 	 * @brief Constructor used to set up the parser
	 * @param values[] This is the storage that is to be used to find the compressed values
	 * @param count The size of the storage array
	 * @param bytesize The size of a single byte is
	 */
    BinArrayParser(T values[], const int count, const int bytesize );

	/**
	 * @brief This is a faux iterator function to allow grabbing pieces out one at a time
	 * @note This does not have any intelligence in it, so do not have two users of the parser calling this
	 * function expecting to get full coverage!
	 */
	int GetNextValue();

    ~BinArrayParser();

	/**
	 * @brief Returns the number if storage elements is used to store the compressed data
	 */
	int GetCount() {return count;}
	
	/**
	 * @brief Returns the number of bits is used to encode a single value
	 */
	int GetFrameSize();

	/**
	 * @brief Return the size a single storage unit (32 if T is a 32 bit integer)
	 */
	int GetWordSize();
protected:
	/**
	 * @brief perform the actual frame capture from within a binary string
	 * @param val This is the array of storage values
	 * @param start This is the position where we want to read a frame of binary data
	 * @param framesize The number of bits we are interested in
	 * @param readlength The number of bits to actually be read in (This is used when a frame exists between 2 members of the val array)
	 */
	int GetBinWindow(T val, int start, int framesize, int readlength);
	void Reset();				///<Reset the <I>cursor</I> back to the start
	
	int count;					///<The size of the storage 

	int framesize;				///<The number of bits used to store a single value
	int wordsize;				///<The size of a storage unit
	int bitPosition;			///<Used to keep up with psuedo iteration
	int lastOffset;				///<Used for psuedo iteration
	T *values;					///<The storage

};


template <class T>
inline
int BinArrayParser<T>::GetFrameSize() { return framesize; }

template <class T>
inline
int BinArrayParser<T>::GetWordSize() { return wordsize; }

template <class T>
inline
BinArrayParser<T>::~BinArrayParser()	{}	

template <class T>
inline
BinArrayParser<T>::BinArrayParser(T values[], const int count, const int bytesize) : 
	count(count), 
	framesize(bytesize), 
	bitPosition(0),
	lastOffset(0),
	values(values) 
{
	wordsize=sizeof(T) * 8;
}

template <class T>
inline
int BinArrayParser<T>::GetBinWindow(T val, int start, int framesize, int readlength) {
	int rFallout=start+readlength;
				
	T v=val;	
	v=v<<start;
	v>>=start+wordsize-rFallout;
	if (start>0)
		v<<=framesize-readlength;
	return (int)v;
}

template <class T>
inline
void BinArrayParser<T>::Reset() {
	bitPosition=0;
	lastOffset=0;
}

template <class T>
inline
int BinArrayParser<T>::GetNextValue() {
	int value=0;
	int wordIdx=bitPosition/wordsize;
	T word=values[wordIdx];
	int startPoint=bitPosition%wordsize;
	if ( (startPoint+framesize) > wordsize) {
		assert(count>wordIdx+1);
		int length=wordsize-startPoint;
		//OK, this is part one of the value
		value=GetBinWindow(word, startPoint, framesize, length);
		//OK, at this point, we want to eat the remainder
		value|=GetBinWindow(values[wordIdx+1], 0, framesize, framesize-length);
	}
	else 
		value=GetBinWindow(word, startPoint, framesize, framesize);
	bitPosition+=framesize;
	return value;
}	


}

#endif
