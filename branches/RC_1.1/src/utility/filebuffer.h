//
// C++ Interface: filebuffer
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYFILEBUFFER_H
#define UTILITYFILEBUFFER_H
#include <iostream>
#include <fstream>
#include <cstdio>				//std::rename(source, dest)	file renaming
#include <string>
#include <assert.h>
#include <set>
#include <stdlib.h>
namespace Utility {

using namespace std;
bool rename(const char *source, const char *dest);
/**
	@brief paging system built to allow a single, ordered report to be stored in memory and paged in and out of memory. 
	@note This assumes that the Type, T, has a valid < operator

	** The following functions are required to be defined by the type
	 * operator<(const T& other) const 				(this is required for STL
	 * operator=(const T& other) 					(basic assignment)
	 * LoadBinary(ifstream& file)					(Load contents from binary)
	 * WriteBinary(ofstream& file)					(write contents to binary)
	 * MergeGroups(const T& other)					(merge contents of one to another. This assumes 
													 that we are populating a massive structure in 
													 piecemeal fashion
	

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
template<class T>
class FileBuffer{
public:
	/**
	 * @brief Instantiate a file buffer object
	 * @param label Used for filenames
	 * @param initBufferSize The default buffer size.
	 * @param maxBufferSize The maximum size the buffer grows before writing to out-cache
	 */
    FileBuffer(const char *label, size_t initBufferSize, size_t maxBufferSize);
    ~FileBuffer();

	/**
	 * @Brief effectively, this closes off the buffer, which triggers rotation and whatnot
	 */
	void Reset();
	void Roll();
	/**
	 * @Brief Add an object to the cache
	 */
	bool Insert(T& newObject);

	/**
	 * @brief writes the contents to a single file and purges buffer/deletes caches
	 */
	void Archive(const char *filename);

	size_t Count() { return lastCount; }

	void DisplayContents(ostream& os);

	void Close(const char *filename);

	template<class F>
	void Run(F fn)  {
		typename set<T>::iterator itr = buffer.begin();
		typename set<T>::iterator end = buffer.end();
		while (itr != end) {
			fn(*itr++);
		}
	}
protected:
	/**
	 * @Brief Add an object to the cache
	 */
	bool Insert(T& newObject, bool doForce);
	/**
	 * @Brief Performs first time use initialization. After this, use reset (this can be skipped by using reset instead)
	 */
	void Init();
	/**
	 * @brief Writes recordCount records to the out cache file and loads the next chunk from the infile
	 */
	bool FlushBuffer(size_t recordCount, bool doWrap = true);

	/**
	 * @brief Reads initBufferSize records from the insource
	 */
	bool LoadBuffer(bool doWrap = true);

	string oFilename;			///<Name of the outgoing cache
	string iFilename;			///<Name of the incoming cache
	ofstream *ofile;			///<outgoing cache (more recent than ifile)
	ifstream *ifile;			///<Cache from previous run
	set<T, less<T> > buffer;	///<Current buffer
	size_t initBufferSize;		///<The size of "chunks" allowed for the buffer when reading in a chunk
	size_t maxBufferSize;		///<The max size the buffer can grow during use
	size_t curCount;			///<current number of models since start of ifile
	size_t lastCount;			///<The last observed size (when closing ofile)

	//If these two are the same, then we are in an open frame, which has no reference file
	T frameBottom;				///<this represents the lower bound (not inclusive). Anything greater will be written to file
	T frameTop;					///<If we are reading from a file, this is the last entry from our current frame
};

inline
bool rename(const char *source, const char *dest) {
	std::string cmd = "mv " + string(source) + " " + string(dest);
	system(cmd.c_str());
	return true;
}

template<class T>
class FileBufferReader {
	FileBufferReader(size_t bufferSize);
	~FileBufferReader();
	

protected:
	ifstream file;
	size_t bufferSize;
};

template<class T>
inline 
FileBuffer<T>::FileBuffer(const char *label, size_t initBufferSize, size_t maxBufferSize) :
		ofile(NULL), ifile(NULL),
		buffer(buffer), initBufferSize(initBufferSize), maxBufferSize(maxBufferSize),
		curCount(0), lastCount(0) { 
	oFilename = string(label) + ".outbuff";
	iFilename = string(label) + ".inbuff";
	unlink(oFilename.c_str());
	unlink(iFilename.c_str());
}

template<class T>
inline
FileBuffer<T>::~FileBuffer() { 
	if (ifile) {
		ifile->close();
		unlink(iFilename.c_str());
		delete ifile;
	}
	if (ofile) {
		ofile->close();
		unlink(oFilename.c_str());
		delete ofile;
	}
}


template<class T>
inline
void FileBuffer<T>::Archive(const char *filename) {
	Close(filename);
}

template<class T>
inline
void FileBuffer<T>::Close(const char *filename) {
	while (FlushBuffer(initBufferSize, false)) { }
	ifile->close();
	delete ifile;
	ifile=NULL;

	ofile->close();
	delete ofile;
	ofile=NULL;
	unlink(filename);
	if (rename(oFilename.c_str(), filename)!=0) 
		cerr<<"Error renaming file."<<filename<<"\n";
	unlink(iFilename.c_str());
}

template<class T>
inline
void FileBuffer<T>::Init() {
	if (ifile != NULL) {
		ifile->close();
		delete ifile;
	}
	if (ofile != NULL) {
		ofile->close();
		delete ofile;
	}
	ifile = new ifstream(iFilename.c_str(), ios::binary);
//cerr<<iFilename<<"\n";
	ofile = new ofstream(oFilename.c_str(), ios::binary|ios::trunc);
//cerr<<oFilename<<"\n";
}

template<class T>
inline
void FileBuffer<T>::Reset() {
	//Basically skip to the end of the souce file so we can start all overa gain
	while (FlushBuffer(initBufferSize, false)) { }
cerr<<"^";cerr.flush();
	Roll();
}

template<class T>
inline 
void FileBuffer<T>::Roll() {
	rename(oFilename.c_str(), iFilename.c_str());
	Init();
	frameTop = T();
}

template<class T>
inline
bool FileBuffer<T>::LoadBuffer(bool doWrap) {
	bool moreToCome = true;
	assert(ifile);
	if (!(ifile->good()) || ifile->eof()) {
		if (doWrap) {
			Roll();
			moreToCome = false;
		} else {
			return false;
 		}
	}
	size_t recordCount = initBufferSize;
	frameBottom = frameTop;

	while (recordCount-- > 0 && !ifile->eof() && ifile->good()) {
		if (frameTop.LoadBinary(*ifile)) {
			Insert(frameTop, true);
		}
		else
			break;
	}
	return moreToCome;
}


template<class T>
inline
void FileBuffer<T>::DisplayContents(ostream& os) {
	//Reset();

	int idx = 0;	
	while (buffer.size() > 0) {
		typename set<T>::iterator itr = buffer.begin();
		typename set<T>::iterator end = buffer.end();

		while (itr != end) {
			os<<idx++<<"\t ";
			itr->Write(os);
			itr++;
		}
		FlushBuffer(initBufferSize, false);
	}
}

template<class T>
inline
bool FileBuffer<T>::FlushBuffer(size_t recordCount, bool doWrap) {
	if (ofile == NULL) {
		Init();
		recordCount = buffer.size();
	}

	assert(ofile->good());


	if (recordCount == 0 || recordCount > buffer.size())
		recordCount = buffer.size();
	typename set<T>::iterator itr = buffer.begin();
	typename set<T>::iterator end = buffer.end();
//frameBottom.Write(cerr);
//itr->Write(cerr);
	//Advance until we have passed the frameBottom
	while (itr != end && (*itr < frameBottom || *itr == frameBottom)) { itr++; }
	typename set<T>::iterator flush=itr;
	//At this point, itr should be the same as the lower frame. We don't want to purge this node, 
	//since we don't know if it's been updated and it belongs to the previous frame

	T curtop;
	if (!ifile->good()) {
		while (itr != end && recordCount-- > 0)	{
			itr->WriteBinary(*ofile);	
			curtop = *(itr++);
			curCount++;
		}
	} else {
	//From here, we just write until we have reached the desired count, or until we've written frameTop
		while (itr != end && (*itr < frameTop || *itr == frameTop || frameTop == frameBottom)) {
			itr->WriteBinary(*ofile);	
			curtop = *(itr++);
			curCount++;
		}
	}
ofile->flush();
	//We want to save the item at itr, but we want to advance it so that it gets deleted. This
	//node replaces frameBottom once we do LoadBuffer()
//	if (itr==end)
//		Roll();
//	else 
	if (flush != itr)
		buffer.erase(flush, itr);


	frameTop = curtop;

//cerr<<"Writing to file, "<<recordCount<<" records ("<<curCount<<")\n";
	//Lets get rid of all of the ones 
	return LoadBuffer(doWrap);
}

template<class T>
inline
bool FileBuffer<T>::Insert(T& newObject) {
	return Insert(newObject, false);
}

template<class T>
inline
bool FileBuffer<T>::Insert(T& newObject, bool doForce) {
	//In case the new entry is outside the paged region
	while (!doForce && buffer.size() >= maxBufferSize && FlushBuffer(initBufferSize, true )) {}

	//At this point, we assume that the buffer matches reasonably well, so we'll insert or update
	typename set<T>::iterator itr = buffer.find(newObject);

	//Check to see if we found a previous copy
	if (itr != buffer.end()) {
		newObject.MergeGroups(*itr);
		buffer.erase(itr);
	}
	if (!buffer.insert(newObject).second) {
		cerr<<"Error inserting in buffer of size("<<buffer.size()<<"): ";newObject.Write(cerr);
	}
/*
if (buffer.find(newObject) == buffer.end()) {
cerr<<"OK, this model went in, but it totally disappeared!!!!!\n";
newObject.Write(cerr);
}*/

	return true;
}



}

#endif
