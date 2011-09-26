//
// C++ Interface: array2d
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYARRAY2D_H
#define UTILITYARRAY2D_H
#include <iostream>
#include <iomanip>

namespace Utility {

/**
	@brief Provides simple wrapping for easy/efficient indexing into 2 dimensional array
	@Note Based on example from http://stackoverflow.com/questions/340943/c-multi-dimensional-arrays-on-the-heap
	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
template <typename T>
class Array2D{
public:
	Array2D() : cols(0), rows(0), data(NULL) { }
    Array2D(const int rows, const int cols) : cols(cols), rows(rows) { data = new T[cols*rows]; memset((void*)data, 0, sizeof(T)*rows*cols); }
    ~Array2D() { delete[] data; }
	T& operator()(const int r, const int c) { return data[r*cols + c]; } 
	Array2D& operator=(const Array2D& other) { if (data) delete[] data; cols=other.cols; rows=other.rows; data=new T[cols*rows]; memcpy((void*)data, (void*)other.data, rows*cols*sizeof(T)); }

	void Print(std::ostream& os)  {
		int width = 8; 
		os<<std::setw(width)<<" ";
		for (int i=0; i<cols; i++) 
			os<<std::setw(width)<<i+1;
		os<<"\n";
		for (int i=0; i<rows; i++) {
			os<<std::setw(width)<<i+1;
			int marg = 0;
			for (int c=0; c<cols; c++) {
				marg+=(*this)(i, c); 
				os<<std::setw(width)<<(*this)(i, c);
			}
			os<<setw(width)<<marg<<"\n";
		}
}
protected:
	int cols;					///<The width of a single row
	int rows;					///<rowcount
	T *data;					///<The actual data
};


}

#endif
