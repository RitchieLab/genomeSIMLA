//Stringmanip.cpp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////

#include "stringmanip.h"

namespace SimPen {



// Use:  Checks that the string passed in is a number
// Arg:  str -- string to check for number
// Ret:  true for number, false otherwise
bool Stringmanip::isnumber(string str){
  for(unsigned int i=0; i<str.length(); i++){
    if((str[i]<'0'||str[i]>'9') && str[i] != '.'){
      return false;
    }
  }
  return true;
}

// Use:  Split string into tokens as defined by the delimiter
// Arg:  input -- string to split
//       delimiter -- string to split on
//       results -- string vector that stores tokens
// Ret:  number of tokens
int Stringmanip::split_string(const string& input, 
  string delimiter, vector<string>& results){
  int iPos = 0;
  int newPos = -1;

  int sizeS2 = 1;
  int isize = input.size();

  vector<unsigned int> positions;
  newPos = input.find_first_of (delimiter, 0);
  if( newPos < 0 ) { return 0; }

  int numFound = 0;

  while( newPos >= iPos ){
    numFound++;
    positions.push_back(newPos);
    iPos = newPos;
    newPos = input.find_first_of (delimiter, iPos+sizeS2);
  }

  for(unsigned int i=0; i <= positions.size(); i++ ){
    string s;
    if( i == 0 ) { s = input.substr( i, positions[i] ); }
    int offset = positions[i-1] + sizeS2;
    if( offset < isize )
    {
      if( i == positions.size() )
      {
        s = input.substr(offset);
      }
      else if( i > 0 )
      {
        s = input.substr( positions[i-1] + sizeS2, 
          positions[i] - positions[i-1] - sizeS2 );
      }
    }
    if( s.size() > 0 )
    {
      results.push_back(s);
    }
  }

  return results.size();
}

// Use: Converts string number to int
// Arg: number -- number as string
// Ret: number as integer
int Stringmanip::stoi(string number){
    int value;
    stringstream ss(number);
    ss >> value;
    return value;
}

// Use: Converts string number to unsigned int
// Arg: number -- number as string
// Ret: number as unsigned integer
unsigned int Stringmanip::stouint(string number){
  unsigned int value;
  stringstream ss(number);
  ss >> value;
  return value;
}

// Use: Converts string number to int
// Arg: number -- number as string
// Ret: number as integer
int Stringmanip::stodata(string number){
  int value;
  stringstream ss(number);
  ss >> value;
  return value;
}

// Use: Converts string number to double
// Arg: number -- number as string
// Ret: number as double
double Stringmanip::stodouble(string number){
  double value;
  stringstream ss(number);
  ss >> value;
  return value;
}

// Use: Converts unsigned integer number to string
// Arg: number -- number as unsigned integer
// Ret: number as string
string Stringmanip::itos(unsigned int number){
  stringstream oss;
  oss << number;
  return oss.str();
}

// Use: Converts long number to string
// Arg: number -- number as long
// Ret: number as string
string Stringmanip::itos(long number){ 
  stringstream oss;
  oss << number;
  return oss.str();
}

// Use: Converts integer number to string
// Arg: number -- number as integer
// Ret: number as string
string Stringmanip::itos(int number){
  stringstream oss;
  oss << number;
  return oss.str();
}

// Use: Converts float number to string
// Arg: number -- number as float
// Ret: number as string
string Stringmanip::itos(float number){
  stringstream oss;
  oss << number;
  return oss.str();
}

// Use: Converts double number to string
// Arg: number -- number as double
// Ret: number as string
string Stringmanip::itos(double number){
  stringstream oss;
  oss << number;
  return oss.str();
}

}
