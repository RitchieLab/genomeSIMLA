//Stringmanip.h

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


///////////////////////////////////////////////////////////////////// 
//
// Class offering functions for handling strings and converting
// data into/from strings.
//
/////////////////////////////////////////////////////////////////////

#ifndef __STRINGMANIP_H__
#define __STRINGMANIP_H__

#include <sstream>
#include <vector>

namespace SimPen {



using namespace std;

class Stringmanip{
  
  public:
    Stringmanip(); 

    static bool isnumber(string str);
    static int split_string(const string& input, 
      string delimiter, vector<string>& results);
    static int stoi(string number); 
    static unsigned stouint(string number);
    static int stodata(string number);
    static string itos(long number);
    static string itos(float number);
    static string itos(int number);
    static string itos(unsigned int number);
    static string itos(double number);
    static double stodouble(string number);
    
};

}

#endif   
