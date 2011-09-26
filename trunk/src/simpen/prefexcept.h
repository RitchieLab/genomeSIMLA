#ifndef __PrefExcept_H__
#define __PrefExcept_H__

#include <string>
#include <iostream>

namespace SimPen {

using namespace std;

class PrefExcept{
        public:
                PrefExcept();
                PrefExcept(string message);
                ~PrefExcept();
                string get_error();
                friend ostream & operator << (ostream & os, const PrefExcept & ke);
        private:
                string error;
};

}

#endif
