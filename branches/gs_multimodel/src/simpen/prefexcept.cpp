#include "prefexcept.h"

namespace SimPen {
// constructor
PrefExcept::PrefExcept(){
        error = "No message set by code.";
}

// constructor
PrefExcept::PrefExcept(string message){
        error = message;
}

// destructor
PrefExcept::~PrefExcept(){
        error = "";
}

string PrefExcept::get_error(){
        return error;
}

ostream & operator << (ostream & os, const PrefExcept & ke){
        os << ke.error;
        return os;
}

}

