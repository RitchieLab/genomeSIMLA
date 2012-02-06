#ifndef SIMULATION_H
#define SIMULATION_H
#include <iostream>
#include <iomanip>
namespace Simulation {


struct Probabilities {
	uint cAB;
	uint cAb;
	uint caB;
	uint cab;
	float mafA;
	float mafB;
	float sum;
	float _pab;
	float _paB;
	float _pAb;
	float _pAB;

	void CalcSum() { sum=(float)(cAB+cAb+caB+cab);_pAB=(float)cAB/sum;_pAb=(float)cAb/sum;_paB=(float)caB/sum;_pab=(float)cab/sum;	}
	float pAB() { return _pAB; }
	float pAb() { return _pAb; }
	float paB() { return _paB; }
	float pab() { return _pab; }


	//This is debug info 
	uint P;
	uint Q;

	bool IsValid() {
		return cAB  && cAb && caB && cab;
	}

	void Report(ostream& os) {
		os<<"Report For Probabilities ("<<P<<"x"<<Q<<")\n";
		os<<"C(AB) "<<cAB<<" C(Ab) "<<cAb<<" C(aB) "<<caB<<" C(ab) "<<cab<<"\n";
		os<<"P(AB): "<<setprecision(2)<<pAB()<<" P(Ab) "<<setprecision(2)<<pAb()<<" P(aB) "<<setprecision(2)<<paB()<<" P(ab) "<<setprecision(2)<<pab()<<"\n";
	}


	Probabilities() : cAB(0), cAb(0), caB(0), cab(0), sum(0), P(0), Q(0) { }
};


}

#endif

