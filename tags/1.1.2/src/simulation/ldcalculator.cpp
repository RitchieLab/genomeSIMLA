//
// C++ Implementation: ldcalculator
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldcalculator.h"

namespace Simulation {

namespace Visualization {

unsigned int maxSnpDistance = 500000;


//Implimentation for LDResult
LdResult::LdResult(): P(NULL), Q(NULL), dprime(0.0), rsquared(0.0), lod(0.0) { }

LdResult::LdResult(Locus* p, Locus* q, int cAB, int cAb, int caB, int cab) : P(p), Q(q)  { 
	float total = 1.0/(float)(cab+caB+cAb+cAB);
	float pab = (float)cab * total;
	float paB = (float)caB * total;
	float pAb = (float)cAb * total;
	float pAB = (float)cAB * total;
	
	//This shouldn't happen except during debug
	assert(CheckMarginals(p, q, pab, paB, pAb, pAB));
	dprime = DPrime(pab, paB, pAb, pAB);
	rsquared = RSquared(pab, paB, pAb, pAB);
	lod = LOD(cab, caB, cAb, cAB, pab, paB, pAb, pAB);

	assert(dprime >= 0 && dprime <= 1.01);
	assert(rsquared >= 0 && rsquared <= 1.01);
}

bool LdResult::CheckMarginals(Locus *p, Locus* q, float pab, float paB, float pAb, float pAB) {
	float PdotB=dotB;
	float Pdotb=dotb;
	float PdotA=dotA;
	float Pdota=dota;

	float marginal = pAB + pAb;
	if (marginal < PdotA - 0.00001 || marginal > PdotA + 0.00001) {
		cerr<<"Marginals aren't adding up for "<<p->GetID()<<"x"<<q->GetID()<<" :\n";
		cerr<<"P(A): "<<PdotA<<" != P(AB) "<<pAB<<" + P(Ab) "<<pAb<<"\n";
		cerr<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
		p->WriteMarkerInfo(cerr, 8);
		q->WriteMarkerInfo(cerr, 8);
		return false;
	}
	marginal = paB + pAB;
	if (marginal < PdotB - 0.00001 || marginal > PdotB + 0.00001) {
		cerr<<"Marginals aren't adding up for "<<p->GetID()<<"x"<<q->GetID()<<" :\n";
		cerr<<"P(B): "<<PdotB<<" != P(aB) "<<paB<<" + P(AB) "<<pAB<<"\n";
		cerr<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
		p->WriteMarkerInfo(cerr, 8);
		q->WriteMarkerInfo(cerr, 8);
		return false;
	}
	marginal = paB + pab;
	if (marginal < Pdota - 0.00001 || marginal > Pdota + 0.00001) {
		cerr<<"Marginals aren't adding up for "<<p->GetID()<<"x"<<q->GetID()<<" :\n";
		cerr<<"P(a): "<<Pdota<<" != P(aB) "<<paB<<" + P(ab) "<<pab<<"\n";
		cerr<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
		p->WriteMarkerInfo(cerr, 8);
		q->WriteMarkerInfo(cerr, 8);
		return false;
	}	
	marginal = pAb + pab;
	if (marginal < Pdotb - 0.00001 || marginal > Pdotb + 0.00001) {
		cerr<<"Marginals aren't adding up for "<<p->GetID()<<"x"<<q->GetID()<<" :\n";
		cerr<<"P(b): "<<Pdotb<<" != P(Ab) "<<pAb<<" + P(ab) "<<pab<<"\n";
		cerr<<"\nP(A): "<<PdotA<<" P(a): "<<Pdota<<", P(B) "<<PdotB<<" P(b) "<<Pdotb<<"\n";
		p->WriteMarkerInfo(cerr, 8);
		q->WriteMarkerInfo(cerr, 8);
		return false;
	}
	return true;
}
float LdResult::DPrime(float pab, float paB, float pAb, float pAB) {
	float dprime = 0.0;
	float delta=(pAB*pab) - (pAb*paB);
	float m=0.0;
	float sign=1.0;
	if (delta>0.0) {
		m=fmin(dotb*dotA, dota*dotB);
		if (m != 0.0)
			dprime=delta/m;
		else
			dprime=0.0;
	assert(dprime >= 0 && dprime <= 1.01);
	}
	else {
		m=fmin(dotA*dotB, dota*dotb);
		if (m != 0.0)
			dprime=delta/m*-1.0;
		else
			dprime=0.0;
	assert(dprime >= 0 && dprime <= 1.01);
	}
//assert(dprime > 0.0);

	return dprime;
}

float LdResult::RSquared(float pab, float paB, float pAb, float pAB) {
	float r=pAB*pab-pAb*paB;
	return r*r/(dotA*dota*dotB*dotb);
}


float LdResult::LOD(int cab, int caB, int cAb, int cAB, float pab, float paB, float pAb, float pAB) {
	float lod = 0.0;
	if (cAB>0)
		lod+=(cAB*log(pAB/(dotA*dotB)));
	if (cAb>0)
		lod+=(cAb*log(pAb/(dotA*dotb)));
	if (caB>0)
		lod+=(caB*log(paB/(dota*dotB)));
	if (cAB>0)
		lod+=(cab*log(pab/(dota*dotb)));

	static float log10 = 1/log(10.0);
	if (lod > 0.0)
		lod*=log10;

	return lod;
}




}

}
