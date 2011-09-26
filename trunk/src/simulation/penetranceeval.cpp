//
// C++ Implementation: penetranceeval
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "penetranceeval.h"
#include <iomanip>
#include <iostream>
namespace Simulation {

/*PenetranceEval::PenCell PenetranceEval::Allele::operator+(const Locus& other) {
	string penID;
	penID.Append(id);
	penID.Append('x');
	penID.Append(other.id);
	return PenCell(penID.c_str(), freq*other.freq);
}

PenetranceEval::Allele& PenetranceEval::Allele::operator=(const PenetranceEval::Allele& other) {
	id = other.id;
	freq = other.freq;
	return &this;
}*/

string PenetranceEval::Locus::GetGenotypeLabel(int gt) {
	string genotype;
	if (gt == 0) {
		genotype.append(1, al1.id);
		genotype.append(1, al1.id);
	}
	else if (gt==1){
		genotype.append(1, al1.id);
		genotype.append(1, al2.id);
	}
	else if (gt ==2){
		genotype.append(1, al2.id);
		genotype.append(1, al2.id);
	}
	return genotype;
}

double PenetranceEval::PenCell::GetPrevalence() {
	return frequency * penetrance;
}

double PenetranceEval::PenCell::GetPrevalence(string genotypeLabel) {
	double prev = 0;
	if (id.find(genotypeLabel) != string::npos)
		prev = GetPrevalence();
	return prev;
}

double PenetranceEval::PenCell::GetAntiPrevalence() {
	return frequency * (1.0 - penetrance);
}

double PenetranceEval::PenCell::GetAntiPrevalence(string genotypeLabel) {
	double prev = 0.0;
	if (id.find(genotypeLabel) != string::npos) 
		prev = GetAntiPrevalence();
	return prev;
}

PenetranceEval::PenCell PenetranceEval::PenCell::Append(const Allele& al1, const Allele& al2) const {
		PenCell newCell(*this);
		newCell.id.append(1, al1.id);
		newCell.id.append(1, al2.id);

		double freq = al1.freq * al2.freq;
		if (al1.freq != al2.freq)
			freq*=2;
		if (newCell.frequency < 0.000000000001) 
			newCell.frequency=freq;
		else
			newCell.frequency*=freq;
		return newCell;
}

void PenetranceEval::NmiCalculator::Report(ostream& os) {
	int cellWidth = 12;
	os<<"<TABLE><TR><TH> </TH><TH>Cases</TH><TH>Controls</TH>"
		"<TH>C+C</TH><TH>Entropy</TH><TH>u</TH><TH>C of AB</TH></TR>\n";
	
	vector<NmiCell>::iterator itr=cells.begin();
	vector<NmiCell>::iterator end=cells.end();
	double sumAB = 0.0;
	while (itr!=end) {
		os<<"<TR><TH>"<<itr->label<<"</TH>"
			<<"<TH>"<<itr->cases<<"</TH>"
			<<"<TH>"<<itr->controls<<"</TH>"
			<<"<TH>"<<itr->GetSum()<<"</TH>"
			<<"<TH>"<<itr->GetEntropy()<<"</TH>"
			<<"<TH>"<<itr->GetT()<<"</TH>"
			<<"<TH>"<<itr->GetU()<<"</TH>"
			<<"<TH>"<<itr->GetCofAB()<<"</TH></TR>\n";
		sumAB+=itr->GetCofAB();
		itr++;			
	}
	os<<"<TR><TH>"<<sumCases<<"</TH><TH>"<<sumControls<<"</TH></TR>\n";
	os<<"<TR><TH>H(C)</TH><TH>"<<((sumCases*log(sumCases))+(sumControls*log(sumControls)))<<"</TH></TR>\n";
	os<<"<TR><TH>H(C|AB)</TH><TH>"<<sumAB<<"</TH></TR>\n";
	os<<"<TR><TH>NMI</TH><TH>"<<setw(cellWidth)<<NMI()<<"</TH></TR>\n</TABLE>\n";
	
}

double PenetranceEval::NmiCalculator::Sum() {
	double sum=0.0;
	
	vector<NmiCell>::iterator itr=cells.begin();
	vector<NmiCell>::iterator end=cells.end();
	
	cout<<setw(12)<<"Label"<<setw(12)<<"Cases"<<setw(12)<<"Controls"<<setw(12)<<"C+C"<<setw(12)<<"Entropy"<<setw(12)<<"T"<<setw(12)<<"u"<<setw(12)<<"C of AB"<<"\n";
	while (itr!=end) {
		sum+=itr->GetCofAB();	
		cout<<setw(12)<<itr->label
			<<setw(12)<<itr->cases
			<<setw(12)<<itr->controls
			<<setw(12)<<itr->GetSum()
			<<setw(12)<<itr->GetEntropy()
			<<setw(12)<<itr->GetT()
			<<setw(12)<<itr->GetU()
			<<setw(12)<<itr->GetCofAB()<<"\n";
		itr++;
	}
	return sum;
}

double PenetranceEval::NmiCalculator::NMI() {
	double ofC = (sumCases*log(sumCases))+(sumControls*log(sumControls));
	double CofAB=Sum();
	return 	(ofC - CofAB) / ofC;
}

void PenetranceEval::AddPenCell(const PenCell &curCell, vector<Locus>::iterator itr) {
	if (itr != loci.end()) {
		PenetranceEval::PenCell c1 = curCell.Append(itr->al1, itr->al1);
		PenetranceEval::PenCell c2 = curCell.Append(itr->al1, itr->al2);
		PenetranceEval::PenCell c3 = curCell.Append(itr->al2, itr->al2);
		AddPenCell(c1, ++itr);
		AddPenCell(c2, itr);
		AddPenCell(c3, itr);
	}
	else	 
		//At this point, we have the penetrance cell we want to add
	 	penCells[curCell.id] = (curCell);
}

int PenetranceEval::InitializeCells() {
	penCells.clear();
	PenCell emptyCell;
	AddPenCell(emptyCell, loci.begin());
	return penCells.size();
}

void PenetranceEval::AddLocus(const double& freq1, const double& freq2, char id1, char id2) {
	loci.push_back(Locus(freq1, freq2, id1, id2));
}

void PenetranceEval::Clear() {
	loci.clear();
	penCells.clear();
}

PenetranceEval::PenCell& PenetranceEval::GetCell(const char *label) {
	return penCells[label];
}


void PenetranceEval::Normalize(map<string, double>& prev, const double& overall){
	map<string, double>::iterator itr = prev.begin();
	map<string, double>::iterator end = prev.end();

	while (itr != end) {
		if (itr->second > 0.00000) 
			itr->second/=overall;
		itr++;
	}
} 

void PenetranceEval::DumpPrevalences(ostream& os) {
	map<string, PenCell>::iterator itr = penCells.begin();
	map<string, PenCell>::iterator end = penCells.end();
	os<<"\tID  \tFreq\tPenetrance\n";
	while (itr != end) {
		os<<"\t"<<itr->second.id<<"\t"<<itr->second.frequency<<"\t"<<itr->second.penetrance<<"\n";
		itr++;
	}
}


void PenetranceEval::EvaluateSingleLocus(ostream& os, const double& prev, const double& anti) {
	const int AA=0, AB=1, BB=2;
	int locCount = loci.size();
	os<<"<P><HR><P><H2>Single Locus Effects</H2>\n";


	for (int i=0; i<locCount; i++) {
		//a=cases(AA),b=cases(AB), c=cases(BB)
		//d=controls(AA),e=controls(AB),f=controls(BB)
		double truePositives=0.0, trueNegatives=0.0, falsePositives=0.0, falseNegatives=0.0;
		double a=0.0, b=0.0, c=0.0, d=0.0, e=0.0, f=0.0;
		Locus &l = loci[i];
		os<<"<P><H3><U>Locus "<<l.al1.id<<"</U></H3>\n";

		string sAA=l.GetGenotypeLabel(AA);
		string sAB=l.GetGenotypeLabel(AB);
		string sBB=l.GetGenotypeLabel(BB);
	
		map<string, PenCell>::iterator itr = penCells.begin();
		map<string, PenCell>::iterator end = penCells.end();
		
		while (itr != end) {
			a+=(itr->second.GetPrevalence(sAA)/prev);
			b+=(itr->second.GetPrevalence(sAB)/prev);
			c+=(itr->second.GetPrevalence(sBB)/prev);
			d+=(itr->second.GetAntiPrevalence(sAA)/anti);
			e+=(itr->second.GetAntiPrevalence(sAB)/anti);
			f+=(itr->second.GetAntiPrevalence(sBB)/anti);

			//Determine risk so we can calculate NMI
			if (itr->second.GetPrevalence()/prev > itr->second.GetAntiPrevalence()/anti) {
				truePositives+=(itr->second.GetPrevalence());
				falseNegatives+=(itr->second.GetAntiPrevalence());
			}
			else {
				falsePositives+=(itr->second.GetPrevalence());
				trueNegatives+=(itr->second.GetAntiPrevalence());
			}
			itr++;
		}	

		NmiCalculator nmi;
		nmi.Append(truePositives, falseNegatives, "High Risk");
		nmi.Append(falsePositives, trueNegatives, "Low Risk");
		//nmi.Report(os);
		
		
		os<<"<TABLE border='1'>\n<TR><TH><B>Genotype</B></TH><TH><B>Prevalence</B></TH><TH><B>Anti Prevalence</B></TH></TR>\n";
		os<<"<TR><TD><B>"<<sAA<<"</B></TD><TD>"<<a<<"</TD><TD>"<<d<<"</TD></TR>\n";
		os<<"<TR><TD><B>"<<sAB<<"</B></TD><TD>"<<b<<"</TD><TD>"<<e<<"</TD></TR>\n";
		os<<"<TR><TD><B>"<<sBB<<"</B></TD><TD>"<<c<<"</TD><TD>"<<f<<"</TD></TR>\n";
		os<<"</TABLE>\n";


		os<<"<P><P>\n";
		os<<"<LEFT><B>Odds Ratio:</B></LEFT>\n";
		os<<"<TABLE border='1'><TR><TH> </TH><TH COLSPAN=2><B>Alleles</B></TH></TR>\n";
		os<<"       <TR><TH> </TH><TH><B>"<<l.al1.id<<"</B></TH><TH><B>"<<l.al2.id<<"</B></TH></TR>\n";
		//Calculate the Odds Ratio
		//Dominant Model: (b+c)d / (e+f)a
		double domOR = (double)(((b+c)*d)/((e+f)*a));
		os<<"       <TR align='right'><TH><B>Dominant Odds Ratio</B></TH><TD>"<<(1.0/domOR)<<"</TD><TD>"<<domOR<<"</TD></TR>\n";

		//Recessive Model: c*(d+e)/f*(a+b)
		double recOR = (double)(((c*(d+e))/((a+b)*f)));
		os<<"       <TR align='right'><TH><B>Recessive Odds Ratio</B></TH><TD>"<<(1.0/recOR)<<"</TD><TD>"<<recOR<<"</TD></TR>\n";
		//Additive Model:  ((b+2c)*(2d+e))/((2d+e)(b+2c))
		double addOR = (double)(((b+2*c)*(2*d+e))/((2*a+b)*(e+2*f)));
		os<<"       <TR align='right'><TH><B>Additive Odds Ratio</B></TH><TD>"<<(1.0/addOR)<<"</TD><TD>"<<addOR<<"</TD></TR>\n";
		os<<"</TABLE>\n";
	
		os<<"<P align='left'><B>NMI(Status|"<<l.al1.id<<") = "<<nmi.NMI()<<"</B></P>\n";
		

		os<<"<P><HR>\n";
	}

	

}

double PenetranceEval::EvaluateMarginals(ostream &os) {
	os<<"<H3>Marginal Variance</H3><TABLE  border='1'>\n";
	os<<"<TR bgcolor=\"#dddddd\"><TH><B>Genotype</B></TH><TH><B>Variance</B></TH>\n";

	double sumMarg = 0.0;
	map<string, double> marginals;
	int locCount = loci.size();

	for (int i=0; i<locCount; i++) {
		Locus &l = loci[i];

		for (int gt=0; gt<3; gt++) {
			string genotype = l.GetGenotypeLabel(gt);
			double margPenetrance = 0.0; 

			map<string, PenCell>::iterator itr = penCells.begin();
			map<string, PenCell>::iterator end = penCells.end();

			while (itr != end) {
				margPenetrance+=itr->second.GetPrevalence(genotype);
				itr++;
			}	
			os<<"<TR bgcolor=\"#ffffff\"><TD>"<<genotype<<"</TD><TD>"<<margPenetrance<<"</TD></TR>\n";
			marginals[genotype]+=margPenetrance;
			
			sumMarg+=margPenetrance;
		}
	}

	os<<"</TABLE>\n";
	sumMarg/=(float)(locCount * 3);
	map<string, double>::iterator itr = marginals.begin();
	map<string, double>::iterator end = marginals.end();
	double variance = 0.0;
	while (itr != end) {
		variance+=((itr->second - sumMarg)*(itr->second-sumMarg));
		itr++;
	}
	return sqrt(variance/(float)(locCount * 3));
}

int PenetranceEval::Evaluate(ostream& os) {

	//First, let's get the pen. cell prevs & anti prev.
	map<string, double> prev;
	map<string, double> antiprev;
	map<string, double> nPrev;
	map<string, double> nAnti;

/*	vector<string> HRgenotypes;
	vector<string> LRgenotypes;
*/	vector<string> genotypes;


	double truePositives = 0.0, falsePositives = 0.0,
		falseNegatives = 0.0, trueNegatives = 0.0;

	double overallPrev = 0.0;
	double overallAnti = 0.0;
	os<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);


	os<<"<CENTER><H2>Multi-Locus Effects</H2>\n";
	double margVar = EvaluateMarginals(os);
	os<<"Marginal Variance Std Deviation: <B>"<<margVar<<"<B>\n<BR><BR>";

	map<string, PenCell>::iterator itr = penCells.begin();
	map<string, PenCell>::iterator end = penCells.end();
	os<<"\n<HR>\n<BR>";
	os<<"<H3>Risk Breakdown</H3><TABLE  border='1'>\n\t<TR bgcolor=\"#dddddd\"> <TH><B>Genotype</B></TH> <TH><B>Prevalence</B></TH> <TH><B>Anti Prevalence</B></TH> <TH><B>Risk</B></TH> <TH><B>Prevalence (Normalized)</B></TH> <TH><B>Anti Prevalence</B></TH></TR>\n";

	while (itr != end) {
		overallPrev+=itr->second.GetPrevalence();
		overallAnti+=itr->second.GetAntiPrevalence();
		itr++;
	}
	itr = penCells.begin();
	while (itr != end) {
		string gt = itr->first;
		double p = itr->second.GetPrevalence();
		double a = itr->second.GetAntiPrevalence();
		//The normalized versions of these two
		double np = 0.0, na = 0.0;
		prev[gt] = p;
		antiprev[gt] = a;

		if (p>0.0) 
			np = p/overallPrev;
		if (a>0.0)
			na = a/overallAnti;

		nPrev[gt] = np;
		nAnti[gt] = na;
	
		string risk;
		//Cell is HR
		if (np > na) {
			os<<"<TR bgcolor=\"#777777\">";
			truePositives+=np;
			falseNegatives+=na;
			risk="HR";
		}
		else	{
			os<<"<TR bgcolor=\"#ffffff\">";
			falsePositives+=np;
			trueNegatives+=na;
			risk="LR";
		}
		os<<"<TD>"<<gt<<"</TD> <TD>"<<p<<"</TD> <TD>"<<a<<"</TD> <TD>"<<risk<<"</TD> <TD>"<<np<<"</TD> <TD>"<<na<<"</TD><TR>\n";

		itr++;
	}	
	os<<"</TABLE>\n";


	os<<"<BR><HR><HR><H3>Odds Ratio (Overall)</H3><TABLE border='1'>\n\t<TR bgcolor=\"#dddddd\"><TH> </TH><TH><B>Cases</B></TH><TH><B>Controls</B></TH></TR>\n";
	
	os<<"\t<TR bgcolor=\"#ffffff\"><TH><B>High Risk</B></TH><TD>"<<truePositives<<"</TD><TD>"<<falseNegatives<<"</TD></TR>\n";
	os<<"\t<TR bgcolor=\"#ffffff\"><TH><B>Low Risk</B></TH><TD>"<<falsePositives<<"</TD><TD>"<<trueNegatives<<"</TD></TR>\n";
	double OR = (truePositives * trueNegatives) / (falsePositives * falseNegatives);

	os<<"\t<TR bgcolor=\"#ffffff\"><TH> </TH><TD><B>Odds Ratio:</B></TD><TD>"<<OR<<"</TD></TR>\n";
	os<<"</TABLE>\n";
	
	EvaluateNMI(os);
	EvaluateSingleLocus(os, overallPrev, overallAnti);

	return -1;
}

void PenetranceEval::EvaluateNMI(ostream& os) {
	map<string, PenCell>::iterator itr = penCells.begin();
	map<string, PenCell>::iterator end = penCells.end();
	NmiCalculator nmi;
	
	while (itr != end) {
		double p=itr->second.GetPrevalence();
		double a=itr->second.GetAntiPrevalence();
		nmi.Append(p, a, itr->first.c_str());
		itr++;
	}
	//nmi.Report(os);
	os<<"<P align='left'><B>NMI(Status|AB) = "<<nmi.NMI()<<"</B></P>\n";
}

double PenetranceEval::EvaluateRisk(map<string, double> &table, vector<string>& genotypes) {
	double risk = 0.0;

	size_t gtCount = genotypes.size();

	for (size_t i=0; i<gtCount; i++)  {
		risk+=table[genotypes[i]];
	}
	return risk;
}

}
