#include "simpen.h"
#include "prefexcept.h"
#include <iomanip>

using namespace std;
using namespace SimPen;


int main(int argc, char** argv)
{
	string filename, output;
	int seed = 0;

	if (argc > 1)
		filename = argv[1];

	if (argc > 2)
		output = argv[2];
	
	vector<DiseaseLocus> loci;
	if (argc > 3) {
		for (int i = 3; i+1<argc; ) {
			int chr = atoi(argv[i++]);
			int loc = atoi(argv[i++]);
			DiseaseLocus l(chr - 1, loc - 1, 0.0, 0.0);
			loci.push_back(l);
		}
		int lociCount = loci.size();
		DiseaseLocus diseaseLoci[lociCount];
		
		cout<<setw(45)<<"Filename: "<<filename<<"\n";
		cout<<setw(45)<<"Output: "<<output<<"\n";
		cout<<setw(45)<<"Model Size: "<<lociCount<<"\n";
		
		cout<<setw(45)<<"Model Loci: ";
		for (int i = 0; i<lociCount; i++) {
			diseaseLoci[i]=loci[i];
			if (i > 0)
				cout<<"x";
			cout<<diseaseLoci[i].chromosome + 1<<":"<<diseaseLoci[i].locusIdx + 1;
		}

		cout<<"\n";
/*  EST Since genomeSIMLA sets this up, this check doesn't apply...need to rethink if we ever create a standalone version of simPEN
		if (lociCount < 2) {
			cout<<"Unable to produce epistatic models with fewer than 2 loci\n";
			return 1;
		}
 */
		return run_simpen(filename.c_str(), output.c_str(), diseaseLoci, lociCount, seed, false);
	}
	else {
		cout<<"Usage: filename output chr loc [chr loc....]\n";
		return 1;
	}
}


