//
// C++ Implementation: genomesimlawrapper
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genomesimlawrapper.h"
#include <sstream>
#include <iomanip>

using namespace std;

uint GenomeSimlaWrapper::curIdx = 0;
bool GenomeSimlaWrapper::UseAdamEve = false;

int main(int argc, char **argv) {
	if (argc < 2) {
		cout<<"genomeSIMLA-Wrapper threadCount params\n";
		abort();
	}
	
	GenomeSimlaWrapper gs(atoi(argv[1]));
	gs.LoadParameters(argv[2]);
	gs.Start();
	return 0;
}

params *GenomeSimlaWrapper::GetNextParam() {
	params *nextParam = NULL;
	LOCK_PARAMLIST;
	if (curIdx < parameterSettings.size())
		nextParam = &parameterSettings[curIdx++];
	cout<<"Working on "<<curIdx<<"\n";
	UNLOCK_PARAMLIST;
	return nextParam;
}

GenomeSimlaWrapper::GenomeSimlaWrapper(uint threadCount) : threadCount(threadCount)	{
/*
	ofstream log("log.txt");
	log<<"Thread Count: "<<threadCount<<"\n";
/*
	int idx = 1;
	for (uint a = 0; a<initPop.size(); a++) 
		for (uint b=0; b<carryCap.size(); b++)
			for (uint c=0; c<M.size(); c++)
				for (uint d=0; d<growth.size(); d++)
					for (uint e=0;e<t.size(); e++) {
						log<<"param"<<idx<<" "<<initPop[a]<<" "<<carryCap[b]<<" "<<M[c]<<" "<<growth[d]<<" "<<t[e]<<"\n";
						parameterSettings.push_back(params(initPop[a], carryCap[b], M[c], growth[d], t[e], idx++));
					}
*/
}

void *GenomeSimlaWrapper::RunThread(void *args) {
	GenomeSimlaWrapper *wrp = (GenomeSimlaWrapper *)args;
	
	params *settings = wrp->GetNextParam();
	while (settings) {
		settings->RunGenomeSIMLA();
		settings = wrp->GetNextParam();
	}
	return NULL;
}

void GenomeSimlaWrapper::Start() {
	pthread_t threads[threadCount];

	cout<<"Threads to be launched: "<<threadCount<<"\n";	
	for (uint i = 1; i<threadCount; i++) {
		cout<<"Launching #"<<i<<"\n";
		pthread_create(&threads[i], NULL, RunThread, (void*)this);
	}
	
	RunThread((void*)this);
	
	for (uint i = 1; i<threadCount; i++) {
		pthread_join(threads[i], NULL);
	}
	
	cout<<"\n\n-------------------------\n-- Done!\n";
	
}

GenomeSimlaWrapper::~GenomeSimlaWrapper()
{
}

/*
INITPOP  500 	750 	1000
CARRYCAP 120000 500000 	900000
M        305 	315 	325 	335 	345 	355
GROWTH   0.005 	0.0075 	0.01
T        0.1 	0.2 	0.3
*/


void GenomeSimlaWrapper::LoadParameters(const char *filename) {

	fstream prm(filename, ios::in);
	
	string temp;
	string var;
	char line[4096];
	while (!prm.eof()) {
		line[0] = '\0';
		prm.getline(line, 4096);
		
		stringstream ss(line);
		ss>>var;

		if (var == "INITPOP") {
			while (!ss.eof()) {
				ss>>temp;
				
				initPop.push_back(atoi(temp.c_str()));
			}
		}
		else if (var == "CARRYCAP") {
			while (!ss.eof() ){
				ss>>temp;
				carryCap.push_back(atoi(temp.c_str()));
			}
		}
		else if (var=="M") {
			while (!ss.eof()) {
				ss>>temp;
				M.push_back(atoi(temp.c_str()));		
			}
		}
		else if (var=="GROWTH") {
			float temp;
			
			while (!ss.eof()) {
				ss>>temp;
				growth.push_back(temp);
			}
		}
		else if (var=="T") {
			float temp;
			while (!ss.eof()) {
				temp = 0.0;
				ss>>temp;
				if (temp > 0.0) 
					t.push_back(temp);
			}
		}
		else if (var=="SEED") {
			uint temp;
			while (!ss.eof() ) {
				temp = 0;
				ss>>temp;
				if (temp>0)
					seeds.push_back(temp);
			}
		}
		else if (var=="OFFSPRING") {
			uint f1, f2;
			while (!ss.eof() ) {
				f1 = f2 = 0;

				ss>>f1>>f2;
				if (f1 > 0 && f2 > 0) {	
					char offspr[32];
					sprintf(offspr, "%d %d", f1, f2);
					offspring.push_back(offspr);
				}
			}
		}
		else if (var == "BW") {
			float f1, f2;
			while (!ss.eof()) {
				f1 = f2 = 0.0;
				ss>>f1>>f2;
					
				if (f1 + f2 > 0.0) {
					char offspr[32];
					sprintf(offspr, "%.10f %.10f", f1, f2);
					cout<<"\nBW "<<offspr<<"\n";
					bw.push_back(offspr);
				}
			}
		}
		else if (var == "USE_ADAM_EVE") {
			UseAdamEve = true;
		}
		else if (var == "FOUNDER_COUNT") {
			int temp;
			while (!ss.eof()) {
				ss>>temp;
				founderCount.push_back(temp);
			}
		}
		else if (var == "DISTORTION") {
			float temp;
			while (!ss.eof()){
				ss>>temp;
				if (ss.str() == "0.0" || temp>0.00)
					distortion.push_back(temp);
			}
		}
		else if (var == "REPEATS") {
			int temp;
			while (!ss.eof()) {
				ss>>temp;
				repeats.push_back(temp);
			}
		}
		else
			cout<<"Skipping line: "<<line<<" due to unkown list name\n";
	}

	char newFilenameYeah[4096];
	sprintf(newFilenameYeah, "%s-log.txt", filename);
    ofstream log(newFilenameYeah);
    log<<"Thread Count: "<<threadCount<<"\n";

	log<<setprecision(10);

    int idx = 1;
	cout<<"Init. Pop : "<<initPop.size()<<"\n";
	cout<<"Carry Cap : "<<carryCap.size()<<"\n";
	cout<<"M         : "<<M.size()<<"\n";
	cout<<"Growth    : "<<growth.size()<<"\n";
	cout<<"t         : "<<t.size()<<"\n";
	cout<<"Seeds	 : "<<seeds.size()<<"\n";
	cout<<"Offspring : "<<offspring.size()<<"\n";

	if (UseAdamEve) {
		cout<<"Founders  : "<<founderCount.size()<<"\n";
		cout<<"Distortion: "<<distortion.size()<<"\n";
		cout<<"Repeats   : "<<repeats.size()<<"\n";
	}
	if (bw.size() == 0) 
		bw.push_back("NO_BLOCKS");
	for (uint a = 0; a<initPop.size(); a++)
		for (uint b=0; b<carryCap.size(); b++)
			for (uint c=0; c<M.size(); c++)
				for (uint d=0; d<growth.size(); d++)
					for (uint e=0;e<t.size(); e++)
						for (uint f=0; f<seeds.size(); f++)  
							for (uint g=0; g<bw.size(); g++) 
								if (UseAdamEve) {
									for (uint h=0; h<founderCount.size(); h++) 
										for (uint i=0; i<distortion.size(); i++) 
											for (uint j=0; j<repeats.size(); j++) {
												log<<"param"<<idx<<" "<<initPop[a]<<" "<<carryCap[b]
													<<" "<<M[c]<<" "<<growth[d]<<" "<<t[e]<<" "
													<<seeds[f]<<" "<<bw[g]<<" "<<founderCount[h]
													<<" "<<distortion[i]<<" "<<repeats[j]<<"\n";
												parameterSettings.push_back(params(initPop[a], carryCap[b], M[c], growth[d], t[e], seeds[f], bw[g], idx++, founderCount[h], distortion[i], repeats[j]));
											}
								}
								else {
									log<<"param"<<idx<<" "<<initPop[a]<<" "<<carryCap[b]<<" "<<M[c]<<" "<<growth[d]<<" "<<t[e]<<" "<<seeds[f]<<" "<<bw[g]<<"\n";
									parameterSettings.push_back(params(initPop[a], carryCap[b], M[c], growth[d], t[e], seeds[f], bw[g], idx++));
								}

}


