//
// C++ Implementation: filebuffer
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "filebuffer.h"

#ifdef FB_TEST

#include <set>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>


using namespace std;
using namespace Utility;
class TestModel{
public:
    TestModel() : value(0) { }
	TestModel(size_t v, int g) : value(v) { groups.insert(g); }

    ~TestModel() { }

	/** The following functions are required to be defined by the type
	 * operator<(const T& other) const 				(this is required for STL
	 * operator=(const T& other) 					(basic assignment)
	 * LoadBinary(ifstream& file)					(Load contents from binary)
	 * WriteBinary(ofstream& file)					(write contents to binary)
	 * MergeGroups(const T& other)					(merge contents of one to another. This assumes 
													 that we are populating a massive structure in 
													 piecemeal fashion)
	 */
	void WriteBinary(std::ostream& file) const {
		file<<value<<"\t";
		set<size_t>::iterator itr = groups.begin();
		set<size_t>::iterator end = groups.end();
		while (itr != end) {
			file<<*itr++<<"\t";
		}
		file<<"\n";
	}
	bool LoadBinary(std::istream& file) {
		char line[4096];
		file.getline(line, 4096);
		stringstream ss(line);
		value = 0;
		groups.clear();
		ss>>value;
		while (!ss.eof()){ 
			string v= "";
			ss>>v;
			if (v!="")
				AddGroup(atoi(v.c_str()));
		}
		return value != 0;
	}	
	bool operator==(const TestModel& other) const { return value == other.value; }	
	void Write(std::ostream& os) const { WriteBinary(os); }
	bool operator>(const TestModel& other) const { return value > other.value; }
	bool operator<(const TestModel& other) const { return value < other.value; }
	void AddGroup(size_t group) { groups.insert(group); }

	void MergeGroups(const TestModel& other) { groups.insert(other.groups.begin(), other.groups.end()); }
protected:
	set<size_t> groups;
	size_t value;
};
int main(int argc, char** argv) {
	vector<size_t> numbers;
	FileBuffer<TestModel> models("test_models", 5000, 20000);
	int groups[] = { 4, 2, 5, 3, 1};
	for (int g=0; g<3; g++) {
		cerr<<"\nWorking on group ("<<groups[g]<<")\n";
		for (size_t i=0; i<1000000; i++)
			numbers.push_back(i+1);
		random_shuffle(numbers.begin(), numbers.end());
		vector<size_t>::iterator itr = numbers.begin();
		vector<size_t>::iterator end = numbers.end();
		while (itr != end) {
			TestModel t(*itr, groups[g]);
			models.Insert(t);
			itr++;
		}
		cerr<<"\n";
	}
	models.Reset();
	models.Close("Fresh_super");
}
#endif
