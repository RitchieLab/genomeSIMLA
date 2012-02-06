// Preferences.h

// Reads preferences file and sets the preferences

#ifndef __NODE_H__
#define __NODE_H__

#include <iostream>

class node{
public:
	node(unsigned int l=0, unsigned int g=0, int i = 0){
		locus = l;
		genotype_num = g;
		index = i;
	}
	
	~node() {	}

	
	unsigned int locus;
	unsigned int genotype_num;
	int index;
};

#endif
