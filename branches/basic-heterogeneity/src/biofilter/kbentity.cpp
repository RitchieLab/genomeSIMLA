//
// C++ Implementation: kbentity.cpp
//
// Description: Basic functionality for the knowledge based entities
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kbentity.h"

using namespace std;

namespace Biofilter {
namespace Knowledge {

SnpManager *KbEntity::snpManager = NULL;

KbEntity::KbEntity()
	:dbID(0), name(""), alias(""), desc(""), processed(false) { }

KbEntity::KbEntity(uint dbID, const string& name, const string& desc)
	: dbID(dbID), name(name), alias(""), desc(desc), processed(false) { }

KbEntity::~KbEntity() { }

uint KbEntity::DbID() {
	return dbID;
}
string KbEntity::Name() { return name; }
void KbEntity::Name(const char *newName) {
	name = newName;
}
std::set<string> KbEntity::Aliases() {
	return aliases;
}

void KbEntity::AddAlias(const char *alias) {
	aliases.insert(alias);
	if (this->alias.length() == 0 && alias[0] >= 'A') {
		string a = alias;
		if (a.find("LOC") == string::npos)
			SetAlias(alias);
	}
}

void KbEntity::SetAlias(const char *alias) {
	if (this->alias.length() == 0 )
		this->alias = alias;
}

string KbEntity::CommonName() {
	if (alias.length() > 0)
		return alias;
	else {
		return name;
	}
}

string KbEntity::Desc() {	return desc; }

void KbEntity::Desc(const char *desc) {
	this->desc = desc;
}

void KbEntity::PrintTabs(int tabCount, ostream& os) {
	for (int i=0; i<tabCount; i++)
		os<<"\t";
}

void KbEntity::MarkAsProcessed(bool isProcessed) {
	processed = isProcessed;
}


}
}
