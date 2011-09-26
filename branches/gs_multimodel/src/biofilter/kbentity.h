//
// C++ Implementation: kbentity.h
//
// Description: Basic functionality for the knowledge based entities
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __KB_ENTITY_H
#define __KB_ENTITY_H

#include <string>
#include <set>
#include "utility/types.h"
#include <soci.h>
#include "snpmanager.h"

namespace soci {
namespace details {

template <>
struct exchange_traits<unsigned int>
{
        typedef basic_type_tag type_family;
        enum { x_type = x_integer };
};
}
}


namespace Biofilter {
namespace Knowledge {

/**
 * @brief Provides basic functionality for knowledge based entities
 */
class KbEntity {
public:
	KbEntity();
	KbEntity(uint, const std::string& name="", const std::string& desc="");///< Basic constructor
	virtual ~KbEntity();
	virtual uint DbID();				///< Return the key
	std::string Name();					///< Return the name
	void Name(const char *);			///< Write to the entities name
	std::string Desc();
	void Desc(const char *desc);		///< Description
	void AddAlias(const char *alias);	///< Adds an alias to the group
	void SetAlias(const char *alias);	///< This will override the real name
	std::string CommonName();		///< Returns the preferred name
	std::set<std::string> Aliases();	///< Returns all aliases
	
	/**
	 * @brief Prints n tabs...helps align text trees
	 */
	void PrintTabs(int tabCount, std::ostream& os);

	virtual void MarkAsProcessed(bool isProcessed = true);
	static SnpManager* snpManager;
protected:
	uint dbID;							///< This is the unique key from the database
	std::string name;					///< Official Name
	std::string alias;					///< Common name
	std::string desc;					///< Description
	std::set<std::string> aliases;
	bool processed;				///< indicates that they have all been processed

};



}
}

#endif //__KB_ENTITY_H
