//
// C++ Interface: appinterface
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIAPPINTERFACE_H
#define GENOMESIM_GUIAPPINTERFACE_H
#include "appcontroller.h"
#include "wxuser.h"
#include <wx/wx.h>
namespace GenomeSIM {

namespace GUI {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/



class AppInterface : public WxUser {
public:
	
	struct Version {
		int maj;
		int min;
		int bug;
		int build;
		string date;

		Version();
		string GetVersion() {
			stringstream ss;
			ss<<"wxGenomeSIMLA "<<maj<<"."<<min<<"."<<bug<<"a ("<<build<<") ";
			return ss.str();
		}
	};

	AppInterface() : appController(NULL) { }
	virtual ~AppInterface() {}

	virtual void InitAppController(AppController *appController) = 0;

	/**
	 * @brief This is called prior to saving the configuration
	 */
	virtual void Commit() = 0;

	/**
	 * @brief This is called just after a configuration has been loaded
	 */
	virtual void RefreshSettings() = 0;

	/**
	 * @brief Compares current state with that of the saved version at least against what is in memory
	 * @return True indicates that something has changed 
	 */
	virtual bool HasChanged() = 0;

	


	//void SetAppController(AppController *controller) { appController = controller; }
protected:
	AppController *appController;
	Version version;
};



}

}

#endif
