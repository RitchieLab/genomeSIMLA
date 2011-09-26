//
// C++ Interface: wxdiseaseconfiguration
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIWXDISEASECONFIGURATION_H
#define GENOMESIM_GUIWXDISEASECONFIGURATION_H
#include "wxuser.h"
#include <string>

namespace GenomeSIM {

namespace GUI {

/**
Base class for various disease configuration pages


	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class wxDiseaseConfiguration : public wxPanel, public WxUser {
public:
    wxDiseaseConfiguration() {}
    ~wxDiseaseConfiguration() {}

	virtual bool CanSave() { return true; }
	virtual bool CanSaveAs() { return true; }
	virtual bool CanImport() { return false; }
	virtual bool CanExport() { return false; }
	virtual bool CanEvaluate() { return true; }

	virtual bool Save(const char *filename, const char *description) = 0;
	virtual bool Import(const char *filename, std::string &description) { assert(0); }
	virtual bool Export(const char *filename, const char *desc) { assert(0); } 
	virtual bool Evaluate() = 0;

	/**
	 * @brief Returns the number of loci associated with the model
	 * @return model size or -1 if the model can accept any number of loci (simpen)
	 */
	virtual int GetLocusCount() { return -1; };

	virtual std::string GetImportFileTypes() { return "File (*.*)|*.*"; }
	virtual std::string GetPreferredExtension()=0;

};

}

}

#endif
