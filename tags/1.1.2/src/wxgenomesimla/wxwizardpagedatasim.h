#ifndef _WXWIZARDPAGEDATASUMMARY_H_
#define _WXWIZARDPAGEDATASUMMARY_H_

#include "appinterface.h"
#include "modellist.h"
#include <iostream>

namespace GenomeSIM {

namespace GUI {

class wxWizardPageDataSummary;
class wxWizardPageDataSim : public wxWizardPageSimple, public AppInterface  {
public:
	wxWizardPageDataSim() { }
	wxWizardPageDataSim(wxWizard *parent) : wxWizardPageSimple(parent) {}
	
	/**
	 * @brief called when we reach the summary page. 
	 * @note It should be restarted each time we arrive at the page, since we might
	 * change things from time to time
	 */
	virtual void PrepSummary(wxWizardPageDataSummary*summary)=0;

	/**
	 * @brief this is called ONLY the last minute before we execute.
	 * @note unlike PrepSummary, this is a one shot deal, and no backups
	 */
	virtual void PrepConfig(std::ostream& config)=0;

	void InitAppController(AppController *ctrl);

	void SetDiseaseLoci(vector<Locus *> loci) {
		this->loci = loci;
	}
	void SetDiseaseModel(ModelList::ModelFileItem &diseaseModel) { this->diseaseModel=diseaseModel; }
protected:
	vector<Locus *> loci;						///<The actual selected loci
	ModelList::ModelFileItem diseaseModel;		///<Model in use

};


class wxWizardPageDataSummary : public wxWizardPageDataSim {
public:
	wxWizardPageDataSummary() {}
	wxWizardPageDataSummary(wxWizard *parent) : wxWizardPageDataSim(parent) { }

	/**
	 * This doesn't really make sense for this one
	 */
	void PrepSummary(wxWizardPageDataSummary *data) { }	

	void PrepConfig(std::ostream& config) { }

	/**
	 * @brief Clear any pre-existing report and set the title. Should be first called
	 */
	virtual void Clear(const char *title) = 0;
	/**
	 * @brief Add a new header chunk to the summary report
	 */
	virtual void AddHeader(const char *txt) = 0;

	/**
	 * @brief Write a portion of the note that is displayed as summary
	 */
	virtual void AddNote(const char *txt) = 0;

	/**
	 * @brief Allow pieces to write segments to the configuration file
	 */
	virtual void AddConfigLine(const char *txt) = 0;
};

inline
void wxWizardPageDataSim::InitAppController(AppController *ctrl) {
	appController = ctrl;
	RefreshSettings();
}

}
}
#endif

