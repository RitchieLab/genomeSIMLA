/* 
 * File:   parislogs.h
 * Author: torstees
 *
 * This is just a single place to keep up with the various files. It makes
 * it easier to find filenames and whatnot
 *
 * Created on January 19, 2010, 12:40 PM
 */

#ifndef PARIS_PARISLOGS_H
#define	PARIS_PARISLOGS_H
#include "utility/strings.h"
#include <fstream>
namespace Paris {

using namespace Utility;

struct FileLogs {
	void Open(const char *prefix) {
		reportPrefix = prefix;

		kbReport.open(GetFilename("KB", "txt").c_str());

	}

	void OpenBinDetails() {
		binDetails.open(GetFilename("bindetails", "txt").c_str());
	}

	void OpenFinalReport(uint idx = 0) {
		if (idx == 0)
			finalReport.open(GetFilename("output", "txt").c_str());
		else {
			if (idx > 1)
				finalReport.close();
			finalReport.open(GetFilename("output", "txt", idx).c_str());
		}
	}

	void OpenFeatureDetails() {
		if (WriteFeatureDetails) {
			featureDetails.open(GetFilename("-features", "txt").c_str());
			emptyFeatureDetails.open(GetFilename("-empty-features", "txt").c_str());
		}
	}

	void OpenDetailedFeatureReport() {
		if (WriteDetailedFeatureReport) {
			detailedFeatureReport.open(GetFilename("-detailed-feature", "txt").c_str());
		}
	}

	std::string GetFilename(const char *suffix, const char *extention, uint idx = 0) {
		std::string fnSuffix = suffix;
		size_t posCol = fnSuffix.find(':', 0);
		
		while (posCol != std::string::npos) {
			fnSuffix = fnSuffix.replace(posCol, 1, "_");
			posCol			= fnSuffix.find(":", 0);
		}

		std::string iteration = "";
		if (idx > 0)
			iteration = "-" + ToString(idx);
		return reportPrefix + "-" + fnSuffix + iteration + "." + extention;
	}

	//static bool WriteBinDetails;
	static bool WriteFeatureDetails;
	static bool WriteDetailedFeatureReport;

	std::ofstream binDetails;
	std::ofstream finalReport;
	std::ofstream kbReport;
	std::ofstream featureDetails;
	std::ofstream emptyFeatureDetails;
	std::ofstream detailedFeatureReport;

	std::string reportPrefix;

	static FileLogs logger;
};
}

#endif	/* PARIS_PARISLOGS_H */

