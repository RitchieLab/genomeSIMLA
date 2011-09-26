/**
 *  Simple Inline Tools to Help With HTML Reporting
 *  Author	: Eric Torstenson (torstenson@chgr.mc.vanderbilt.edu)
 */
#include <string>
#include "utility/strings.h"
namespace Biofilter {
namespace HTML {

using namespace std;

/**
 * @brief constructs a URL to query ensembl for a snp
 */
std::string EnsemblSnpReference(uint rsID);

/**
 * @brief Build a URL reference to a SNP at Ensembl's website
 * @param rsID The integer portion of the RS Number
 */
std::string LinkSnpReference(uint rsID);

/**
 * @brief Constructs a URL to query Ensembl for a gene (by ensembl ID)
 */
std::string EnsemblGeneReference(const char *ensID);
/**
 * @Brief Build a URL reference to a gene at Ensembl's Website
 * @param ensID The Stable ID String 
 * @param alias (optional) What the user will see (defaults to ensID)
 */
std::string LinkGeneReference(const char *ensID, const char *alias = NULL);


inline
std::string EnsemblGeneReference(const char *ensID) {
	return "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + std::string(ensID);
}
inline
std::string EnsemblSnpReference(uint rsID) {
	return "http://www.ensembl.org/Homo_sapiens/Variation/Summary?source=dbSNP;v=rs" + Utility::ToString(rsID);
}

inline
std::string LinkSnpReference(uint rsID) {
	return "<A HREF='"+EnsemblSnpReference(rsID)+"'>"+Utility::ToString(rsID)+"</A>";
}

inline
std::string LinkGeneReference(const char *ensID, const char *alias) {
	if (alias == NULL)
		alias = ensID;
	return "<A HREF='" + EnsemblGeneReference(ensID) + "'>" + std::string(alias) + "</A>";
}

}
}
