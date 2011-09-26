/* 
 * File:   bin.h
 * Author: torstees
 *
 * Bins are key to the permutation process. Each bin will hold similar 
 *
 * Created on January 5, 2010, 12:15 PM
 */

#ifndef PARIS_BIN_H
#define	PARIS_BIN_H

#include <vector>
#include "feature.h"
#include "utility/random.h"
namespace Paris {

class Bin {
public:
	Bin();
	Bin(uint binIndex);
	Bin(const Bin& orig);
	virtual ~Bin();

	/**
	 * @brief Add a feature to the collection
	 * @param feature the feature being added to the collection
	 * @return number of features after the new one has been added
	 */
	uint AddFeature(Feature* feature);
	/**
	 * @brief Draw count features into features, preventing any from ignore from being drawn
	 * @param count The number of features to draw
	 * @param ignore set of features that we don't want to add
	 * @param features the set where these will be added
	 **/
	void Draw(uint count, std::set<Feature*>& ignore, std::set<Feature*>& features);

	static Utility::Random &rnd;
private:
	uint binIndex;											///< For lookup purposes
	std::vector<Feature*> features;					///< Array of features
};



}
#endif	/* _BIN_H */

