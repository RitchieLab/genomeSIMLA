/* 
 * File:   bin.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 12:15 PM
 */

#include "bin.h"
#include <iostream>
namespace Paris {
Utility::Random &Bin::rnd = Utility::Random::globalGenerator;

Bin::Bin() : binIndex((uint)-1) {
}

Bin::Bin(const Bin& other) : binIndex(other.binIndex), features(other.features) {
}

Bin::Bin(uint binIndex) : binIndex(binIndex)  {}

Bin::~Bin() {
}

uint Bin::AddFeature(Feature* feature) {
	feature->BinIndex(binIndex);
	features.push_back(feature);
	return features.size();
}

void Bin::Draw(uint count, std::set<Feature*>& ignore, std::set<Feature*>& features) {
	std::set<Feature*> localFeatures;

	std::set<Feature*>::iterator notIgnored = ignore.end();
	int featureCount = this->features.size();

	while (localFeatures.size() < count) {
		uint fID = rnd(featureCount);
		Feature *f = this->features[fID];
		if (ignore.find(f) == notIgnored)
			localFeatures.insert(f);
	}
	features.insert(localFeatures.begin(), localFeatures.end());
}


#ifdef TEST_APP
TEST(BinTesting, AddFeatureToBin) {
	Bin bin(5);

	Feature f1(1, "1", 1, 100);
	Feature f2(2, "1", 101, 200);

	EXPECT_EQ(1, bin.AddFeature(&f1));
	EXPECT_EQ(2, bin.AddFeature(&f2));
	EXPECT_EQ(5, f1.BinIndex());
	EXPECT_EQ(5, f1.BinIndex());
}

TEST(BinTesting, DrawFeatureFromBin) {
	Bin bin(20);
	
	Feature f1(1, "1", 1, 100);
	Feature f2(2, "1", 101, 200);
	Feature f3(3, "1", 201, 300);

	bin.AddFeature(&f1);
	bin.AddFeature(&f2);

	std::set<Feature*> used;
	used.insert(&f1);

	std::set<Feature*> drawn;
	bin.Draw(1, used, drawn);
	ASSERT_EQ(1, drawn.size());
	EXPECT_EQ(&f2, *(drawn.begin()));

	//Make sure that we really are using the ignored list
	bin.AddFeature(&f3);
	used.insert(&f2);
	drawn.clear();
	bin.Draw(1, used, drawn);
	ASSERT_EQ(1, drawn.size());
	EXPECT_EQ(&f3, *(drawn.begin()));
	
	drawn = used;
	bin.Draw(1, used, drawn);
	ASSERT_EQ(3, drawn.size());


}

#endif //TEST_APP
}
