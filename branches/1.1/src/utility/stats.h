#include <vector>


template<class T>
class Accum {
public:
	std::vector<T> data;
	T sum;

	Accum(int size = 0) : sum(0.0) {
		data.reserve(size);
	}

	void AddValue(T val) {
		data.push_back(val);
		sum+=val;
	}

	float GetMean() {
		return (float)sum/(float)data.size();
	}

	T Sum() {
		return sum;
	}

	float GetStdDev() {
		typename std::vector<T>::iterator itr = data.begin();
		typename std::vector<T>::iterator end = data.end();
		
		float mean = GetMean();
		
		float vsum = 0;	
		while (itr != end) {	
			T variance = *itr++ - mean;
			
			vsum += variance*variance;
		}

		return sqrt(vsum/(float)data.size());
	}
};

