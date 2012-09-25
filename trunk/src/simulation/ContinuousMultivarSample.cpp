/*
 * ContinuousMultivarSample.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: jrw32
 */

#include "ContinuousMultivarSample.h"
#include "locus.h"

#include "utility/exception.h"
#include "utility/random.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using boost::algorithm::split;
using boost::algorithm::join;
using boost::algorithm::is_any_of;
using boost::to_upper;

using boost::lexical_cast;
using boost::bad_lexical_cast;

using std::vector;
using std::ifstream;
using std::stringstream;

namespace Simulation {

ContinuousMultivarSample::ContinuousMultivarSample(float genoError,
		float phenoError, float missingData, int n_in, const string& e_fn,
		const string& c_fn, const string& name_in) :
		Sample(genoError, phenoError, missingData), effect_fn(e_fn), covar_fn(c_fn), name(name_in), n_indiv(n_in) {

	// 1st, open the covariance matrix and check to ensure it is formatted
	// properly (square, positive definite)
	vector<double> row_one;

	// Open the covariance
	ifstream covar_file(covar_fn.c_str());
	if (!covar_file.is_open()){
		std::cerr<<"WARNING: cannot find " << covar_fn <<", exiting.";
		//throw an exception here!
		throw Utility::Exception::FileIO(covar_fn.c_str());
	}

	string line_one;
	getline(covar_file, line_one);
	stringstream ss(line_one);

	while(ss.good()){
		double f;
		ss >> f;
		row_one.push_back(f);
	}

	// OK, now we know the size of our matrix
	n_effects = row_one.size();
	covar_matrix = gsl_matrix_alloc(n_effects,n_effects);

	// Go through the first row and put it into the matrix
	for(int i=0; i<n_effects; i++){
		gsl_matrix_set(covar_matrix, 0, i, row_one[i]);
	}

	// Now, add all the other rows
	int row = 1;
	int col = 0;
	while(covar_file.good()){
		string next_line;
		getline(covar_file, next_line);
		stringstream l(next_line);

		if(next_line.size() && row >= n_effects){
			//Throw an exception, b/c we're too big!
			throw Utility::Exception::General("Too Many Rows in covariance matrix");
		}

		while(next_line.size() && l.good()){

			if(col >= n_effects){
				// Throw exception
				throw Utility::Exception::General("Uneven number of columns in covariance matrix (too many)");
			}

			double f;
			l >> f;
			gsl_matrix_set(covar_matrix, row, col, f);
			++col;
		}

		if(col < n_effects){
			//Throw exception
			throw Utility::Exception::General("Uneven number of columns in covariance matrix (too few)");
		}
		++row;
	}

	// OK, now we've read the covariance matrix, we need to check it for
	// positive definiteness using the cholesky factorization.
	gsl_error_handler_t* err_h = gsl_set_error_handler_off();
	int err = gsl_linalg_cholesky_decomp(covar_matrix);
	// restore the old error handler
	gsl_set_error_handler(err_h);

	if (err == GSL_EDOM){
		// Raise not positive definite
		throw Utility::Exception::General("Covariance matrix invalid (not positive definite)");
	}else if(err != GSL_SUCCESS){
		// Raise unknown error
		throw Utility::Exception::General("Unknown error checking covariance matrix");
	}

	// Now, read the effect matrix
	ifstream eff_file(effect_fn.c_str());
	if (!eff_file.is_open()){
		std::cerr<<"WARNING: cannot find " << effect_fn <<", exiting.\n";
		//throw an exception here!
		throw Utility::Exception::FileIO(effect_fn.c_str());
	}
	vector<string> result;

	int line_no = 0;
	while(eff_file.good()){
		string next_line;
		result.clear();
		getline(eff_file, next_line);
		++line_no;

		split(result, next_line, is_any_of(" \n\t"), boost::token_compress_on);

		if(result.size() && result[0][0] != '#'){

			int eff_num;
			double eff_size;
			if(result.size() == 3){
				// SNP label format
				try{
					// NOTE: effect index is 1-based
					eff_num = lexical_cast<int>(result[0]);
					if(eff_num <= 0 || eff_num > n_effects){
						std::cerr << "WARNING: Effect index out of range on line " << line_no << std::endl;
					}else{
						eff_size = lexical_cast<double>(result[2]);
						effect_matrix_lbl[eff_num-1][result[1]] = eff_size;
					}
				}catch(bad_lexical_cast&){
					// Badly formatted file
					std::cerr << "WARNING: Badly formatted effect matrix on line " << line_no << std::endl;
				}
			}else if(result.size() == 4){
				// SNP index format
				int chr;
				int snp_idx;
				try{
					// NOTE: effect index is 1-based
					eff_num = lexical_cast<int>(result[0]);
					if(eff_num <= 0 || eff_num > n_effects){
						std::cerr << "WARNING: Effect index out of range on line " << line_no << std::endl;
					}else{
						eff_size = lexical_cast<double>(result[3]);
						chr = lexical_cast<int>(result[1]);
						snp_idx = lexical_cast<int>(result[2]);
						effect_matrix_idx[eff_num-1][std::make_pair(chr-1, snp_idx-1)] = eff_size;
					}
				}catch(bad_lexical_cast&){
					std::cerr << "WARNING: Badly formatted effect matrix on line " << line_no << std::endl;
				}
			}else if(next_line.size()){
				std::cerr << "WARNING: Wrong number of arguments in effect matrix on line " << line_no << std::endl;
				// Bad formatting!
			}
		}
	}
}

void ContinuousMultivarSample::BuildSample(PoolManager& pools, PenetranceModel* model){

	// First off, we don't need the Penetrance model, so issue a warning if we
	// have one
	if(model){
		std::cerr << "WARNING: Penetrance Models are ignored with the Continuous Multivariate Sample, ignoring\n";
	}

	// The gsl vector for the mean vector calculation
	gsl_vector* mean_vec = gsl_vector_calloc(n_effects);

	// The gsl vector to hold the std. normals in
	gsl_vector* norm_vec = gsl_vector_calloc(n_effects);

	int id = 0;
	for(int i=0; i<n_indiv; i++){
		Individual* newPerson = new Individual(id, id, pools.GetPoolCount());
		++id;
		pools.DrawIndividual(*newPerson);
		people.push_back(newPerson);

		// I think this is the way to draw a sample.  I saw it in BasicSample.
		PoolManager::Iterator itr = pools.GetIterator();
		ChromPool* pool = itr.GetNext();
		while(pool){
			pool->ReserveIndividual(newPerson);
			pool=itr.GetNext();
		}

		// OK, now do the sampling of the continuous Multivariate

		// Set the normal vector and mean vector to 0
		gsl_vector_set_zero(mean_vec);
		gsl_vector_set_zero(norm_vec);

		// First, we will need the mean vector calculated from the effect matrix

		// Walk through the effect matrix indexed by labels first
		map<int, map<string, double> >::const_iterator eff_lbl_itr = effect_matrix_lbl.begin();
		map<string, double>::const_iterator lbl_itr;
		Locus tmp_loc;
		while(eff_lbl_itr != effect_matrix_lbl.end()){
			// NOTE: I don't need to check bounds of the effect index; we did
			// that when we parsed the effect matrix
			lbl_itr = (*eff_lbl_itr).second.begin();
			while(lbl_itr != (*eff_lbl_itr).second.end()){
				bool success = pools.GetLocusReference((*lbl_itr).first.c_str(), tmp_loc);
				if(success){
					double init_val = gsl_vector_get(mean_vec, (*eff_lbl_itr).first);
					init_val += (*lbl_itr).second * newPerson->GetGenotype(tmp_loc.GetChromID(), tmp_loc.GetID());
					gsl_vector_set(mean_vec, (*eff_lbl_itr).first, init_val);
				}else if(i==0){
					std::cerr << "WARNING: Unknown locus " << (*lbl_itr).first << ", ignoring.\n";
				}

				++lbl_itr;
			}

			++eff_lbl_itr;
		}

		// Now walk through the effect matrix indexed by positions
		map<int, map<pair<int, int>, double> >::const_iterator eff_idx_itr = effect_matrix_idx.begin();
		map<pair<int, int>, double>::const_iterator idx_itr;
		ChromPool* tmp_cpool = 0;
		while(eff_idx_itr != effect_matrix_idx.end()){
			idx_itr = (*eff_idx_itr).second.begin();
			while(idx_itr != (*eff_idx_itr).second.end()){
				unsigned int chr_idx = (*idx_itr).first.first;
				unsigned int snp_idx = (*idx_itr).first.second;
				tmp_cpool = pools.GetChromosome(chr_idx);
				if(!tmp_cpool || snp_idx >= tmp_cpool->GetLociCount()){
					if(i == 0){
						std::cerr << "WARNING: Unknown locus at Chrom " << chr_idx +1 << ", Pos " << snp_idx+1 << ", ignoring.\n";
					}
				}else{
					double init_val = gsl_vector_get(mean_vec, (*eff_idx_itr).first);
					init_val += (*lbl_itr).second * newPerson->GetGenotype(chr_idx, snp_idx);
					gsl_vector_set(mean_vec, (*eff_idx_itr).first, init_val);
				}

				++idx_itr;
			}

			++eff_idx_itr;
		}

		// Now, generate a vector of N(0,1) R.V.s and put it into the norm_vec
		// Start by generating ceil(n/2) pairs of U[0,1) R.V.s and convert them
		// using the Box-Muller transform
		bool is_odd = n_effects % 2;
		for (int j=0; j<(n_effects+1)/2; j++){
			double r1 = Utility::Random::globalGenerator.drand();
			double r2 = Utility::Random::globalGenerator.drand();

			gsl_vector_set(norm_vec, 2*j, sqrt(-2*log(r1))*cos(2*M_PI*r2) );

			if(!is_odd || 2*j+1 >= n_effects){
				gsl_vector_set(norm_vec, 2*j+1, sqrt(-2*log(r1))*sin(2*M_PI*r2) );
			}
		}

		// Now, the easy part (believe it or not!)
		// A vector x = \mu + L * z should have the desired properties where:
		//   \mu is the mean vector
		//   L is the lower triangular of the cholesky decomposition
		//   z is a vector of IID N(0,1) RVs

		// First, calculate L*z and put the result into z
		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, covar_matrix, norm_vec);

		// Now, add \mu to L*z and put the result into z
		gsl_blas_daxpy(1, mean_vec, norm_vec);

		// Convert norm_vec to something a little more usable... perhaps a vector?
		effect_map[newPerson].reserve(n_effects);
		vector<double>& eff_vector = effect_map[newPerson];
		eff_vector.clear();
		for(int j=0; j<n_effects; j++){
			eff_vector.push_back(gsl_vector_get(norm_vec, j));
		}

	}


}

void ContinuousMultivarSample::ApplyPhenocopyError(PoolManager* mgr, PenetranceModel* model){
	// A no-op function because phenotypic error is nonsensical here
	// We may need to figure out what "phenotypic error" is
}

void ContinuousMultivarSample::Purge(){
	Sample::Purge();
	effect_map.clear();
}

void ContinuousMultivarSample::Reset(){
	Purge();
}

string ContinuousMultivarSample::Details(){
	return "Continuous Multivariable Effects: " + name;
}
string ContinuousMultivarSample::GetDescriptor(){
	stringstream ss;
	ss.precision(0);
	if (genoError > 0.0)
		ss<<"-"<<setprecision(2)<<(genoError*100.0)<<"GE";
	if (phenoError > 0.0)
		ss<<"-"<<setprecision(2)<<(phenoError*100.0)<<"PE";
	if (missingData > 0.0)
		ss<<"-"<<setprecision(2)<<(missingData*100.0)<<"M";
	ss<<setprecision(2)<<"-"<<name<<".mdr";
	return ss.str();
}
string ContinuousMultivarSample::GetType(){
	return "CMV";
}
string ContinuousMultivarSample::GetLabel(){
	return name;
}
string ContinuousMultivarSample::GetSummary(){
	stringstream ss;
	ss << n_indiv << " Individuals";
	return ss.str();
}

string ContinuousMultivarSample::GetDetails(){
	stringstream ss;
	ss << "DATASET " << GetType() << " " << name << " " << n_indiv << " "
			<< effect_fn << " " << covar_fn << " "
			<< genoError<<" "<<phenoError<<" "<<missingData;
	return ss.str();
}
int ContinuousMultivarSample::WriteDataset(std::ostream& os, uint* gtCounts){
	uint count = people.size();
	map<Individual*, vector<double> >::const_iterator eff_itr;
	int skipped = 0;
	for (uint i=0; i<count; i++){
		eff_itr = effect_map.find(people[i]);
		if(eff_itr == effect_map.end()){
			std::cerr << "WARNING: cannot find the effects for an individual.\n";
			++skipped;
		}else{
			vector<double>::const_iterator d_itr = (*eff_itr).second.begin();
			while(d_itr != (*eff_itr).second.end()){
				os << *d_itr << " ";
				++d_itr;
			}
			people[i]->WriteMDR(os, gtCounts, false);
		}
	}

	return count-skipped;
}
int ContinuousMultivarSample::WriteBinaryDataset(std::ostream&, std::ostream&){
	std::cerr << "WARNING: Writing binary dataset not implemented for Continuous Multivariable.\n";
	return 0;
}


void ContinuousMultivarSample::GenerateReport(std::ostream& os, uint padding){
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(genoError*100.0)<<"% Genotype Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(phenoError*100.0)<<"% Phenocopy Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(missingData*100.0)<<"% Missing Data "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<n_effects<<" Effect Variables"<<endl;
	os<<setw(padding - 5)<<"Continuous Multivariable Sample : "<<n_indiv<<" Individuals" << endl;
}

void ContinuousMultivarSample::LoadBinarySample(ifstream* genotypes, ifstream* meta, uint peopleCount, uint chromCount, vector<uint>& chrId){
	std::cerr << "WARNING: Loading binary dataset not implemented for Continuous Multivariable.\n";
}


void ContinuousMultivarSample::WriteSnpDetails(const char* project, PoolManager* pools){
	// Lifted directly from BasicSample.  This probably needs to be in a better place

	stringstream filename;
	filename<<project<<".map";
	ofstream mapFile(filename.str().c_str());
	PoolManager::Iterator itr = pools->GetIterator();
	ChromPool *pool = itr.GetNext();

	while (pool) {
		LocusArray &loci = pool->GetLoci();
		LocusArray::iterator loc = loci.begin();
		LocusArray::iterator end = loci.end();

		while (loc != end) {
			mapFile<<pool->GetLabel()<<" "<<loc->GetLabel()<<" "<<loc->GetLocation()<<"\n";
			loc++;
		}
		pool = itr.GetNext();
	}
}

ContinuousMultivarSample::~ContinuousMultivarSample() {
	gsl_matrix_free(covar_matrix);
}

} /* namespace Simulation */
