#ifndef DATAGENUTIL_H
#define DATAGENUTIL_H

#include "gmmutil.h"
#include <boost/concept_check.hpp>

namespace datagenutil
{
	
	struct BoundingBox
	{
		Vector min;
		Vector max;
	};
	
	/**
	 * generates a random Gaussian mixture model.
	 * For each generated covariance matrix this method ensures that the proportion between its maximum and minimum sqrt(eigenvalue) is in [minSqrtEWProportion,maxSqrtEWProportion]
	 * and that the minimum eigenvalue is in [minSqrtEW, maxSqrtEW]. 
	 * The entries of the means are uniformly distributed in a cube and then scaled such that the separation is equal to the given separation value.
	 */
	GMMDesc generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type separation,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp = -1);
		
	/**
	 * generates a random Gaussian mixture model.
	 * For each generated covariance matrix this method ensures that the proportion between its maximum and minimum sqrt(eigenvalue) is in [minSqrtEWProportion,maxSqrtEWProportion]
	 * and that the minimum eigenvalue is in [minSqrtEW, maxSqrtEW]. 
	 * The entries of the means are uniformly distributed in [minMean, maxMean].
	 */
	GMMDesc generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
												fp_type minMean, fp_type maxMean,
												fp_type minSqrtEW, fp_type maxSqrtEW, 
												fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
												fp_type weightExp,
												fp_type minSqrtEWExp
 											);
	
	/**
	 * computes the separation as defined by Dasgupta, i.e.
	 * min_ij { ||mu_i - mu_j|| / sqrt( max(trace(Sigma_i), trace(Sigma_j) ) ) }.
	 */
	fp_type getSeparation(GMMDesc const& gmmdesc);
	
	/**
	 * generates a dataset from according to a given Gaussian mixture model.
	 */
	void generateInputFromGMM(commonutil::DataSet&, GMMDesc const&, idx_type, std::mt19937& gen);
			
		/**
	 * adds noise points uniformly at random in an enlarged bounding box of the given data.
	 * The enlarged bounding box has sidelengths which are 1.2 times as large as the side lengths of the actual bounding box
	 * and has the same center.
	 */
	void addUniformNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen);
	
	/**
	 * draws points uniformly at random from the box that is defined by anchor+{0,1}^d*interval, where d is the dimension of the interval and center.
	 */
	void addUniformBox(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen,  Vector const& anchor, Vector const& interval);
	
	/**
	 * returns the minimum and maximum entries that points in the given dataset have.
	 */
	BoundingBox boundingBox(commonutil::DataSet const& dataset);
	
	/**
	 * multiplies the vector returned by datagenutil::getExpIncrWeights() times the given number and rounds the resulting entries.
	 * It is ensured that all vector entries sum up to the number by setting the last entry of the vector to  number - sum of the remaining entries.
	 */
	Vector getExpIncrNumbers(idx_type k, fp_type exponent, idx_type number);
	
	/**
	 * returnes a vector (2^0, 2^(exponent), 2^(2*exponent),...,2^(k*exponent)) normalized by the sum of its entries.
	 * I.e., the entries of the returned vector are exponentially increasing, are in [0,1], and sum up to 1.
	 */
	Vector getExpIncrWeights(idx_type k, fp_type exponent);
};

#endif // GENUTIL_H
