#ifndef GMMUTIL_H
#define GMMUTIL_H

#include "base.h"
#include "gmmdesc.h"
#include "commonutil.h"

namespace gmmutil
{

	/**
	 * returns the minimum squared (pairwise) distance between the given points.
	 */
	const fp_type minSquaredDistance(Matrix const&);

	/**
	 * returns the density of each of the points as defined by the gaussian distribution (which is defined by the given mean and covariance).
	 */
	Vector gaussianDensity(Matrix const& points, Vector const& mean, Matrix covariance);
	
	/**
	 * returns the densitiy of each of the points as defined by the given mixture model.
	 */
	Vector gmmDensity(Matrix const& points, GMMDesc const& desc);
	
	/**
	 * returns the negative-log-likelihood of each of the points w.r.t. a gaussian distribution.
	 */
	Vector gaussianNLL(Matrix const&, Vector const&, Matrix);
	
	/**
	 * returns the negative-log-likelihood of each of the points w.r.t. the gaussian mixture model.
	 */
	Vector gmmNLL(Matrix const&, GMMDesc const&);
	
	
	fp_type nll(commonutil::DataSet const&, GMMDesc const&);
	
	/**
	 * returns the mahalanobis distance between the given gaussian and each of the points.
	 */
	Vector gaussianMahalanobis(Matrix const& points, Vector const& mean, Matrix covariance);
		
	/**
	 * returns the minimal Mahalanobis distance to a mean w.r.t. the respective covariance matrix
	 * for each of the given points.
	 */
	Vector gmmMinMahalanobis(Matrix const&, GMMDesc const&);
	
	/**
	 * returns densities defined by factor * gmmMinMahalanobis + (1-factor) * uniform.
	 */
	Vector adaptiveDensities(commonutil::DataSet const&, GMMDesc const&, fp_type const ratio);


	/**
	 * returns the optimal (k=1)-solution. 
	 * Guarantees that the returned covariance matrix is positive definite: In case the computed optimal covariance is not positive definite it
	 * returns an approximated spherical covariance or the identity matrix.
	 */
	GMMDesc optimalGaussian(commonutil::DataSet const& input);

	
	/**
	 * returnes a matrix whose (i,j)-th entry contains w(x_j) * w_i * N(x_j | mu_i, Sigma_i).
	 * 
	 */
	Matrix pmatrix(commonutil::DataSet const& input, GMMDesc gmmdesc);
		
	
	/**
	 * returns the assignment of each point to its likeliest component, i.e. argmax_j=1..k{ w_j N(x | \mu_j,\sigma_j) }.
	 */
	std::vector<idx_type> gaussPartition(Matrix const& points, GMMDesc const& gmmdesc);
	

	idx_type minNLL(commonutil::DataSet const&, std::vector<GMMDesc> const&);
}


#endif
