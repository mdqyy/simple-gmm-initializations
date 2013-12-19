#ifndef KMEANSUTIL_H
#define KMEANSUTIL_H

#include "base.h"
#include "gmmdesc.h"
#include "gmmutil.h"

#include <iostream>

namespace kmeansutil
{

	/**
	 * creates a GMM based on a set of k means and the clustering induced by these means 
	 * (i.e. each point is assigned to its nearest mean) as follows:
	 * - The mean of a component is set to the centroid of the corresponding  cluster.
	 * - The covariance of a component is set to the optimal covariance of the (k=1)-problem
	 *   that is given by the corresponding cluster and its mean.
	 * - The weights are estimated by the fraction of points assigned to the corresponding
	 *   mean.
	 */
	GMMDesc meansToGMM(commonutil::DataSet const& input, Matrix const& means, bool verbose = false);
	
	/**
	 * wrap some means into a GMMDesc. That is, the means are used as centers of the components,
	 * while the weights and covariances are estimated by considering the kmeans clusters given by the means.
	 */
	GMMDesc wrapMeans(commonutil::DataSet const& input, Matrix const& means, bool computeWeightAndCovar);
  
	/**
	 * computes the sum of the squared distances between a point and its closest mean (weighted by the weight of the point).
	 */
	fp_type kmeanscost(commonutil::DataSet const& data, Matrix const& means);
		
	Vector squaredDistances(Matrix const& points, Vector const& point);
	fp_type minSquaredDistance(Matrix const& points, Vector const& point);
	Vector minSquaredDistances(Matrix const& points, Matrix const& means);
	
	/**
	 * returns the assignment of each point to its nearest mean.
	 */
	std::vector<idx_type> kmeansPartition(Matrix const& points, Matrix const& means);
	
	/**
	 * returns the index of the mean that has the smalles squared distance to the point.
	 */
	idx_type nearest(Vector const& point, Matrix const& means);
	
	/**
	 * the density of a point is given by its squared distance to its nearest center (normalized by 
	 * the sum of all squared distances * 2) plus ratio/2 times uniform distribution.
	 */
	Vector adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio);
}


#endif
