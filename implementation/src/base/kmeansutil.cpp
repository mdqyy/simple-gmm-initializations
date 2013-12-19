#include "kmeansutil.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <ctime>

#include "gmmutil.h"
#include "linalgutil.h"

#include "../settings/settings.h"

GMMDesc kmeansutil::meansToGMM(commonutil::DataSet const& input, Matrix const& means, bool verbose)
{
	GMMDesc desc;

	std::size_t k = means.cols();

	if (k==0)
		gmmlab_throw("initutil::kmeansToGMM() - No means given.");

	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("initutil::kmeansToGMM() - Input is empty.");

	fp_type totalWeight = input.weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("initutil::kmeansToGMM() - Total input weight less or equal to zero.");

	Matrix pmatrix(k,n);
	for (idx_type i=0; i<k; ++i)
		pmatrix.row(i) = (input.points.colwise()-means.col(i)).colwise().squaredNorm();

	std::vector<idx_type> indices(n);
	for (idx_type i=0; i<n; ++i)
		pmatrix.col(i).minCoeff(&indices[i]);

	desc.weights = Vector::Zero(k);
	desc.means = Matrix::Zero(d,k);
	desc.covariances = std::vector<Matrix>(k, Matrix::Zero(d,d));
	std::vector<fp_type> sphericals(k, 0);

	for (std::size_t i=0; i<n; ++i)
	{
		desc.weights[indices[i]] += input.weights[i];
		desc.means.col(indices[i]).noalias() += input.weights[i]*input.points.col(i);
	}

	for (std::size_t i=0; i<k; ++i)
		if (desc.weights[i] > 0)
			desc.means.col(i) = desc.means.col(i)/desc.weights[i];
		else
			desc.means.col(i) = means.col(i);

	for (std::size_t i=0; i<n; ++i)
	{
		Vector y = input.points.col(i) - desc.means.col(indices[i]);
		desc.covariances[indices[i]].noalias() += input.weights[i]*(y*y.transpose());
		sphericals[indices[i]] += input.weights[i]*y.squaredNorm();
	}

	for (std::size_t i=0; i<k; ++i)
	{
		if (desc.weights[i] > 0)
		{
			desc.covariances[i] = desc.covariances[i]/desc.weights[i];
			if (!linalg::spd(desc.covariances[i]))
			{
				desc.covariances[i] = Matrix::Identity(d,d)*(sphericals[i]/(desc.weights[i]*d));
				if (!linalg::spd(desc.covariances[i]))
				{
					desc.covariances[i] = Matrix::Identity(d,d);
					if (verbose)
						std::cout << "initutil::kmeansToGMM() - replaced non-spd covariance " << i+1 << " with unit sphere covariance." << std::endl;
				}
				else
					if (verbose)
						std::cout << "initutil::kmeansToGMM() - replaced non-spd covariance " << i+1 << " with spherical covariance." << std::endl;
			}
		}
		else
		{
			desc.weights[i] = totalWeight/k;
			desc.covariances[i] = Matrix::Identity(d,d);
			if (verbose)
				std::cout << "initutil::kmeansToGMM() - replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
		}

	}

	// finally normalize the weights
	desc.weights /= desc.weights.sum();
	
	return desc;
}


GMMDesc kmeansutil::wrapMeans(commonutil::DataSet const& input, Matrix const& means, bool computeWeightAndCovar)
{
	GMMDesc desc;

	std::size_t k = means.cols();

	if (k==0)
		gmmlab_throw("kmeansutil::wrapMeans() - No means given.");

	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("kmeansutil::wrapMeans() - Input is empty.");

	fp_type totalWeight = input.weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("kmeansutil::wrapMeans() - Total input weight less or equal to zero.");

	desc.means = means;
	
	if(computeWeightAndCovar){
	
	  // assign each point to its nearest center
	  Matrix pmatrix(k,n);
	  for (idx_type i=0; i<k; ++i)
		  pmatrix.row(i) = (input.points.colwise()-means.col(i)).colwise().squaredNorm();
	  std::vector<idx_type> indices(n);
	  for (idx_type i=0; i<n; ++i)
		  pmatrix.col(i).minCoeff(&indices[i]);
	  
	  // estimate weights and covariances
	  desc.weights = Vector::Zero(k);
	  desc.covariances = std::vector<Matrix>(k, Matrix::Zero(d,d));
	  std::vector<fp_type> sphericals(k, 0);
	  for (std::size_t i=0; i<n; ++i)
	  {
		  desc.weights[indices[i]] += input.weights[i];
		  Vector y = input.points.col(i) - desc.means.col(indices[i]);
		  desc.covariances[indices[i]].noalias() += input.weights[i]*(y*y.transpose());
		  sphericals[indices[i]] += input.weights[i]*y.squaredNorm();
	  }
	  for (std::size_t i=0; i<k; ++i)
	  {
		  if (desc.weights[i] > 0)
		  {
			  desc.covariances[i] = desc.covariances[i]/desc.weights[i];
			  if (!linalg::spd(desc.covariances[i]))
			  {
				  desc.covariances[i] = Matrix::Identity(d,d)*(sphericals[i]/(desc.weights[i]*d));
				  if (sphericals[i]/(desc.weights[i]*d) <= 0)
				  {
					  desc.covariances[i] = Matrix::Identity(d,d);
					  if (commonSettings().verbose)
						  std::cout << "kmeansutil::wrapMeans() - replaced non-spd covariance " << i+1 << " with unit sphere covariance." << std::endl;
				  }
				  else
					  if (commonSettings().verbose)
						  std::cout << "kmeansutil::wrapMeans() - replaced non-spd covariance " << i+1 << " with spherical covariance." << std::endl;
			  }
		  }
		  else
		  {
			  desc.weights[i] = totalWeight/k;
			  desc.covariances[i] = Matrix::Identity(d,d);
			  if (commonSettings().verbose)
				  std::cout << "kmeansutil::kmeansToGMM() - replaced empty covariance " << i+1 << " with unit sphere covariance." << std::endl;
		  }

	  }

	  // finally normalize the weights
	  desc.weights /= desc.weights.sum();
	
	}
	
	return desc;
}

fp_type kmeansutil::kmeanscost(commonutil::DataSet const& data, Matrix const& means)
{
	Vector v = kmeansutil::minSquaredDistances(data.points, means);
	idx_type n = data.points.cols();

	assert(v.size() == n && data.weights.size() == n);

	for (idx_type i = 0; i < n; ++i)
		v[i] *= data.weights[i];
	
	return v.sum();
}


Vector kmeansutil::squaredDistances(Matrix const& points, Vector const& point)
{
	return (points.colwise() - point).colwise().squaredNorm();
}

				
fp_type kmeansutil::minSquaredDistance(Matrix const& points, Vector const& point)
{
	Vector dists = kmeansutil::squaredDistances(points, point);
	idx_type index;
	dists.minCoeff(&index);
	return dists[index];
}

Vector kmeansutil::minSquaredDistances(Matrix const& points, Matrix const& means)
{
	const idx_type n = points.cols();
	Vector dists = Vector::Zero(n);
	for (idx_type i = 0; i < n; ++i)
		dists[i] = minSquaredDistance(means, points.col(i));
	return dists;
}

std::vector<idx_type> kmeansutil::kmeansPartition(Matrix const& points, Matrix const& means)
{
	const idx_type n = points.cols();
	std::vector<idx_type> partition(n);
	for(idx_type i=0; i<n; ++i)
		partition.at(i) = (kmeansutil::nearest(points.col(i), means));
	return partition;
}


idx_type kmeansutil::nearest(Vector const& point, Matrix const& means)
{
	idx_type index;
	kmeansutil::squaredDistances(means, point).minCoeff(&index);
	return index;
}


Vector kmeansutil::adaptiveDensities(commonutil::DataSet const& data, Matrix const& means, fp_type const ratio)
{
	assert (ratio>=0 && ratio<=1);

	Vector densities = kmeansutil::minSquaredDistances(data.points, means);
	idx_type n = data.points.cols();

	for (idx_type i=0; i<n; ++i)
		if (!(densities[i]>0)) // test if <0 or NaN
			densities[i] = 0;

	fp_type sum = densities.sum();
	if (sum>0)
		densities /= sum;

	fp_type iidDens = fp_type(1)/n;
	for (idx_type j=0; j<n; ++j)
	{
		densities[j] = ratio*densities[j]+(1-ratio)*iidDens;
		densities[j] *= data.weights[j];
	}

	sum = densities.sum();
	if (sum>0)
		densities /= sum;

	return densities;
}


