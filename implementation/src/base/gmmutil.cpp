#include "gmmutil.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <ctime>
#include <cmath>

#include "linalgutil.h"
#include "gaussian.h"
#include "kmeansutil.h"
#include "commonutil.h"

#include "../settings/settings.h"


fp_type gmmutil::nll(commonutil::DataSet const& data, GMMDesc const& desc)
{
	assert(desc.means.cols()==desc.weights.size() && desc.covariances.size()==desc.weights.size());
	
	Vector v = gmmNLL(data.points, desc);
	idx_type n = data.points.cols();
	
	assert(v.size()==n && data.weights.size()==n);
	
	for (idx_type i=0; i<n; ++i)
		v[i] *= data.weights[i];
	
	return v.sum();
}

const fp_type gmmutil::minSquaredDistance(Matrix const& points)
{
	fp_type min = 0;
	idx_type n = points.cols();
	for (idx_type i=1; i<n; ++i)
	{
		Vector sn = (points.rightCols(n-i).colwise()-points.col(i-1)).colwise().squaredNorm();
		idx_type m = sn.size();
		for (idx_type j=0; j<m; ++j)
			if (sn[j]>0&&(min==0||sn[j]<min))
				min = sn[j];
	}
	return min;
}

Vector gmmutil::gaussianDensity(Matrix const& points, Vector const& mean, Matrix covariance)
{
	linalg::cholesky(covariance);
	Vector v = covariance.diagonal();
	idx_type d = v.size();
	fp_type logSqrt = d*log(2*M_PI)/2;
	for (idx_type i=0; i<d; ++i)
		logSqrt += log(v(i));
	linalg::ltrinv(covariance);
	v.noalias() = (covariance*(points.colwise()-mean)).colwise().squaredNorm();
	
	idx_type n = v.size();
	for (idx_type i=0; i<n; ++i)
		v(i) = exp(-0.5*v(i)-logSqrt);
	
	return v;
}

Vector gmmutil::gmmDensity(Matrix const& points, GMMDesc const& desc)
{
	std::size_t k = desc.weights.size();
	
	assert(desc.means.cols()==k && desc.covariances.size()==k);
	
	Vector sum = Vector::Zero(points.cols());
	for (std::size_t i=0; i<k; ++i)
		sum.noalias() += desc.weights[i]*gaussianDensity(points, desc.means.col(i), desc.covariances.at(i));
	
	return sum;
}


Vector gmmutil::gaussianNLL(Matrix const& points, Vector const& mean, Matrix covariance)
{
	linalg::cholesky(covariance);
	Vector v = covariance.diagonal();
	idx_type d = v.size();
	fp_type logSqrt = d*log(2*M_PI)/2;
	for (idx_type i=0; i<d; ++i)
		logSqrt += log(v(i));
	linalg::ltrinv(covariance);
	v.noalias() = (covariance*(points.colwise()-mean)).colwise().squaredNorm();
	
	d = v.size();
	for (idx_type i=0; i<d; ++i)
		v(i) = 0.5*v(i)+logSqrt;
	
	return v;
}


Vector gmmutil::gmmNLL(Matrix const& points, GMMDesc const& desc)
{
	idx_type k = desc.weights.size();
	idx_type n = points.cols();
	if (k==0)
		return Vector::Zero(n);
	
	Vector densities = gmmDensity(points, desc);
	for (size_t i=0; i<n; ++i)
		if (densities[i]>0 && densities[i]<std::numeric_limits<fp_type>::infinity())
			densities[i] = -log(densities[i]);
		else
		{
			Vector nlls(k);
			for (idx_type j=0; j<k; ++j)
				nlls[j] = gaussianNLL(points.col(i), desc.means.col(j), desc.covariances.at(j))[0]
				-log(desc.weights[j]);
			idx_type index;
			if (densities[i]<=0)
			{
				densities[i] = nlls.minCoeff(&index);
				if(commonSettings().verbose)
					std::cout << "gmmutil::gmmNLL() - cost of point in gmm was zero!!! Approximated by cost in component " << index+1 << std::endl;
			}
			else
			{
				densities[i] = nlls.minCoeff(&index);
				if(commonSettings().verbose)
					std::cout << "gmmutil::gmmNLL() - cost of point in gmm was infinity!!! Approximated by cost in component " << index+1 << std::endl;
			}
		}
		
	return densities;
}


Vector gmmutil::gaussianMahalanobis(Matrix const& points, Vector const& mean, Matrix covariance)
{
	linalg::cholesky(covariance);
	linalg::ltrinv(covariance);
	return (covariance*(points.colwise()-mean)).colwise().squaredNorm();
}


Vector gmmutil::gmmMinMahalanobis(Matrix const& points, GMMDesc const& desc)
{
	idx_type n = points.cols();
	idx_type k = desc.weights.size();
	
	assert(k>0 && n>0);
	
	Matrix mahalanobis(k,n);
	
	for (idx_type i=0; i<k; ++i)
		mahalanobis.row(i) = gmmutil::gaussianMahalanobis(points, desc.means.col(i),desc.covariances.at(i));
	
	return mahalanobis.colwise().minCoeff();
}

Vector gmmutil::adaptiveDensities(commonutil::DataSet const& data, GMMDesc const& desc, fp_type const ratio)
{
	assert (ratio>=0 && ratio<=1);
	
	Vector densities = gmmutil::gmmMinMahalanobis(data.points, desc);
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
		densities[j] = ratio*densities[j]+(1.-ratio)*iidDens;
		densities[j] *= data.weights[j];
	}
	
	sum = densities.sum();
	if (sum>0)
		densities /= sum;
	
	return densities;
}

std::vector<idx_type> gmmutil::gaussPartition(Matrix const& points, GMMDesc const& gmmdesc)
{
	const idx_type n = points.cols();
	const idx_type k = gmmdesc.means.cols();
	std::vector<idx_type> partition;
	// initial solution: component 0
	Vector highestDensities = gmmutil::gaussianDensity(points, gmmdesc.means.col(0), gmmdesc.covariances[0]);
	for(idx_type i=0; i<n; ++i)
		partition.push_back((idx_type)0);
	Vector tmpDensities = Vector(n);
	// iterate over components 1 to k-1
	for(idx_type j=1; j<k; ++j)
	{
		tmpDensities = gmmutil::gaussianDensity(points, gmmdesc.means.col(j), gmmdesc.covariances[j]);
		for(idx_type i=0; i<n; ++i)
		{
			if(highestDensities(i) < tmpDensities(i))
			{
				partition.at(i) = j;
				highestDensities(i) = tmpDensities(i);
			}
		}
	}
	return partition;
}

Matrix gmmutil::pmatrix(commonutil::DataSet const& input, GMMDesc gmmdesc)
{
	const idx_type k = gmmdesc.means.cols();
	const idx_type n = input.points.cols();
	Matrix pmatrix(k,n);
	for (idx_type i=0; i<k; ++i)
		pmatrix.row(i) = gmmdesc.weights[i]*gmmutil::gaussianDensity(input.points, gmmdesc.means.col(i), gmmdesc.covariances[i]);
	for (idx_type j=0; j<n; ++j)
		pmatrix.col(j) *= input.weights[j];
	return pmatrix;
}

GMMDesc gmmutil::optimalGaussian(commonutil::DataSet const& input)
{
	const idx_type n = input.points.cols();
	const idx_type d = input.points.rows();
	
	const fp_type totalWeight = input.weights.sum();
	
	GMMDesc output;
	
	output.means = Matrix::Zero(d,1);	
	for(idx_type i=0; i<n; ++i)
		output.means.col(0) += input.weights(i)*input.points.col(i);
	output.means.col(0) = output.means.col(0) / totalWeight;
	
	Matrix covar = Matrix::Zero(d,d);
	Matrix diffs = (input.points.colwise()-output.means.col(0));
	fp_type spherical = 0;
	for(idx_type i=0; i<n; ++i)
	{
		covar.noalias() += (input.weights(i)*diffs.col(i))*diffs.col(i).transpose();
		spherical += input.weights[i]*diffs.col(i).squaredNorm();
	}
	covar = covar/totalWeight;
	
	if (!linalg::spd(covar))
	{
		fp_type factor = spherical/(totalWeight*d);
		if(std::isfinite(factor) && factor>0)
			covar = Matrix::Identity(d,d)*(spherical/(totalWeight*d));
		else
			covar = Matrix::Identity(d,d);
		
	}
	
	output.covariances.push_back(covar);
	output.weights = Vector(1);
	output.weights(0) = 1.;
	return output;
}


idx_type gmmutil::minNLL(commonutil::DataSet const& data, std::vector<GMMDesc> const& models)
{
	const idx_type num = models.size();
	Vector nlls(num);
	for (idx_type i=0; i<num; ++i)
		nlls[i] = nll(data, models[i]);
	idx_type index;
	nlls.minCoeff(&index);
	return index;
}



