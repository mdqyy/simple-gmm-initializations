#ifndef GAUSSIANMIXTURE_H
#define GAUSSIANMIXTURE_H

#include "base.h"
#include "gmmdesc.h"
#include "gaussian.h"

#include <random>
#include <boost/concept_check.hpp>

/**
 * @brief Data drawn from a gaussian mixture model.
 */
struct CompleteData
{
	Vector point;
	idx_type source;
};


/**
 * @brief Gaussian mixture model distribution.
 */
class GaussianMixture
{
public:

	GaussianMixture(GMMDesc const&);

	/**
	 * evaluates the density of the GMM distribution at the given point x.
     */
	virtual fp_type density(Vector const& x) const;

	/**
	 * computes the negative log-likelihood of the density at the given point x.
     */
	virtual fp_type nll(Vector const& x) const;
	
	virtual fp_type minNLL(Vector const& x) const;
	
	virtual fp_type minSquaredMahalanobis(Vector const& x) const;
	
	template<typename RndEngine> Vector draw(RndEngine&) const;
	
	template<typename RndEngine> CompleteData drawCompleteData(RndEngine& re) const;

private:
	GMMDesc gmmdesc;
	std::vector<SingleGaussian> gaussians;
};

template<typename RndEngine>
Vector GaussianMixture::draw(RndEngine& re) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	if (k==0)
		throw "There are no gaussians in this mixture!";

	std::uniform_real_distribution<> urd(0, 1);
	fp_type r = urd(re);
	
	std::size_t i = 0;
	fp_type w = this->gmmdesc.weights[0];
	while (r>w&&i<k)
	{
		r -= w;
		w = this->gmmdesc.weights[++i];
	}
				
	return this->gaussians[i].draw(re);
}

template<typename RndEngine>
CompleteData GaussianMixture::drawCompleteData(RndEngine& re) const
{
	std::size_t k = this->gmmdesc.weights.size();
	
	assert(this->gaussians.size()==k);	
	
	if (k==0)
		throw "There are no gaussians in this mixture!";

	std::uniform_real_distribution<> urd(0, 1);
	fp_type r = urd(re);
	
	std::size_t i = 0;
	fp_type w = this->gmmdesc.weights[0];
	while (r>w&&i<k)
	{
		r -= w;
		w = this->gmmdesc.weights[++i];
	}
				
	CompleteData cd;
	cd.point = this->gaussians[i].draw(re);
	cd.source = i;
	return cd;
}

#endif 
