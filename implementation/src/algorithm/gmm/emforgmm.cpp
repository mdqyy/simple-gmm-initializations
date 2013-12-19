#include "emforgmm.h"

#include <algorithm>
#include <iostream>
#include <cmath>

#include "../../base/linalgutil.h"
#include "../../base/initutil.h"
#include "../../base/gmmutil.h"

#include "../../settings/settings.h"


EMforGMM::EMforGMM(commonutil::DataSet const& ds, bool v, uint32_t s)
	: RandomAlgorithm(ds,v,s)
{
}

EMforGMM::EMforGMM(commonutil::DataSet const& ds, bool v, std::mt19937& gen)
	: RandomAlgorithm(ds,v,0)
{
	this->gen = gen;
}

void EMforGMM::init(unsigned int k)
{
	std::mt19937 gen(this->seed++);
	this->desc = initutil::adaptiveSphericalGMM(*this->input, k, gen);
	this->change = true;
}

void EMforGMM::init(GMMDesc const& desc)
{
	this->desc = desc;
	this->change = true;
}

void EMforGMM::run(unsigned int numSteps)
{
	clock_t start, end;
	start = clock();

	std::size_t k = this->desc.weights.size();

	assert(this->desc.means.cols()==k && this->desc.covariances.size()==k);

	if (k==0)
		gmmlab_throw("No or empty initial solution given.");

	idx_type n = this->input->points.cols();
	idx_type d = this->input->points.rows();

	if (n==0 || d==0)
		gmmlab_throw("Input is empty.");

	fp_type totalWeight = this->input->weights.sum();
	if (totalWeight<=0)
		gmmlab_throw("Total input weight less or equal to zero.");

	fp_type minWeight = totalWeight;
	for (idx_type i=0; i<n; ++i)
		if (this->input->weights[i]>0 && this->input->weights[i]<minWeight)
			minWeight = this->input->weights[i];

	Matrix pmatrix(k,n);
	for (unsigned int step=0; step<numSteps; ++step)
	{
		for (idx_type i=0; i<k; ++i)
			pmatrix.row(i) = this->desc.weights[i]*gmmutil::gaussianDensity(this->input->points,
				this->desc.means.col(i), this->desc.covariances[i]);

		for (idx_type j=0; j<n; ++j)
		{
			fp_type sum = pmatrix.col(j).sum();

			if (sum<=0)
			{
				if (this->verbose)
					std::cout << "EMforGMM::run() - no gaussian is responsible for point " << j << std::endl;
				pmatrix.col(j) = this->desc.weights;
			}
			else
				for (idx_type i=0; i<k; ++i)
					pmatrix(i,j) /= sum;
		}

		for (idx_type j=0; j<n; ++j)
			pmatrix.col(j) *= this->input->weights[j];

		for (idx_type i=0; i<k; ++i)
		{
			fp_type resp = pmatrix.row(i).sum();
			if (resp>minWeight)	// only update mean and covariance, if the gaussian has
								// non-negligible responsibility
			{
				this->desc.weights[i] = resp/totalWeight; // may become zero, no problem

				Vector p = Vector::Zero(d);
				for (idx_type j=0; j<n; ++j)
					p.noalias() += pmatrix(i,j)*this->input->points.col(j);
				this->desc.means.col(i) = p/resp;
				
				Matrix m = Matrix::Zero(d,d);
				for (idx_type j=0; j<n; ++j)
				{
					Vector y = this->input->points.col(j) - this->desc.means.col(i);
					m.noalias() += pmatrix(i,j)*(y*y.transpose());
				}
				m = m/resp;
				if (linalg::spd(m))
					this->desc.covariances[i] = m;
				else
				{
					m = (this->desc.covariances[i]+m)/2.0;
					if (linalg::spd(m))
					{
						this->desc.covariances[i] = m;
						if (this->verbose)
							std::cout << "EMforGMM::run() - matrix not spd, mixing with old covariance of gaussian " << i+1 << std::endl;
					}
					else
						if (this->verbose)
							std::cout << "EMforGMM::run() - matrix not spd, keeping old covariance " << i+1 << std::endl;
				}
			}
			else
			{
				this->desc.weights[i] = 1.0/k;

				this->desc.means.col(i) = this->input->points.col(commonutil::randomIndex(this->input->weights, this->gen));
				if (this->verbose)
					std::cout << "EMforGMM::run() - replaced empty mean " << i+1 << " by random input point" << std::endl;

				Vector sqDist = (this->desc.means.colwise()-this->desc.means.col(i)).colwise().squaredNorm();
				fp_type maxd = sqDist.maxCoeff();
				for (std::size_t j=0; j<k; ++j)
					if (this->desc.means.col(i)==this->desc.means.col(j))
						sqDist[j] = maxd;

				fp_type factor = (1.0/(2*d))*sqDist.minCoeff();
				if (!std::isfinite(factor) || factor <= 0) // i.e., factor*I_d would not be positive definite
				{
					this->desc.covariances[i].noalias() = (1.0/(2*d))*Matrix::Identity(d,d);
					if (this->verbose)
						std::cout << "EMforGMM::run() - replaced covariance " << i+1 << " by scaled identity matrix, factor = 1.0/(2*d)" << std::endl;
				}
				else
				{
					this->desc.covariances[i].noalias() = factor*Matrix::Identity(d,d);
					if (this->verbose)
						std::cout << "EMforGMM::run() - replaced covariance " << i+1 << " by scaled identity matrix, factor = " << sqDist.minCoeff() << "/(2*d)" << std::endl;
				}
				
			}
		}

		// finally normalize the weights (maybe weights were replaced and the sum is not 1)
		this->desc.weights /= this->desc.weights.sum();
	}

	this->change = true;

	end = clock();
	this->runtime += double(end-start)/CLOCKS_PER_SEC;
}

EMforGMM* EMforGMM::toEMforGMM(GMMAlgorithm* a)
{
	return dynamic_cast<EMforGMM*>(a);
}

