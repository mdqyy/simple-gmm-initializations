#include "datagenutil.h"

#include "base.h"
#include "gmmdesc.h"
#include "gmmutil.h"
#include "gaussian.h"
#include "gaussianmixture.h"
#include "linalgutil.h"
#include "commonutil.h"
#include "../settings/settings.h"
#include <boost/iterator/iterator_concepts.hpp>

GMMDesc datagenutil::generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type separation,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp)
{
	
	if(minSqrtEWExp < 0)
	  gmmlab_throw("datagenutil::generateRandomGMMWithUniformMeans() - Use minSqrtEWExp > 0.");
  
	GMMDesc truth = datagenutil::generateRandomGMMWithUniformMeans(d,k,gen,0,maxSqrtEW*maxSqrtEWProportion,minSqrtEW,maxSqrtEW,minSqrtEWProportion,maxSqrtEWProportion,weightExp,minSqrtEWExp);
	
	fp_type separation_tmp  = datagenutil::getSeparation(truth);
	truth.means = separation/separation_tmp*truth.means;
	
	if (commonSettings().verbose)
		std::cout << "generated " << truth << std::endl;
		
	return truth;
}

GMMDesc datagenutil::generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type minMean, fp_type maxMean,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp)
{
	
	if(minSqrtEWExp < 0)
	  gmmlab_throw("datagenutil::generateRandomGMMWithUniformMeans() - Use minSqrtEWExp > 0.");
  
	GMMDesc truth;
	
	std::uniform_real_distribution<> urd(0, 1);	

	// (1) create means uniformly at random in a given area
	Vector mean = Vector(d);
	truth.means = Matrix::Zero(d,k);
	for(idx_type l=0; l<k; ++l)
	{	
		commonutil::fill(mean, gen, minMean, maxMean);
		truth.means.col(l) = mean;
	}
		
	// (2) create random covariances with restricted eigenvalues
		
	// create factors that determine the growth of the minimum eigenvalues
	Vector minSqrtEWExpVector = datagenutil::getExpIncrWeights(k, minSqrtEWExp);
	minSqrtEWExpVector = minSqrtEWExpVector * 1./minSqrtEWExpVector.maxCoeff();
	// shuffle these factors (unless all factors = 1/k) to avoid a 
	// correlation between minimum eigenvalues and weights (which are determined
	// deterministically
	if(minSqrtEWExp != 0)
	{
		std::vector<fp_type> shuffleVector(minSqrtEWExpVector.size());
		for(idx_type i=0; i<shuffleVector.size(); ++i)
			shuffleVector.at(i) = minSqrtEWExpVector(i);
		std::shuffle<>(shuffleVector.begin(), shuffleVector.end(), gen);
		for(idx_type i=0; i<shuffleVector.size(); ++i)
			minSqrtEWExpVector(i) = shuffleVector.at(i);
	}
	// create covariance matrizes
	Vector eigenvalues = Vector(d);
	Matrix randMatrix = Matrix(d,d);
	Matrix randOrthonormalMatrix = Matrix(d,d);
	for(idx_type l=0; l<k; ++l)
	{			
		// create random matrix M
		for (idx_type i=0; i<d; ++i)
			for (idx_type j=0; j<d; ++j)
				randMatrix(i,j)=urd(gen);
			
		// obtain random orthonormal matrix Q by the QR-decomposition of the random matrix M
		Eigen::HouseholderQR<Matrix> qr(randMatrix);
		randOrthonormalMatrix = qr.householderQ();
		
		// create random eigenvalues
		fp_type minEW, maxEW;
		minEW = minSqrtEW + minSqrtEWExpVector(l) * (maxSqrtEW - minSqrtEW);	
		maxEW = (minSqrtEWProportion + urd(gen) * (maxSqrtEWProportion-minSqrtEWProportion))*minEW;
		minEW = pow(minEW, 2);
		maxEW = pow(maxEW, 2);
		commonutil::fill(eigenvalues, gen, minEW, maxEW);
		// make sure that minEW and maxEW appear in the EW's
		eigenvalues(0) = minEW;
		eigenvalues(d-1) = maxEW;
		// Note that it is not necessary to shuffle the eigenvalues, since we multiply the resulting
		// diagonal matrix containing these eigenvalues with a completely random matrix (Q).
		
		// random covariance = Q^T * Diag(eigenvalues) * Q
		truth.covariances.push_back(randOrthonormalMatrix.transpose()*eigenvalues.asDiagonal()*randOrthonormalMatrix);
	}
	
	// (3) weights
	truth.weights = datagenutil::getExpIncrWeights(k, weightExp);
	
	return truth;
}

fp_type datagenutil::getSeparation(GMMDesc const& gmmdesc)
{
	idx_type k = gmmdesc.means.cols();
	idx_type sqrtd = sqrt(gmmdesc.means.rows());
	
	//Vector maxSqrtEigenvalues = Vector::Zero(k);
	Vector traces = Vector::Zero(k);
	fp_type separation;
	
	for(idx_type i=0; i<k; ++i)
	{
		//Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(gmmdesc.covariances.at(i));
		//if (eigensolver.info() != Eigen::Success)
		//	gmmlab_throw("datagetutil::getSeparation() - Eigensolver failed.");
	    //maxSqrtEigenvalues(i) = std::sqrt(eigensolver.eigenvalues().maxCoeff());
		traces(i) = gmmdesc.covariances.at(i).trace();
		
		for(idx_type j=0; j<i; ++j)
		{
			//fp_type tmp = (gmmdesc.means.col(i)-gmmdesc.means.col(j)).norm() / (sqrtd * std::max(maxSqrtEigenvalues(i),maxSqrtEigenvalues(j)));
			fp_type tmp = (gmmdesc.means.col(i)-gmmdesc.means.col(j)).norm() / std::sqrt(std::max(traces(i),traces(j)));
			if( (i==1&&j==0) || tmp < separation)
				separation = tmp;
		}
	}
	return separation;	
}

void datagenutil::generateInputFromGMM(commonutil::DataSet& target, GMMDesc const& desc, idx_type n, std::mt19937& gen)
{
	assert(desc.means.cols()==desc.weights.size() && desc.covariances.size()==desc.weights.size());

	const idx_type d = desc.means.rows();
	if (desc.means.size()==0)
		gmmlab_throw("datagenutil::generateNoisyInputFromGMM() - Cannot draw data from an empty gmmdesc.");

	// resize the input matrix so that each data point will form one column
	target.points.resize(d, n);
	target.weights.resize(n);

	GaussianMixture gm(desc);
	for (std::size_t i=0; i<n; ++i)
	{
		target.weights[i] = 1;
		target.points.col(i) = gm.draw(gen);
		fp_type density = gm.density(target.points.col(i));
		if (density<=0)
			std::cout << "datagenutil::generateInputFromGMM - WARNING!!! Drew point " << target.points.col(i) << " having density " << density << "." << std::endl;
	}
}



void datagenutil::addUniformNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen)
{	
	const idx_type old_n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	
	if(old_n==0||d==0)
		gmmlab_throw("datagenutil::addUniformNoise() - Given dataset must not be empty.");
	
	const BoundingBox boundingBox = datagenutil::boundingBox(dataset);
	Vector anchor = boundingBox.min;
	Vector interval = boundingBox.max-boundingBox.min;
	anchor = anchor - 0.1*interval;
	interval = interval + 0.2*interval;
	
	datagenutil::addUniformBox(dataset, d, n, gen, anchor, interval);
}



void datagenutil::addUniformBox(commonutil::DataSet& dataset, idx_type d,  idx_type n,  std::mt19937& gen, Vector const& anchor, Vector const& interval)
{		
	const idx_type old_n = dataset.points.cols();
	
	if(interval.rows()!=d || anchor.rows()!=d || interval.rows()!=d)
		gmmlab_throw("datagenutil::addPointsUniformlyAtRandom() - Interval or anchor have the wrong dimension.");
	
	// resize
	dataset.points.conservativeResize(d,old_n+n);
	dataset.weights.conservativeResize(old_n+n);
	
	// add points that are drawn uniformly at random
	std::uniform_real_distribution<> urd(0, 1);	
	for(idx_type i=0; i<n; ++i)
	{
		Vector newPoint = Vector(d);
		for(idx_type j=0; j<d; ++j)
			newPoint(j) = urd(gen)*interval(j);
		dataset.points.col(old_n+i) = newPoint+anchor;
		dataset.weights(old_n+i)=1.;
	}		
}

Vector datagenutil::getExpIncrNumbers(idx_type k, fp_type exponent, idx_type number)
{
	Vector output = datagenutil::getExpIncrWeights(k, exponent);
	output = output*number;
	for(idx_type i=0; i<k-1; ++i)
		output(i) = round(output(i));
	output(k) = 0;
	output(k) = number - output.sum();
	return output;
}

Vector datagenutil::getExpIncrWeights(idx_type k, fp_type exponent)
{
	Vector output = Vector(k);
	fp_type factor = pow(2,exponent);
	output(0) = 1.;
	for(idx_type i=1; i<k; ++i)
	  output(i) = output(i-1)*factor;	
	output = output/output.sum();
	return output;
}

datagenutil::BoundingBox datagenutil::boundingBox(commonutil::DataSet const& dataset)
{
	const idx_type n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	
	if(n==0||d==0)
		gmmlab_throw("datagenutil::boundingBox() - Given dataset must not be empty.");
		
	BoundingBox boundingBox;
	boundingBox.min = dataset.points.col(0);
	boundingBox.max = boundingBox.min;
	for(idx_type i=1; i<n; ++i)
	{
		boundingBox.min = boundingBox.min.cwiseMin(dataset.points.col(i));
		boundingBox.max = boundingBox.max.cwiseMax(dataset.points.col(i));
	}
	return boundingBox;
}
