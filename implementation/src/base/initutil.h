#ifndef INITUTIL_H
#define INITUTIL_H

#include "base.h"
#include "gmmdesc.h"
#include "commonutil.h"

#include <iostream>

namespace initutil
{


	/**
	 * chooses k means uniformly at random from the given input set.
	 */
	Matrix uniformMeans(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen);

	
	/**
	 * computes means according to initutil::uniformMeans, then uses these means to
	 * compute a GMM by using kmeans::meansToGMM.
	 */
	GMMDesc uniformMeansToGMM(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen, bool computeWeightAndCovar = 1);
	
	/**
	 * computes means according to initutil::uniformMeans, then uses these means
	 * as centers of a GMM by using kmeans::wrapMeans. 
	 */
	GMMDesc wrappedUniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar = 1);
	
	/**
	 * chooses k means according to kmeans++.
	 */
	Matrix adaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	
	/**
	 * computes means according to initutil::adaptiveMeans, then uses these means to
	 * compute a GMM by using kmeansutil::meansToGMM.
	 */
	GMMDesc adaptiveMeansToGMM(commonutil::DataSet const&, idx_type, std::mt19937& gen);
	
	/**
	 * computes means according to initutil::adaptiveMeans, then uses these means
	 * as centers of a GMM by using kmeansutil::wrapMeans.
	 */
	GMMDesc wrappedAdaptiveMeans(commonutil::DataSet const&, idx_type, std::mt19937& gen, bool computeWeightAndCovar = 1);

	/**
	 * initialization according to Dasgupta&Schulman but without pruning.
	 */
	GMMDesc uniformSphericalGMM(commonutil::DataSet const&, idx_type, std::mt19937& gen);
	
	/**
	 * initialization according to Dasgupta&Schulman.
	 * 
	 * @param oversamplingFactor determines how many components are sampled before pruning, i.e., oversamplingFactor*k*ln(oversamplingFactor*k)
	 */
	GMMDesc uniformSphericalGMMwithPruning(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int oversamplingFactor = 2);

	/**
	 * Gonzales-like pruning according to Dasgupta&Schulman.
	 */
	GMMDesc prune(GMMDesc& desc, idx_type k, std::mt19937& gen);
  
	/**
	 * adaptive initialization similar to Dasgupta&Schulman.
	 */
	GMMDesc adaptiveSphericalGMM(commonutil::DataSet const&, idx_type, std::mt19937& gen);

	
	/**
	 * alternately performs two steps: Chooses a point wrt gmmutil::adaptiveDensities and
	 * calls kmeansutil::meansToGMM wrt the sampled point and the means of the current gmm.
	 */
	GMMDesc alternatingAdaptiveMeanAndMeanToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor = 1.);

	
	/**
	 * runs agglomerativeMeans and calls kmeansuitl::meansToGMM wrt the input set and the returned means.
	 */
	GMMDesc sampleAgglomerativeMeansToGMM(const commonutil::DataSet& input, const idx_type k, const fp_type sampleSizeFactor, const bool precompute, std::mt19937& gen);

	
	/**
	 * computes average linkage clustering. 
	 * NOTE: only works with unweighted points
	 */
	Matrix agglomerativeMeans(commonutil::DataSet const&, idx_type, bool = false);
	
	/**
	 * computes 2*d adaptive initializations according to kmeans++ (independent of each other),
	 * then samples k centers from these 2*d*k sampled centers according to kmeans++
	 */
	Matrix resampleAdaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	
	/**
	 * first computes 2*d adaptive spherical initializations (independent from each other), 
	 * then computes an adaptive spherical initialization with respect to the means of these
	 * initializations.
	 */
	GMMDesc resampleGMMInit(commonutil::DataSet const& input, InitMethod initmethod, idx_type k, std::mt19937& gen);
	
	/**
	 * restarts the EM-Algorithm as many times as defined by subruns (if subruns == 0, it will be set to k). 
	 * Each EM instance is initialized via initmethod and runs for substeps steps.
	 * Returns the resulting GMMDesc with the smallest cost. 
	 */
	GMMDesc restartInit(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, InitMethod initmethod, idx_type subruns, idx_type substeps, fp_type factor = 1.);
	
	/**
	 * implements Gonazalez algorithm.
	 */
	Matrix gonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	
	GMMDesc gonzalezForGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool use2GMM, fp_type sampleFactor);
	
	GMMDesc gonzalezKwedlo(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor);
	
	/**
	 * runs Gonzalez algorithm and wraps the Means using kmeansutil::wrapMeans.
	 */
	GMMDesc wrappedGonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar);
	
	/**
	 * runs Gonzalez algorithm and calls kmeansutil::meansToGMM.
	 */
	GMMDesc gonzalezToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	
	
	/**
	 * throws an exception if k==0 || d==0 || n==0.
	 */
	void check(idx_type k, idx_type d, idx_type n);
	
}

#endif