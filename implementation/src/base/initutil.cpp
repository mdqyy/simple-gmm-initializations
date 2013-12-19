#include "initutil.h"

#include "gmmutil.h"
#include "linalgutil.h"
#include "kmeansutil.h"
#include "initutil.h"
#include "commonutil.h"
#include "similarities.h"

#include "../algorithm/gmm/emforgmm.h"
#include "../algorithm/gmmalgorithm.h"

void initutil::check(idx_type k, idx_type d, idx_type n)
{
	if(k==0)
		gmmlab_throw("initutil::check() - Number of means is zero.")
	if(d==0 || n==0)
		gmmlab_throw("initutil::check() - Input is empty.");
}


GMMDesc initutil::restartInit(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, InitMethod initmethod, idx_type subruns, idx_type substeps, fp_type factor)
{
	initutil::check(k, input.points.rows(), input.points.cols());
	
	if(subruns == 0)
		subruns = k;
	
	GMMAlgorithm* em = new EMforGMM(input, 0, gen);
	GMMDesc tmp;
	GMMDesc bestGMMDesc;
	fp_type minCosts;
	for(idx_type t=0; t < subruns; ++t)
	{
		switch(initmethod){
			case UNIFORM_MEANS2GMM:
				tmp = initutil::uniformMeansToGMM(input, k, gen);
				break;
			case ADAPTIVE_MEANS2GMM:
				tmp = initutil::adaptiveMeansToGMM(input, k, gen);
				break;
			case GONZALEZ2GMM:
				tmp = initutil::gonzalezToGMM(input, k, gen);
				break;
			case UNIFORM_SPHERICAL:
				tmp = initutil::uniformSphericalGMM(input, k, gen);
				break;
			case ADAPTIVE_SPHERICAL:
				tmp = initutil::adaptiveSphericalGMM(input, k, gen);
				break;
			case ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM:
				tmp = initutil::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor);
				break;
			default:
				gmmlab_throw("initutil::restartInit() - Unsupported initmethod.");
		}
		
		em->init(tmp);
		em->run(substeps);
		tmp = em->getGMMDesc();
		fp_type tmpCosts = gmmutil::nll(input, tmp);
		if(t==0 || (tmpCosts < minCosts))
		{
			minCosts = tmpCosts;
			bestGMMDesc = tmp;
		}
	}
	delete em;
	
	//std::cout << "restartUniformMeans2GMMEM = \n" << bestGMMDesc << std::endl;
	
	return bestGMMDesc;
}


Matrix initutil::uniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	initutil::check(k, d, n);
	
	Matrix means(d,k);
	for (idx_type i=0; i<k; ++i)
		means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));

	return means;
}

GMMDesc initutil::uniformMeansToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{
	return kmeansutil::meansToGMM(input, uniformMeans(input,k,gen), computeWeightAndCovar);
}

GMMDesc initutil::wrappedUniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{
	return kmeansutil::wrapMeans(input, uniformMeans(input,k,gen), computeWeightAndCovar);
}


Matrix initutil::adaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	initutil::check(k, d, n);
	
	Matrix means(d,k);
	Vector sqNorms;
	for (idx_type i=0; i<k; ++i)
		if (i==0)
			means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));
		else
		{
			if (i==1)
				sqNorms = (input.points.colwise()-means.col(0)).colwise().squaredNorm();
			else
			{
				for (idx_type j=0; j<n; ++j)
				{
					fp_type sqn = (input.points.col(j)-means.col(i-1)).squaredNorm();
					if (sqn<sqNorms[j])
						sqNorms[j]=sqn;
				}
			}

			means.col(i) = input.points.col(commonutil::randomIndex(sqNorms, gen));
		}

	return means;
}


Matrix initutil::gonzalez(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen)
{	
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();
	
	initutil::check(k, d, n);
	
	Matrix means(d,k);
	Vector sqNorms;
	for (idx_type i=0; i<k; ++i)
		if (i==0)
			means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));
		else
		{
			if (i==1)
				sqNorms = (input.points.colwise()-means.col(0)).colwise().squaredNorm();
			else
			{
				for (idx_type j=0; j<n; ++j)
				{
					fp_type sqn = (input.points.col(j)-means.col(i-1)).squaredNorm();
					if (sqn<sqNorms[j])
						sqNorms[j]=sqn;
				}
			}
			
			idx_type index;
			sqNorms.maxCoeff(&index);
			means.col(i) = input.points.col(index);
		}
		
		return means;
}

GMMDesc initutil::gonzalezToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{	
	return kmeansutil::meansToGMM(input, initutil::gonzalez(input, k, gen));
}

GMMDesc initutil::wrappedGonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{	
	return kmeansutil::wrapMeans(input, initutil::gonzalez(input, k, gen), computeWeightAndCovar);
}

GMMDesc initutil::adaptiveMeansToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	return kmeansutil::meansToGMM(input, adaptiveMeans(input,k,gen));
}

GMMDesc initutil::wrappedAdaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar)
{
	return kmeansutil::wrapMeans(input, adaptiveMeans(input,k,gen), computeWeightAndCovar);
}


GMMDesc initutil::alternatingAdaptiveMeanAndMeanToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor)
{  
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initutil::check(k, d, n);

	GMMDesc desc;
	Vector sample;
	
	for (idx_type i=0; i<k; ++i)
	{
		if (i==0)
			// draw first point uniformly
			sample = input.points.col(commonutil::randomIndex(input.weights, gen)); 
			// ... yeah, it's useless, but now it stays here because otherwise the seeds are screwed up and we get different results
		else
		{
			// draw next point w.r.t. current mixture
			Vector densities = gmmutil::adaptiveDensities(input, desc, factor);
			sample = input.points.col(commonutil::randomIndex(densities, gen)); 
		}
		
		Matrix tmpMeans = Matrix::Zero(d,i+1);
		if(i>0)
			tmpMeans.block(0,0,d,i) = desc.means;
		tmpMeans.col(i) = sample;
		
		desc = kmeansutil::meansToGMM(input, tmpMeans);
	}

	return desc;
}



GMMDesc initutil::gonzalezForGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool use2GMM, fp_type sampleFactor)
{  
	
// 	std::cout << "sampleFactor = " << sampleFactor << std::endl;	
// 	std::cout << "use2GMM = " << use2GMM << std::endl;	
	
	idx_type sampleSize = sampleFactor * input.points.cols();
	
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initutil::check(k, d, n);

	GMMDesc desc;
	Vector sample;
	
	fp_type trace;
	if(!use2GMM)
	{
		GMMDesc gmmdesc = gmmutil::optimalGaussian(input);
		Matrix covar = gmmdesc.covariances.at(0);
		trace = covar.trace()/(10.*d*k);
		
		desc.means.resize(0,0);
		desc.weights.resize(0);

	}
		
	std::uniform_real_distribution<> urd(0,1);
	std::uniform_int_distribution<> uid(0,input.points.cols()-1);
	
	Matrix randMatrix(d,d);
	Matrix randOrthonormalMatrix(d,d);
	Vector eigenvalues(d);
	
	for (idx_type i=0; i<k; ++i)
	{
		
		if (i==0)
			// draw first point uniformly
			sample = input.points.col(commonutil::randomIndex(input.weights, gen)); 
		else
		{
			// draw next point w.r.t. current mixture
			
			if(sampleFactor < 1)
			{
				commonutil::DataSet samples;
				samples.points.resize(d,sampleSize);
				samples.weights.resize(sampleSize);
				idx_type index;
				for(idx_type j=0; j<sampleSize; ++j){
					index = uid(gen);
				
// std::cout << "index = " << index << std::endl;
					
					samples.points.col(j) = input.points.col(index);
					samples.weights(j) = input.weights(index);
				}
				Vector densities = gmmutil::adaptiveDensities(samples, desc, 1.);	
				
// std::cout << "densities.size() = " << densities.size() << std::endl;

				densities.maxCoeff(&index);
				sample = samples.points.col(index); 

			}
			else
			{
				idx_type index;
				Vector densities = gmmutil::adaptiveDensities(input, desc, 1.);	

// std::cout << "densities.size() = " << densities.size() << std::endl;

				densities.maxCoeff(&index);
				sample = input.points.col(index); 
			}
			
		}
		
		if(use2GMM)
		{
			Matrix tmpMeans = Matrix::Zero(d,i+1);
			if(i>0)
				tmpMeans.block(0,0,d,i) = desc.means;
			tmpMeans.col(i) = sample;
			
			desc = kmeansutil::meansToGMM(input, tmpMeans);
		}
		else
		{
			desc.means.conservativeResize(d,i+1);
			desc.means.col(i) = sample;
			
			desc.weights.conservativeResize(i+1);
			
			// random covariance
			for (idx_type i=0; i<d; ++i)
				for (idx_type j=0; j<d; ++j)
					randMatrix(i,j)=urd(gen);
			Eigen::HouseholderQR<Matrix> qr(randMatrix);
			randOrthonormalMatrix = qr.householderQ();
			commonutil::fill(eigenvalues,gen, 1.,10.);
			fp_type tmp = trace/eigenvalues.sum();
			eigenvalues = tmp * eigenvalues;
			randMatrix = randOrthonormalMatrix.transpose()*eigenvalues.asDiagonal()*randOrthonormalMatrix;
			
			desc.covariances.push_back(randMatrix);
		}
	}
	
	if(!use2GMM)
	{
		// estimate weights
		desc.weights = Vector::Zero(k);
		std::vector<idx_type> partition = gmmutil::gaussPartition(input.points, desc);
		for(idx_type n=0; n<input.points.cols(); ++n)
			++desc.weights(partition.at(n));		
		desc.weights /= input.points.cols();
	}

	return desc;
}



GMMDesc initutil::gonzalezKwedlo(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor)
{  
	
// 	std::cout << "sampleFactor = " << sampleFactor << std::endl;	
// 	std::cout << "use2GMM = " << use2GMM << std::endl;	
	
	idx_type sampleSize = sampleFactor * input.points.cols();
	
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initutil::check(k, d, n);

	GMMDesc desc;
	Vector sample;

	
	GMMDesc gmmdesc = gmmutil::optimalGaussian(input);
	Matrix covar = gmmdesc.covariances.at(0);
	fp_type trace = covar.trace()/(10.*d*k);	
	desc.means.resize(0,0);
	desc.weights.resize(0);

		
	std::uniform_real_distribution<> urd(0,1);
	std::uniform_int_distribution<> uid(0,input.points.cols()-1);
	
	Matrix randMatrix(d,d);
	Matrix randOrthonormalMatrix(d,d);
	Vector eigenvalues(d);
	
	for (idx_type i=0; i<k; ++i)
	{
		
		if (i==0)
			// draw first point uniformly
			sample = input.points.col(commonutil::randomIndex(input.weights, gen)); 
		else
		{
			// draw next point w.r.t. current mixture
			
			if(sampleFactor < 1)
			{
				commonutil::DataSet samples;
				samples.points.resize(d,sampleSize);
				samples.weights.resize(sampleSize);
				idx_type index;
				for(idx_type j=0; j<sampleSize; ++j){
					index = uid(gen);
				
// std::cout << "index = " << index << std::endl;
					
					samples.points.col(j) = input.points.col(index);
					samples.weights(j) = input.weights(index);
				}
				Vector densities = gmmutil::adaptiveDensities(samples, desc, 1.);	
				
// std::cout << "densities.size() = " << densities.size() << std::endl;

				densities.maxCoeff(&index);
				sample = samples.points.col(index); 

			}
			else
			{
				idx_type index;
				Vector densities = gmmutil::adaptiveDensities(input, desc, 1.);	

// std::cout << "densities.size() = " << densities.size() << std::endl;

				densities.maxCoeff(&index);
				sample = input.points.col(index); 
			}
			
		}
		
		desc.means.conservativeResize(d,i+1);
		desc.means.col(i) = sample;
		
		desc.weights.conservativeResize(i+1);
		desc.weights(i) = urd(gen);
			
		// random covariance
		for (idx_type i=0; i<d; ++i)
			for (idx_type j=0; j<d; ++j)
				randMatrix(i,j)=urd(gen);
		Eigen::HouseholderQR<Matrix> qr(randMatrix);
		randOrthonormalMatrix = qr.householderQ();
		commonutil::fill(eigenvalues,gen, 1.,10.);
		fp_type tmp = trace/eigenvalues.sum();
		eigenvalues = tmp * eigenvalues;
		randMatrix = randOrthonormalMatrix.transpose()*eigenvalues.asDiagonal()*randOrthonormalMatrix;
			
		desc.covariances.push_back(randMatrix);
		
	}
	

	desc.weights /= desc.weights.sum();

	return desc;
}

GMMDesc initutil::uniformSphericalGMMwithPruning(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int oversamplingFactor)
{
	initutil::check(k, input.points.rows(), input.points.cols());
	
    idx_type oversamplingk = (idx_type)(round(oversamplingFactor*1.*k*log(oversamplingFactor*1.*k)));
		
	if(oversamplingk < k)
	  gmmlab_throw("initutil::uniformSphericalGMMwithPruning() - oversamplingk should be larger than k.");
	
	GMMDesc oversampled = initutil::uniformSphericalGMM(input, oversamplingk, gen);
	return initutil::prune(oversampled, k, gen);
}


GMMDesc initutil::prune(GMMDesc& desc, idx_type k, std::mt19937& gen)
{
	if(desc.weights.size() == k)
		return desc;
	if(desc.weights.size() < k)
		gmmlab_throw("gmmutil::prune - gmm has less than k components");
	
	idx_type d = desc.means.rows();

	// 1. remove center estimates whose mixing weights are below w_T = 1/(4k) 
	idx_type desc_k = desc.weights.size();
	idx_type index;
	while(desc.weights.size()>k)
	{
		if(desc.weights.minCoeff(&index) < 1/(4.*desc_k))
		{
			commonutil::erase(desc.weights,index);
			desc.covariances.erase(desc.covariances.begin()+index);
			commonutil::eraseColumn(desc.means, index);
		}
		else
			break;
	}

	desc_k = desc.weights.size();
	if (desc_k>k)
	{
		GMMDesc prunedGMM;
		prunedGMM.weights = Vector::Constant(k, 1.0/k);
		prunedGMM.means = Matrix::Zero(d, k);
		prunedGMM.covariances.clear();

		// 2. prune the center estimates (adaption of Gonzales (1985))
		//    and set the mixing weights to w_i = 1/k

		// 2.a Choose one of these centers arbitrarily
		std::uniform_int_distribution<> uid(0, desc.weights.size()-1);
		index = uid(gen);

		prunedGMM.means.col(0) = desc.means.col(index);
		prunedGMM.covariances.push_back(desc.covariances.at(index));

		// 2.b Pick the remaining k-1 centers iteratively: pick the center farthest from the ones picked so far,
		//     i.e., pick he point x with $min_{y\in S} dist(x,y)$
		fp_type diffcovar, dist;
		for(idx_type i=1; i<k; ++i)
		{
			Vector distances = Vector(desc_k);
			for (idx_type j=0; j<desc_k; ++j)
			{
				diffcovar = desc.covariances.at(j)(0,0) - prunedGMM.covariances.at(0)(0,0);
				if(diffcovar <= 0)
				  diffcovar = FP_INFINITE;
				  //gmmlab_throw("initutil::prune() - difference of the covariances is <= 0.");
				distances[j] = (desc.means.col(j)-prunedGMM.means.col(0)).norm()/diffcovar;
				
				for (idx_type l=1; l<i; ++l)
				{
				  	diffcovar = desc.covariances.at(j)(0,0) - prunedGMM.covariances.at(l)(0,0);
					if(diffcovar <= 0)
					  diffcovar = FP_INFINITE;
					  //gmmlab_throw("initutil::prune() - difference of the covariances is <= 0.");
					dist = (desc.means.col(j)-prunedGMM.means.col(l)).norm()/diffcovar;

					if (dist<distances[j])
						distances[j]=dist;
				}
			}

			// pick the farthest
			distances.maxCoeff(&index);

			// insert into pruned GMM
			prunedGMM.means.col(i)=desc.means.col(index);
			prunedGMM.covariances.push_back(desc.covariances.at(index));
		}

		
		//std::cout << "prunedGMM = \n" << prunedGMM << std::endl;
		
		return prunedGMM;
	}
	return desc;
}


GMMDesc initutil::uniformSphericalGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type n = input.points.cols();
	idx_type d = input.points.rows();
	
	initutil::check(k, d, n);

	GMMDesc desc;
	desc.weights = Vector::Constant(k, 1.0/k);
	desc.covariances = std::vector<Matrix>(k, (1.0/(2*d))*Matrix::Identity(d, d));

	desc.means = Matrix(d,k);
	for (size_t i=0; i<k; ++i)
		desc.means.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));

	for (size_t i=0; i<k; ++i)
	{
		Vector sq = (desc.means.colwise()-desc.means.col(i)).colwise().squaredNorm();

		fp_type min = 0;
		for (size_t j=0; j<k; ++j)
			if (j!=i && sq[j]>0 && (min==0 || sq[j]<min))
				min = sq[j];
		if (min<=0)
			min = 1;

		desc.covariances[i] *= min;
	}

	return desc;
}

GMMDesc initutil::adaptiveSphericalGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initutil::check(k, d, n);

	Matrix samples(d,k+1);	// two samples for the first gaussian, one for each of the others
	Vector sqDist(k);		// minimum distance to nearest sample (for each sample)

	GMMDesc desc;
	for (idx_type i=0; i<k+1; ++i)
	{
		if (i<2)
			// draw the first two points uniformly
			samples.col(i) = input.points.col(commonutil::randomIndex(input.weights, gen));
		else
		{
			// draw next point w.r.t. current mixture
			Vector densities = gmmutil::adaptiveDensities(input, desc, 1);
			samples.col(i) = input.points.col(commonutil::randomIndex(densities, gen));
		}

//		std::cout << "sampled points = " << samples.col(i).transpose() << std::endl;

		if (i>0)
		{
			desc.weights = Vector::Constant(i, 1.0/i);
			desc.means = samples.block(0,1,d,i);
			desc.covariances = std::vector<Matrix>(i, (1.0/(2*d))*Matrix::Identity(d,d));

			if (i==1)	// after uniformly drawing two samples create the first model with one component
			{
				sqDist[0] = (samples.col(0)-samples.col(1)).squaredNorm();
				if (sqDist[0]<=0)
					sqDist[0] = 1;
			}
			else if (i==2)
			{
				sqDist[0] = (samples.col(1)-samples.col(2)).squaredNorm();
				if (sqDist[0]<=0)
					sqDist[0] = 1;
				sqDist[1] = sqDist[0];
			}
			else if (i>2)
			{
				Vector mins = (samples.block(0,1,d,i-1).colwise()-samples.col(i)).colwise().squaredNorm();

				// set sqDist to the minimum non-zero distance to one of the other samples
				// this is due to the possibility of sampling the same point multiple times
				sqDist[i-1] = 0;
				for (idx_type j=0; j<i-1; ++j)
					if (mins[j]>0 && (sqDist[i-1]==0 || mins[j]<sqDist[i-1]))
						sqDist[i-1] = mins[j];
				if (sqDist[i-1]<=0)
					sqDist[i-1] = 1;	// prevent zero covariances in pathological cases

				for (idx_type j=0; j<i-1; ++j)
					if (mins[j]>0 && mins[j]<sqDist[j])
						sqDist[j] = mins[j];
			}

			for (idx_type j=0; j<i; ++j)
				desc.covariances[j] *= sqDist[j];
		}
	}

	return desc;
}


GMMDesc initutil::sampleAgglomerativeMeansToGMM(const commonutil::DataSet& input, const idx_type k, const fp_type sampleSizeFactor, const bool precompute, std::mt19937& gen)
{
	idx_type d = input.points.rows();
	idx_type sampleSize = (int) round(input.points.cols()*sampleSizeFactor);
	
	// sample points uniformly at random
	commonutil::DataSet sample;
	sample.points.resize(d,sampleSize);
	sample.weights.resize(sampleSize);
	std::uniform_int_distribution<> uid(0, input.points.cols()-1);
	for(idx_type i=0; i<sampleSize; ++i)
	{
		idx_type rndIndex = uid(gen);
		sample.points.col(i) = input.points.col(rndIndex);
		sample.weights(i) = input.weights(rndIndex);
	}
	
	// start agglomerativeMeans
	Matrix means = initutil::agglomerativeMeans(sample, k, precompute);
	
	return kmeansutil::meansToGMM(input, means);
}


// NOTE: only works with unweighted points
Matrix initutil::agglomerativeMeans(commonutil::DataSet const& input, const idx_type k, const bool precompute)
{
	GMMDesc desc;

	if (k==0)
		gmmlab_throw("initutil::agglomerativeToGMM() - Empty mixture requested.");

	idx_type n = input.points.cols();
	idx_type d = input.points.rows();

	if (n==0 || d==0)
		gmmlab_throw("initutil::agglomerativeToGMM() - Input is empty.");


	// initialize partition with trivial n-clustering
	std::vector<std::vector<idx_type>> partition;
	std::vector<Vector> partition_sums;
	for(idx_type i=0; i<n; ++i){
	  std::vector<idx_type> cluster;
	  cluster.push_back(i);
	  partition.push_back(cluster);
	  partition_sums.push_back(input.points.col(i));
	}

	// distance measure
	AverageLinkage averageLinkageDis;
	

	// precompute dissimilarities
	Matrix dis;
	if (precompute)
	{
		dis = Matrix::Zero(n,n);
		for (idx_type i=0; i<n; ++i)	
			for (idx_type j=0; j<=i; ++j)
				dis(i,j) = averageLinkageDis(partition_sums.at(i), partition.at(i).size(), partition_sums.at(j), partition.at(j).size());
	}

	for (idx_type r=n; r>k; --r)
	{				
		idx_type first = -1;
		idx_type second = -1;
		fp_type min = FP_INFINITE;
		for (idx_type i=0; i<r; ++i)
			for (idx_type j=0; j<i; ++j) // j < i
			{
				//std::cout << "i=" << i << ", j=" << j << std::endl;
				if (i!=j)
				{
					fp_type d;
					if (precompute)
						d = dis(i,j);
					else
						d = averageLinkageDis(partition_sums.at(i), partition.at(i).size(), partition_sums.at(j), partition.at(j).size());

					if ((i==1 && j==0) || d < min)
					{
						min = d;
						first = j;
						second = i;
					}
				}
			}
		
		// merge clusters (note: first < second)
		// 1. update sufficient statistics 
		partition_sums.at(first) = partition_sums.at(first)+partition_sums.at(second);
		// 2. move points from second cluster to the first
		while(!partition.at(second).empty())
		{
		  idx_type point = partition.at(second).back();
		  partition.at(second).pop_back();
		  partition.at(first).push_back(point);
		} 
		// overwrite second cluster with the last stored cluster and remove the last cluster  (note: first < second)
		partition.at(second) = partition.back();
		partition_sums.at(second) = partition_sums.back();
		partition.pop_back();
		partition_sums.pop_back();
		
		if (precompute)
		{
			// cluster with index second has been overwritten by the last cluster (with index r-1)
			for(idx_type i=0; i<r-1; i++)
			{
				if(i!=second)
				{
					fp_type d = dis(r-1,i);
					if(i < second)
						dis(second, i) = d;
					else if(second < i)
						dis(i,second) =  d;
				}
			}
			
			// compute distances wrt the newly formed cluster which is stored at index first
			for (idx_type i=0; i<r-1; i++)
			{
				if(i!=first){
					fp_type d = averageLinkageDis(partition_sums.at(first), partition.at(first).size(), partition_sums.at(i), partition.at(i).size());
					if(i < first)
						dis(first, i) = d;
					else if(first < i)
						dis(i,first) =  d;
				}
			}
		}
	}
	
	Matrix means(d,k);
	for(idx_type i=0; i<k; ++i)
		means.col(i) = partition_sums.at(i) / partition.at(i).size();
	
	
	return means;
}
