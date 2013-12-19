#include "alternating_adaptivemeans_means2gmm_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

AlternatingAdaptiveMeansAndMeans2GMMID::AlternatingAdaptiveMeansAndMeans2GMMID(uint32_t s, fp_type factor) : seed(s), factor(factor)
{
	std::stringstream sstream;
	sstream << "Alt_AdaptMean+Means2GMM"<< "_f" << factor << "_i" << s ;
	this->nametag = sstream.str();
}

GMMDesc AlternatingAdaptiveMeansAndMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor);
}
