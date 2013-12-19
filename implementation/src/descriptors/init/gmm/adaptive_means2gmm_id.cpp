#include "adaptive_means2gmm_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

AdaptiveMeans2GMMID::AdaptiveMeans2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "AdaptiveMeans2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc AdaptiveMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::adaptiveMeansToGMM(input, k, gen);
}
