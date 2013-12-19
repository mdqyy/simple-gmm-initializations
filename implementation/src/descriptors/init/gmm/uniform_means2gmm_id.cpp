#include "uniform_means2gmm_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

UniformMeans2GMMID::UniformMeans2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "UniformMeans2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc UniformMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::uniformMeansToGMM(input, k, gen);
}
