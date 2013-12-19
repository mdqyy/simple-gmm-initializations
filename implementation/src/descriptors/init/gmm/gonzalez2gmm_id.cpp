#include "gonzalez2gmm_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

Gonzalez2GMMID::Gonzalez2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "Gonzalez2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc Gonzalez2GMMID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::gonzalezToGMM(input, k, gen);
}
