#include "agglomerative_init_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"
#include "../../../settings/configparser.h"
#include <boost/algorithm/string.hpp>

AgglomerativeInitID::AgglomerativeInitID(fp_type sampleSizeFactor, uint32_t s) : seed(s), sampleSizeFactor(sampleSizeFactor)
{
	std::stringstream sstream;
	sstream << "AgglomerativeInit_" << "_sample" << sampleSizeFactor << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc AgglomerativeInitID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	fp_type sampleSizeFactor = this->sampleSizeFactor;
	if(this->sampleSizeFactor == 0)
		sampleSizeFactor = 1./k;
	return initutil::sampleAgglomerativeMeansToGMM(input, k, sampleSizeFactor, 0, gen);
}
