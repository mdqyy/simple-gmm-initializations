#include "gonzalez_forgmm_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

GonzalezForGMMID::GonzalezForGMMID(uint32_t s, bool use2GMM, fp_type sampleSizeFactor) : seed(s), use2GMM(use2GMM), sampleSizeFactor(sampleSizeFactor)
{
	std::stringstream sstream;
	sstream << "GonzalezForGMM" << "_i" << s << "_2gmm" << use2GMM <<"_sample" << sampleSizeFactor ;
	this->nametag = sstream.str();
}

GMMDesc GonzalezForGMMID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	//(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool use2GMM, fp_type sampleFactor)
	return initutil::gonzalezForGMM(input, k, gen, use2GMM, sampleSizeFactor);
}
