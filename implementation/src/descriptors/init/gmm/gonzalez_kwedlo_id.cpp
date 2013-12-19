#include "gonzalez_kwedlo_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

GonzalezKwedlo::GonzalezKwedlo(uint32_t s, fp_type sampleSizeFactor) : seed(s), sampleSizeFactor(sampleSizeFactor)
{
	std::stringstream sstream;
	sstream << "GonzalezKwedlo" << "_i" << s <<"_sample" << sampleSizeFactor ;
	this->nametag = sstream.str();
}

GMMDesc GonzalezKwedlo::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	//(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor)
	return initutil::gonzalezKwedlo(input, k, gen, sampleSizeFactor);
}
