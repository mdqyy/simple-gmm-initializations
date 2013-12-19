#include "adaptive_spherical_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

AdaptiveSphericalID::AdaptiveSphericalID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "AdaptiveSpherical" << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc AdaptiveSphericalID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::adaptiveSphericalGMM(input, k, gen);
}
