#include "uniform_spherical_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

UniformSphericalID::UniformSphericalID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "UniformSpherical" << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc UniformSphericalID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::uniformSphericalGMM(input, k, gen);
}
