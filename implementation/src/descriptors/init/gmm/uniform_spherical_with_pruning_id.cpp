#include "uniform_spherical_with_pruning_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"

UniformSphericalWithPruningID::UniformSphericalWithPruningID(unsigned int oversamplingFactor, uint32_t s) : seed(s), oversamplingFactor(oversamplingFactor)
{
	std::stringstream sstream;
	sstream << "UniformSphericalWithPruning " << "_o" << oversamplingFactor << "_i" << s ;
	this->nametag = sstream.str();
}

GMMDesc UniformSphericalWithPruningID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::uniformSphericalGMMwithPruning(input, k, gen, oversamplingFactor);
}
