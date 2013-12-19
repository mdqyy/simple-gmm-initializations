#ifndef UNIFORM_SPHERCIAL_WITH_PRUNING_ID_H
#define UNIFORM_SPHERCIAL_WITH_PRUNING_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization by Dasgupta & Schulman with pruning!
 */
class UniformSphericalWithPruningID : public InitDescriptor
{
public:
	UniformSphericalWithPruningID(unsigned int oversamplingFactor, uint32_t seed = 1);

	virtual ~UniformSphericalWithPruningID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	unsigned int oversamplingFactor;
};

#endif
