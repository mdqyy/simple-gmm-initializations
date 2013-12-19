#ifndef UNIFORM_SPHERICAL_ID_H
#define UNIFORM_SPHERICAL_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization by Dasgupta & Schulman
 */
class UniformSphericalID : public InitDescriptor
{
public:
	UniformSphericalID(uint32_t = 1);

	virtual ~UniformSphericalID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
};

#endif
