#ifndef ADAPTIVE_SPHERICAL_ID_H
#define ADAPTIVE_SPHERICAL_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization inspired by Dasgupta & Schulman
 */
class AdaptiveSphericalID : public InitDescriptor
{
public:
	AdaptiveSphericalID(uint32_t = 1);

	virtual ~AdaptiveSphericalID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
};

#endif
