#ifndef UNIFORM_MEANS2GMM_ID_H
#define UNIFORM_MEANS2GMM_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class UniformMeans2GMMID : public InitDescriptor
{
public:
	UniformMeans2GMMID(uint32_t = 1);

	virtual ~UniformMeans2GMMID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
};

#endif
