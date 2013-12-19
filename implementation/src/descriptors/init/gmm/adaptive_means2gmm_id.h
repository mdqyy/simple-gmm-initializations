#ifndef ADAPTIVE_MEANS2GMM_ID_H
#define ADAPTIVE_MEANS2GMM_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class AdaptiveMeans2GMMID : public InitDescriptor
{
public:
	AdaptiveMeans2GMMID(uint32_t = 1);

	virtual ~AdaptiveMeans2GMMID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
};

#endif
