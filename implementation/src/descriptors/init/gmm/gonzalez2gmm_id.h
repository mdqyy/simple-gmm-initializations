#ifndef GONZALEZ2GMM_ID_H
#define GONZALEZ2GMM_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class Gonzalez2GMMID : public InitDescriptor
{
public:
	Gonzalez2GMMID(uint32_t = 1);

	virtual ~Gonzalez2GMMID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
};

#endif
