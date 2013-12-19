#ifndef AGGLOMERATIVE_INIT_ID_H
#define AGGLOMERATIVE_INIT_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class AgglomerativeInitID : public InitDescriptor
{
public:
	AgglomerativeInitID(fp_type sampleSizeFactor, uint32_t = 1);

	virtual ~AgglomerativeInitID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	fp_type sampleSizeFactor;
};

#endif
