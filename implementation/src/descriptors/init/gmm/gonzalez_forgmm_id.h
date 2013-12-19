#ifndef GONZALEZ_FORGMM_ID_H
#define GONZALEZ_FORGMM_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class GonzalezForGMMID : public InitDescriptor
{
public:
	GonzalezForGMMID(uint32_t seed, bool useGMM, fp_type sampleSizeFactor);

	virtual ~GonzalezForGMMID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	bool use2GMM;
	fp_type sampleSizeFactor;
};

#endif
