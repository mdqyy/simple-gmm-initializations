#ifndef ALTERNATING_ADAPTIVEMEANS_MEANS2GMM_ID_H
#define ALTERNATING_ADAPTIVEMEANS_MEANS2GMM_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeansAndMeans2GMMID : public InitDescriptor
{
public:
	AlternatingAdaptiveMeansAndMeans2GMMID(uint32_t seed, fp_type factor);

	virtual ~AlternatingAdaptiveMeansAndMeans2GMMID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	fp_type factor;
};

#endif
