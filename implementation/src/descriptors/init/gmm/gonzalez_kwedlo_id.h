#ifndef GONZALEZ_KWEDLO_ID_H
#define GONZALEZ_KWEDLO_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class GonzalezKwedlo : public InitDescriptor
{
public:
	GonzalezKwedlo(uint32_t seed, fp_type sampleSizeFactor);

	virtual ~GonzalezKwedlo()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	fp_type sampleSizeFactor;
};

#endif
