#ifndef EMALGODESCRIPTOR_H
#define EMALGODESCRIPTOR_H

#include "../algodescriptor.h"

#include <vector>

/**
 * @brief Interface for describing algorithms to be run
 */
class EMAlgoDescriptor : public AlgoDescriptor
{
public:
	EMAlgoDescriptor(unsigned int, uint32_t);

	virtual ~EMAlgoDescriptor()
	{
	}
	
	GMMAlgorithm* create(commonutil::DataSet const&) const;

protected:
	uint32_t seed;
};

#endif
