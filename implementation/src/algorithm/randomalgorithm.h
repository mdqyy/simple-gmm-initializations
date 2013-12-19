#ifndef RANDOMALGORITHM_H
#define RANDOMALGORITHM_H

#include "gmmalgorithm.h"

#include "../base/gmmdesc.h"
#include "../base/commonutil.h"

#include <set>
#include <ctime>
#include <vector>

/**
* @brief BBKSS algorithm for GMMs
*
* @ingroup learning_algorithms
*/
class RandomAlgorithm : public GMMAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	RandomAlgorithm(commonutil::DataSet const&, bool, uint32_t);

	virtual ~RandomAlgorithm()
	{
	}

	/**
	* @brief does a dynamic cast of the given GMMAlgorithm to SamplingAlgorithm
	* @return NULL if the GMMAlgorithm is not a SamplingAlgorithm instance
	*/
	static RandomAlgorithm* toRandomAlgorithm(GMMAlgorithm* a);

protected:
	uint32_t seed;
	std::mt19937 gen;
};

#endif
