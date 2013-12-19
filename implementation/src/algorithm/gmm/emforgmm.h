#ifndef EMFORGMM_H
#define EMFORGMM_H

#include "../randomalgorithm.h"

#include "../../base/commonutil.h"
#include "../../base/gmmutil.h"

/**
* @brief EM algorithm for GMMs
*/
class EMforGMM : public RandomAlgorithm
{
public:
	EMforGMM(commonutil::DataSet const& input, bool verbose, uint32_t seed);
	
	EMforGMM(commonutil::DataSet const& input, bool verbose, std::mt19937& gen);

	virtual ~EMforGMM()
	{
	}

	virtual void init(unsigned int);
	virtual void init(GMMDesc const& desc);
	virtual void run(unsigned int numSteps = 1);

	/**
	* @brief does a dynamic cast of the given GMMalgorithm to EMforGMM
	* @return NULL if the Algorithm is not a EMforGMM instance
	*/
	static EMforGMM* toEMforGMM(GMMAlgorithm* a);
};

#endif
