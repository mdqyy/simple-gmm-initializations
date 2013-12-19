#include "randomalgorithm.h"

RandomAlgorithm::RandomAlgorithm(commonutil::DataSet const& ds, bool v, uint32_t s)
	: GMMAlgorithm(ds,v), seed(s), gen(s)
{
	this->seed++;
}

RandomAlgorithm* RandomAlgorithm::toRandomAlgorithm(GMMAlgorithm* a)
{
	return dynamic_cast<RandomAlgorithm*>(a);
}

