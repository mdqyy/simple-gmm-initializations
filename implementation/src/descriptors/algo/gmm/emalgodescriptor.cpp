#include "emalgodescriptor.h"

#include "../../../settings/settings.h"
#include "../../../algorithm/gmm/emforgmm.h"


EMAlgoDescriptor::EMAlgoDescriptor(unsigned int st, uint32_t sd) : seed(sd)
{
	this->numSteps = st;
	std::stringstream sstream;
	sstream << "EM" << "_s" << st << "_a" << sd;
	this->nametag = sstream.str();
}

GMMAlgorithm* EMAlgoDescriptor::create(commonutil::DataSet const& input) const
{
	return new EMforGMM(input, commonSettings().verbose, this->seed);
}
