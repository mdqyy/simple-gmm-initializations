#ifndef RESTART_INIT_ID_H
#define RESTART_INIT_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class RestartInitID : public InitDescriptor
{
public:
	RestartInitID(std::string initmethod, unsigned int subruns, unsigned int substeps, uint32_t seed, fp_type factor = 1.);

	virtual ~RestartInitID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);

private:
	uint32_t seed;
	InitMethod initmethod;
	unsigned int subruns;
	unsigned int substeps;
	fp_type factor;
};

#endif
