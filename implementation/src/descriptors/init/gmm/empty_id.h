#ifndef EMPTY_ID_H
#define EMPTY_ID_H

#include "../../../base/gmmutil.h"
#include "../../../base/gmmdesc.h"
#include "../initdescriptor.h"

#include <vector>

/**
 * @brief Descriptor for Initialization
 */
class EmptyID : public InitDescriptor
{
public:
	EmptyID();

	virtual ~EmptyID()
	{
	}
	
	GMMDesc compute(commonutil::DataSet const&, unsigned int);
};

#endif
