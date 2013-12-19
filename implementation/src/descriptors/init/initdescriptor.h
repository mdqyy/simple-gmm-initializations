#ifndef INITDESCRIPTOR_H
#define INITDESCRIPTOR_H

#include "../../base/gmmutil.h"

#include <vector>

/**
 * @brief Interface for describing data sets to be loaded or generated
 */
class InitDescriptor
{
public:
	virtual GMMDesc compute(commonutil::DataSet const&, unsigned int) = 0;
	std::string tag();
	
protected:
	std::string nametag;
};

#endif
