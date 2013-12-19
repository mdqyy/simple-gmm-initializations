#ifndef DATADESCRIPTOR_H
#define DATADESCRIPTOR_H

#include "../../base/gmmutil.h"
#include "../../base/gmmdesc.h"

#include <vector>

/**
 * @brief Interface for describing data sets to be loaded or generated
 */
class DataDescriptor
{
public:
	typedef std::vector<unsigned int>::const_iterator kIterator;

	virtual GMMDesc retrieve(commonutil::DataSet&) = 0;
	std::string tag();
	kIterator first_k() const;
	kIterator last_k() const;
	
protected:
	std::string nametag;
	std::vector<unsigned int> kList;	
};

#endif
