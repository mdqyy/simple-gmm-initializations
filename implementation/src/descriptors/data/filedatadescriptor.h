#ifndef FILEDATADESCRIPTOR_H
#define FILEDATADESCRIPTOR_H

#include "datadescriptor.h"

#include <vector>

/**
 * @brief Descriptor for data sets to be loaded from file
 */
class FileDataDescriptor : public DataDescriptor
{
public:
	FileDataDescriptor(std::string const&, unsigned int, unsigned int, unsigned int, bool);

	virtual ~FileDataDescriptor()
	{
	}
	
	void add_k(unsigned int);
	
	GMMDesc retrieve(commonutil::DataSet& input);

private:
	std::string prefix;
	unsigned int firstCoord = 0;
	unsigned int numCoords = 0;
	unsigned int blockLength = 1;
	bool normalize = false;
};

#endif
