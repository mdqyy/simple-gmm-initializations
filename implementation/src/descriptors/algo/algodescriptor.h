#ifndef ALGODESCRIPTOR_H
#define ALGODESCRIPTOR_H

#include "../../algorithm/gmmalgorithm.h"

#include <vector>

/**
 * @brief Interface for describing algorithms to be run
 */
class AlgoDescriptor
{
public:
	virtual GMMAlgorithm* create(commonutil::DataSet const&) const = 0;
	std::string tag() const;

	void setNext(AlgoDescriptor*);
	AlgoDescriptor* getNext() const;
	
	unsigned int steps() const;
	
protected:
	std::string nametag;
	AlgoDescriptor* next = NULL;
	unsigned int numSteps = 1;
};

#endif
