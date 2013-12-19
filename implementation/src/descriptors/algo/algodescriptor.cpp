#include "algodescriptor.h"

std::string AlgoDescriptor::tag() const
{
	return this->nametag;
}

void AlgoDescriptor::setNext(AlgoDescriptor* ad)
{
	this->next = ad;
}

AlgoDescriptor* AlgoDescriptor::getNext() const
{
	return this->next;
}
	
unsigned int AlgoDescriptor::steps() const
{
	return this->numSteps;
}
	