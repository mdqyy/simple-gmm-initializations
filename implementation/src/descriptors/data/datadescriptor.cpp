#include "datadescriptor.h"

std::string DataDescriptor::tag()
{
	return this->nametag;
}

DataDescriptor::kIterator DataDescriptor::first_k() const
{
	return this->kList.begin();
}

DataDescriptor::kIterator DataDescriptor::last_k() const
{
	return this->kList.end();
}



