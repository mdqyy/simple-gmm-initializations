#include "empty_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/gmmdesc.h"
#include "../../../base/initutil.h"

EmptyID::EmptyID()
{
	std::stringstream sstream;
	sstream << "Empty";
	this->nametag = sstream.str();
}

GMMDesc EmptyID::compute(commonutil::DataSet const& input, unsigned int k)
{
	GMMDesc gmmdesc;
	gmmdesc.weights.resize(0);
	gmmdesc.means.resize(0,0);
	gmmdesc.covariances.resize(0);	
	return gmmdesc;
}
