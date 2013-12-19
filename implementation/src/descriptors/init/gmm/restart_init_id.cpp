#include "restart_init_id.h"

#include "../../../base/base.h"
#include "../../../base/ioutil.h"
#include "../../../base/initutil.h"
#include "../../../settings/configparser.h"
#include <boost/algorithm/string.hpp>

RestartInitID::RestartInitID(std::string initmethod, unsigned int subruns, unsigned int substeps, uint32_t s, fp_type factor) : seed(s), subruns(subruns), substeps(substeps), factor(factor)
{
	std::stringstream sstream;
	getInitMethod(initmethod, this->initmethod);
	boost::to_lower(initmethod);
	if(this->initmethod == ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM)
		sstream << "RestartInit_" << initmethod << "_f" << factor << "_subruns" << subruns << "_substeps" << substeps << "_i" << s;
	else
		sstream << "RestartInit_" << initmethod << "_subruns" << subruns << "_substeps" << substeps << "_i" << s;
	this->nametag = sstream.str();
}

GMMDesc RestartInitID::compute(commonutil::DataSet const& input, unsigned int k)
{
	std::mt19937 gen(seed);
	return initutil::restartInit(input, k, gen, this->initmethod, subruns, substeps, factor);
}
