#include "gen_uniform_means_and_restricted_covar_dd.h"

#include "../../base/datagenutil.h"
#include "../../settings/configparser.h"

GenUniformMeansAndRestrictedCovarDD::GenUniformMeansAndRestrictedCovarDD(unsigned int n, unsigned int d, unsigned int k, fp_type weightFloatExp, fp_type separation, fp_type sqrtEWFloatExp, fp_type minMinSqrtEW, fp_type maxMinSqrtEW, fp_type minSqrtEWProp, fp_type maxSqrtEWProp, MeasurementError measurementError, uint32_t seed) :
size(n), dim(d), comp(k), weightFloatExp(weightFloatExp), separation(separation), sqrtEWFloatExp(sqrtEWFloatExp), minMinSqrtEW(minMinSqrtEW), maxMinSqrtEW(maxMinSqrtEW), minSqrtEWProp(minSqrtEWProp), maxSqrtEWProp(maxSqrtEWProp), measurementError(measurementError), genSeed(seed)
{
	std::stringstream sstream;
	sstream << "gen_umrc_n" << n << "_d" << d << "_k" << k  << "_sep" << separation << "_w" << weightFloatExp << "_ewexp" << sqrtEWFloatExp << "_ewmin" << minMinSqrtEW << "_ewmax" << maxMinSqrtEW << "_ewpmin" << minSqrtEWProp << "_ewpmax" << maxSqrtEWProp << "_me" << getName(measurementError) << "_g" << seed;
	this->nametag = sstream.str();
	
	this->kList.push_back(k);
}

GMMDesc GenUniformMeansAndRestrictedCovarDD::retrieve(commonutil::DataSet& input)
{
	std::mt19937 gen(this->genSeed);
	
	GMMDesc truth = datagenutil::generateRandomGMMWithUniformMeans(dim, comp, gen, separation,
													   minMinSqrtEW, maxMinSqrtEW, 
													   minSqrtEWProp, maxSqrtEWProp,
													   weightFloatExp, sqrtEWFloatExp);
	unsigned int uniformNoisePoints = round(0.1 * this->size);
		
	switch(this->measurementError)
	{
		case NO_ME: 
			datagenutil::generateInputFromGMM(input, truth, this->size, gen);
			break;
		case UNIFORM_NOISE:
			datagenutil::generateInputFromGMM(input, truth, this->size - uniformNoisePoints, gen);
			datagenutil::addUniformNoise(input, uniformNoisePoints, gen);
			break;
		default:
			gmmlab_throw("GenUniformMeansAndRestrictedCovarDD::retrieve() - Unknown type of measurement error");
	}
			
	return truth;
}



