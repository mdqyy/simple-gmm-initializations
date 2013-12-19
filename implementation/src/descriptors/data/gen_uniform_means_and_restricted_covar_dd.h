#ifndef GEN_UNIFORM_MEANS_AND_RESTRICTED_COVAR_DD_H
#define GEN_UNIFORM_MEANS_AND_RESTRICTED_COVAR_DD_H

#include "datadescriptor.h"

#include <vector>
#include <boost/concept_check.hpp>

/**
 * 
 */
class GenUniformMeansAndRestrictedCovarDD: public DataDescriptor
{
public:
	GenUniformMeansAndRestrictedCovarDD(unsigned int n, unsigned int d, unsigned int k, fp_type weightFloatExp, 
										fp_type separation,
										fp_type sqrtEWFloatExp, fp_type minMinSqrtEW, fp_type maxMinSqrtEW, fp_type minSqrtEWProp, fp_type maxSqrtEWProp, 
										MeasurementError measurementError,   
										uint32_t seed);

	virtual ~GenUniformMeansAndRestrictedCovarDD()
	{
	}
	
	GMMDesc retrieve(commonutil::DataSet& input);

private:
	unsigned int size, dim, comp;
	fp_type weightFloatExp;
	fp_type separation;
	fp_type sqrtEWFloatExp, minMinSqrtEW, maxMinSqrtEW, minSqrtEWProp, maxSqrtEWProp;
	MeasurementError measurementError;
	
	uint32_t genSeed;
};

#endif

