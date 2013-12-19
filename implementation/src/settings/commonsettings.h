#ifndef COMMONSETTINGS_H
#define COMMONSETTINGS_H

#include "../base/base.h"

#include <string>

struct CommonSettings
{
	bool verbose = false;
	CostMeasure costmeasure = NLL;
};

#endif
