#ifndef TESTLABSETTINGS_H
#define TESTLABSETTINGS_H

#include "../base/base.h"

#include <string>

struct TestLabSettings
{
	TestLabMode mode = DEFAULT_TEST;
	std::string outputDir, datafile, initfile, algofile; 
};

#endif
