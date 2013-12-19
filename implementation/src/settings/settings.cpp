#include "settings.h"

CommonSettings& commonSettings()
{
	static CommonSettings cs;
	return cs;
}

TestLabSettings& testlabSettings()
{
	static TestLabSettings ts;
	return ts;
}

