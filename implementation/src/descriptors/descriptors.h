#ifndef DESCRIPTORS_H
#define DESCRIPTORS_H

#include <string>
#include <vector>
#include <boost/concept_check.hpp>

#include "data/datadescriptor.h"
#include "init/initdescriptor.h"
#include "algo/algodescriptor.h"

namespace ioutil
{
	/**
	 * converts strings of the form a-b.c to three unsigned ints a, b and c, where b and c are optional
	 */
	void readRangeSpec(std::string const& spec, unsigned int& from, unsigned int& to, unsigned int& step);
	
	/**
	 * converts strings of the form <range spec>,<range spec>,... to a std::vector<unsigned int>
	 */
	std::vector<unsigned int> readIntList(std::istringstream& iss);
	
	/**
	 * converts strings of the form fp_type,fp_type, ... to a std::vector<fp_type>
	 */
	std::vector<fp_type> readFpTypeList(std::istringstream& iss);
	
	fp_type readFpType(std::istringstream& iss);
	std::string readString(std::istringstream& iss);
	std::string readUpperString(std::istringstream& iss);
	unsigned int readInt(std::istringstream& iss);

	std::vector<DataDescriptor*> createDataDescriptors(std::string const&);
	std::vector<InitDescriptor*> createInitDescriptors(std::string const&);
	std::vector<AlgoDescriptor*> createAlgoDescriptors(std::string const&);
}

#endif