#include "ioutil.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdint>
#include <dirent.h>
#include <boost/algorithm/string.hpp>
//#include <boost/regex.hpp>

#include "base.h"
#include "gmmdesc.h"
#include "gmmutil.h"

#include "../settings/settings.h"


bool ioutil::endsWith(std::string s, std::string postfix)
{
	return s.size()>=postfix.size() && s.substr(s.size()-postfix.size(), std::string::npos) == postfix;
}

void ioutil::setPrecision(std::ofstream& filestream){
	filestream << std::setiosflags(std::ios::scientific)
		<< std::setiosflags(std::ios::showpoint)
		<< std::setprecision(std::numeric_limits<fp_type>::digits10);
}

void ioutil::openStatisticFile(std::string const& path)
{
	std::string name = path;
	name.append(CSV_FILE_EXTENSION);

	std::cout << "-- writing statistics to file " << name << std::endl;
	statistic_filestream.open(name.c_str());
	setPrecision(statistic_filestream);
	statistic_filestream << "round" << SEPARATOR << "runtime" << SEPARATOR << "cost" << std::endl;
	statistic_first_entry_done = false;
}


bool ioutil::StatisticFileIsOpen()
{
	return ioutil::statistic_filestream.is_open();
}


void ioutil::storeStatisticEntry(std::size_t round, double runtime, fp_type costs)
{	
	statistic_filestream << round << SEPARATOR << runtime << SEPARATOR << costs << std::endl;
}


void ioutil::storeStatisticMarker(std::size_t mainrnd, double maintime, std::size_t minrnd, double mintime)
{
	statistic_filestream << mainrnd << SEPARATOR << maintime << std::endl;
	statistic_filestream << minrnd << SEPARATOR << mintime << std::endl;
}

void ioutil::storeStatisticMarkerVector(std::vector<std::size_t> endRounds,
	std::vector<double> endTimes, std::size_t minrnd, double mintime)
{
	if (endRounds.size()!=endTimes.size())
		gmmlab_throw("ioutil::storeStatisticMarkerVector() - round and time vectors have different length!");
		
	for (std::size_t i=0; i<endRounds.size(); ++i)
		statistic_filestream << endRounds[i] << SEPARATOR << endTimes[i] << std::endl;
	statistic_filestream << minrnd << SEPARATOR << mintime << std::endl;
}


void ioutil::closeStatisticFile()
{
	statistic_filestream.close();
}


void ioutil::openDataFile(std::string const& path)
{
	std::string name = path;
	name.append(CSV_FILE_EXTENSION);

	std::cout << "-- writing data to file " << name << std::endl;
	data_filestream.open(name.c_str());
	setPrecision(statistic_filestream);
}


bool ioutil::DataFileIsOpen()
{
	return ioutil::data_filestream.is_open();
}


void ioutil::storeDataEntry(Vector entry)
{	
	idx_type d = entry.size();
	for(idx_type i=0; i<d; ++i){
		data_filestream << entry(i);
		if(i<d-1)
			data_filestream << ioutil::SEPARATOR;
	}
	data_filestream << std::endl;
}


void ioutil::closeDataFile()
{
	data_filestream.close();
}


void ioutil::closeAllFiles()
{
	if (ioutil::StatisticFileIsOpen())
		ioutil::closeStatisticFile();
	if(ioutil::DataFileIsOpen())
		ioutil::closeDataFile();
}


void ioutil::loadDataFromCSV(commonutil::DataSet& target, std::string const& filename,
	unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize)
{
	std::ifstream filestream;

	std::string line;
	idx_type num = 0;
	idx_type dim = 0;

	idx_type maxdim = numCoords==0?0:firstCoord+numCoords-1;
	// in the first pass count the number of lines and check, if all input points have the
	// same dimension or if maxdim > 0, if all have at least dimension maxdim
	if (commonSettings().verbose)
		std::cout << "analyzing file " << filename << std::endl;
	filestream.open(filename.c_str());
	while(getline(filestream, line))
	{
		num++;
		std::istringstream iss(line);
		idx_type tmpdim = 0;
		std::string entry;
		while (getline(iss, entry, ','))
		{
			++tmpdim;
			if (tmpdim==maxdim)
				break;
		}

		if(dim == 0)
			dim = tmpdim;
		if (dim != tmpdim)
		{
			std::cout << "ioutil::loadDataFromCSV - found data points with different dimensions "
					  << dim << " and " << tmpdim << std::endl;
			filestream.close();
			gmmlab_throw("ioutil::loadDataFromCSV - loading failed!!!");
		}
	}
	filestream.close();

	if (firstCoord>dim)
	{
		if (commonSettings().verbose)
			std::cout << "ioutil::loadDataFromCSV - selected coordinates out of range!" << std::endl;
		return;
	}

	// dimension
	unsigned int coordsToBeLoaded = numCoords;
	if (numCoords==0 || firstCoord+numCoords-1>dim)
		coordsToBeLoaded = dim-firstCoord+1;

	idx_type shiftlen = (blockLength-1)*coordsToBeLoaded;

	unsigned int targetDim = blockLength*coordsToBeLoaded;
	if (target.weights.size()==0)
		target.points.resize(targetDim, 0);
	if (target.points.rows()!=targetDim)
		gmmlab_throw("ioutil::loadDataFromCSV() - target dataset has incompatible dimension.");

	// resize the input matrix so that each data point will form one column
	unsigned int blocks = num-blockLength+1;
	idx_type oldsize = target.weights.size();
	idx_type newsize = target.weights.size()+blocks;
	target.weights.conservativeResize(newsize);
	target.weights.tail(blocks) = Vector::Ones(blocks);
	target.points.conservativeResize(targetDim, newsize);

	// in the second pass blockwise read the coordinates into the input matrix
	std::cout << "reading " << blocks << " data points in " << targetDim << " dimensions from "
			  << num << " lines of " << filename << std::endl;
	Vector oneblock(targetDim);
	filestream.open(filename.c_str());
	for (idx_type i=0; i<num; ++i)
	{
		getline(filestream, line);
		std::istringstream iss(line);
		std::string entry;
		oneblock.head(shiftlen) = oneblock.tail(shiftlen);
		for (idx_type j=0; j<dim; ++j)
		{
			getline(iss, entry, ',');
			if(j>=firstCoord-1 && j<firstCoord-1+coordsToBeLoaded)
				oneblock(targetDim-coordsToBeLoaded+j-firstCoord+1)=atof(entry.c_str());
		}
		if (i>=blockLength-1)
			target.points.col(oldsize+i) = oneblock;
//		std::cout << "new block: " << oneblock << std::endl;
	}
	filestream.close();

	// normalization
	if (normalize)
	{
		std::cout << "   normalizing data..." << std::endl;
		for (idx_type i=0; i<targetDim; ++i)
		{
			fp_type coeff = target.points.row(i).tail(blocks).minCoeff();
			target.points.row(i).tail(blocks) -= Vector::Constant(blocks, coeff);
			coeff = target.points.row(i).tail(blocks).maxCoeff();
			if (coeff>0)
				target.points.row(i).tail(blocks) /= coeff;
		}
	}
//	std::cout << std::endl;
}


void ioutil::getFilenames(std::string const& pattern, std::vector<std::string>& list)
{	
	// divide pattern into path and prefix
	std::string dirname, filepat;
	std::size_t index = pattern.find_last_of('/');
	if (index==std::string::npos)
		dirname = "./";
	else
		dirname = pattern.substr(0,index+1);
	filepat = pattern.substr(index+1,std::string::npos);

// 	dirent* thisentry;
// 	DIR *thisdir = opendir(".");
// 	std::cout << "------------------" << std::endl
// 	while(thisentry = readdir(thisdir))
// 		std::cout << thisentry->d_name << std::endl;
	
	// open directory
	if (commonSettings().verbose)
		std::cout << "searching directory " << dirname << std::endl;
	DIR* dir;
	dirent* entry;
	if ( (dir = opendir(dirname.c_str())) != NULL)
	{
		while ( (entry = readdir(dir)) != NULL)
			if (entry->d_name[0] != '.')
			{
				std::string filename(entry->d_name);
//				if (matchesWildcard(filename, filepat))
				if (filename.substr(0,filepat.size())==filepat)
				{
					if (commonSettings().verbose)
						std::cout << "   found " << filename << std::endl;
					std::string path(dirname);
					path += filename;
					list.push_back(path);
				}
			}
		closedir(dir);
		if (commonSettings().verbose)
			std::cout << list.size() << " files found." << std::endl;
	}
	else
	{
		std::cout << "ioutil::getFilenames() - Trying to open directory " << dirname << std::endl;
		gmmlab_throw("ioutil::getFilenames() - Failed to open directory.");
	}
}

const std::string ioutil::commonFileExtension(std::vector<std::string> const& list)
{
	if (list.empty())
		return "";

	std::size_t index = list.begin()->find_last_of('.');
	if (index==std::string::npos)
		return "";

	std::string extension = list.begin()->substr(index+1,std::string::npos);
	for (std::vector<std::string>::const_iterator it = list.begin(); it!=list.end(); ++it)
	{
		index = it->find_last_of('.');
		if (index==std::string::npos || it->substr(index+1,std::string::npos)!=extension)
			return "";
	}
	return extension;
}


void ioutil::loadData(commonutil::DataSet& target, std::string const& pattern,
	unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize)
{
	std::cout << "Loading datapoints from files starting with: " << pattern << std::endl;
	idx_type presize = target.weights.size();

	// get filenames from pattern and sort alphabetically
	std::vector<std::string> filelist;
	getFilenames(pattern, filelist);
	std::sort(filelist.begin(), filelist.end());

	if (filelist.empty())
		std::cout << "   WARNING - No files found!" << std::endl;
	else
	{
		// check for common file extension to apply proper load procedure
		std::string common = commonFileExtension(filelist);
		boost::to_upper(common);
		FileFormat ff = UNKNOWN;
		if (common == "CSV")
			ff = CSV;
		else
			std::cout << "   WARNING - Unknown or diverse file extensions!" << std::endl;

		// iterate over list and load files
		idx_type maxdrop = 0;
		for (std::vector<std::string>::iterator it = filelist.begin(); it!=filelist.end(); ++it)
			switch (ff)
			{
				case CSV:
					loadDataFromCSV(target, *it, firstCoord, numCoords, blockLength, normalize);
					break;
			}
		std::cout << "   " << target.weights.size()-presize << " data points read." << std::endl;
	}
}


