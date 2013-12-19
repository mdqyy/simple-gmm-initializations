#ifndef IOUTIL_H
#define IOUTIL_H

#include "gmmutil.h"

#include <fstream>

namespace ioutil
{
	enum FileFormat { UNKNOWN, CSV };

	static std::ofstream data_filestream;
	static std::ofstream statistic_filestream;
	static fp_type statistic_last_cost; 
	static bool statistic_first_entry_done;

	const std::string CSV_FILE_EXTENSION = ".csv";

	const char SEPARATOR = ',';


	bool endsWith(std::string , std::string);
	void setPrecision(std::ofstream&);
	std::vector<std::string> readDirectory(const std::string&);

	void openStatisticFile(std::string const&);
	bool StatisticFileIsOpen();
	void storeStatisticEntry(std::size_t, double, fp_type);
	void storeStatisticMarker(std::size_t, double, std::size_t, double);
	void storeStatisticMarkerVector(std::vector<std::size_t>, std::vector<double>, std::size_t, double);
	void closeStatisticFile();
	
	void openDataFile(std::string const&);
	bool DataFileIsOpen();
	void storeDataEntry(Vector entry);
	void closeDataFile();

	void closeAllFiles();
	
	void loadDataFromCSV(commonutil::DataSet&, std::string const&,
		unsigned int = 1, unsigned int = 0, unsigned int = 1, bool = false);

	void loadData(commonutil::DataSet& target, std::string const& pattern,
		unsigned int firstCoord= 1, unsigned int numCoords = 0, unsigned int blockLength = 1,
		bool = false);

	void getFilenames(std::string const&, std::vector<std::string>&);
	const std::string commonFileExtension(std::vector<std::string> const&);
}

#endif