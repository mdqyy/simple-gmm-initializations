#include <iostream>
#include <iomanip>
#include <climits>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <map>
#include <boost/filesystem.hpp>

#include "base/base.h"
#include "base/similarities.h"
#include "base/ioutil.h"
#include "base/gmmutil.h"
#include "base/kmeansutil.h"

#include "settings/configparser.h"
#include "settings/settings.h"

#include "descriptors/descriptors.h"


GMMDesc truth;
commonutil::DataSet input;
std::vector<double> runtimes;
std::vector<fp_type> costs;

fp_type mincost;
double mincost_runtime = 0;
std::size_t mincost_round = 0;



inline double secondsSince(clock_t t)
{
	return double(clock()-t)/CLOCKS_PER_SEC;
}

void init() {

	std::cout << std::endl << std::endl
		<< "TestLab - Run automated tests of initialization methods" << std::endl
		<< "-------------------------------------------------------" << std::endl << std::endl;

	std::cout << std::setiosflags(std::ios::scientific)
		<< std::setiosflags(std::ios::showpoint)
		<< std::setprecision(std::numeric_limits<fp_type>::digits10);
}

fp_type getCosts(commonutil::DataSet const& data, GMMDesc const& gmmdesc)
{  		
	if(commonSettings().costmeasure == NLL)
		return gmmutil::nll(input, gmmdesc);
	gmmlab_throw("main_gmmlab::getCosts() - undefined costmeasure")
}

void handleSolution(std::size_t round, GMMDesc const& gmmdesc, double time, bool hasChanged)
{
	
	if (ioutil::StatisticFileIsOpen())
	{
		// compute costs (only if the solution has changed)
		if(!costs.empty() && !hasChanged)
		    costs.push_back(costs.back());
		else
		    costs.push_back(getCosts(input, gmmdesc));		
		
		// output
		ioutil::storeStatisticEntry(round, time, costs.back());

		// minimum?
		if(costs.size()==1 || (mincost > costs.back()))
		{
			mincost = costs.back();
			mincost_runtime = time;
			mincost_round = round;
		}
	}
}


double runAlgoSequence(AlgoDescriptor const* descriptor, GMMDesc initialSolution, double initialTime,
	std::string const& outputDir, std::string const& outputFile)
{
	// Initialize output
	std::stringstream pathStream;
	pathStream << outputDir << "/csv";
	boost::filesystem::path outPath(pathStream.str());
	boost::filesystem::create_directories(outPath);
	std::stringstream filenameStream;
	filenameStream << pathStream.str() << "/" << outputFile;
	ioutil::openStatisticFile(filenameStream.str());


	// compute maximum spread
	Vector spreads = (input.points.rowwise().maxCoeff()
					 -input.points.rowwise().minCoeff());
					 
	if (commonSettings().verbose)
		std::cout << "maximum spread is " << spreads.maxCoeff() << std::endl;

	// Initialisation
	double totaltime = initialTime;
	GMMDesc last = initialSolution;

	// handle initial solution
	handleSolution(0, last, totaltime, true);

	// start computation of new algorithm sequence
	std::cout << "starting new algorithm sequence:" << std::endl;
			  
	AlgoDescriptor const* nextAlgo = descriptor;
	std:size_t roundOffset = 0;
	std::vector<std::size_t> finalRounds;
	std::vector<double> finalTimes;
	while (nextAlgo!=NULL)
	{
		std::cout << "creating algorithm from descriptor: "
				  << nextAlgo->tag() << std::endl;
				  
		GMMAlgorithm* algorithm = nextAlgo->create(input);
		algorithm->init(last);
		
		// run current algorithm
		for (std::size_t i=0; i<nextAlgo->steps(); ++i)
		{
			algorithm->run();
			handleSolution(roundOffset+i+1, algorithm->getGMMDesc(),
				algorithm->getRuntime()+totaltime, algorithm->hasChanged());
		}
		
		// Save results and delete algorithm
		last = algorithm->getGMMDesc();
		roundOffset += nextAlgo->steps();
		totaltime += algorithm->getRuntime();
		finalRounds.push_back(roundOffset);
		finalTimes.push_back(totaltime);
		std::cout << "   computation took " << algorithm->getRuntime() << " seconds." << std::endl;

		delete algorithm;

		nextAlgo = nextAlgo->getNext();
	}
	std::cout << std::endl;

	if (ioutil::StatisticFileIsOpen())
	{
		// Store statistics markers
		ioutil::storeStatisticMarkerVector(finalRounds, finalTimes, mincost_round, mincost_runtime);

		// Print overall results
		if (!truth.empty())
			std::cout << "         cost of truth = " << getCosts(input, truth) << std::endl;
		std::cout << "cost of final solution = " << costs.back() << std::endl;
		std::cout << "     cheapest solution = " << mincost << " after " << mincost_runtime << " seconds" << std::endl;
	}

	std::cout << std::endl;
	
	// Finalize output
	ioutil::closeAllFiles();
	
	return totaltime-initialTime;
}

void runDefaultTests()
{
	std::vector<DataDescriptor*> dataDescriptors = ioutil::createDataDescriptors(testlabSettings().datafile);
	std::vector<InitDescriptor*> initDescriptors = ioutil::createInitDescriptors(testlabSettings().initfile);
	std::vector<AlgoDescriptor*> algoDescriptors = ioutil::createAlgoDescriptors(testlabSettings().algofile);

	// loop over all data sets
	for (unsigned int dataIndex=0; dataIndex<dataDescriptors.size(); ++dataIndex)
	{
		// clear input
		input.weights.resize(0);
		input.points.resize(0,0);
	
		std::cout << ">" << std::endl
				  << ">> retrieving new data set from descriptor: "
				  << dataDescriptors[dataIndex]->tag() << std::endl
				  << ">" << std::endl;

		// load or generate data from descriptor		
		truth = dataDescriptors[dataIndex]->retrieve(input);
		if (input.weights.size()==0)
		{
			std::cout << "   skipping empty data set." << std::endl << std::endl;
			continue;
		}
		
		idx_type n = input.points.cols();
		idx_type d = input.points.rows();
//		std::cout << std::endl;
		
		
		std::cout << "Using " << d << " dimensions of the " << n << " input points." << std::endl;
		
		// loop over all values of k
		for (DataDescriptor::kIterator kIter=dataDescriptors[dataIndex]->first_k();
				kIter!=dataDescriptors[dataIndex]->last_k(); ++kIter)
		{
			unsigned int k = *kIter;
			std::cout << std::endl << "> computing solutions with k=" << k << std::endl << std::endl;
										  
			std::stringstream outDirStream;
			outDirStream << testlabSettings().outputDir << "/"
						 << dataDescriptors[dataIndex]->tag()
						 << "_k" << k;
			
		    // for average time stats
			double totalAlgotime = 0;
			unsigned int totalAlgoruns = 0;
			
			// loop over all initial solutions
			for (unsigned int initIndex=0; initIndex<initDescriptors.size(); ++initIndex)
			{
				std::cout << "computing initial solution from descriptor: "
						  << initDescriptors[initIndex]->tag() << std::endl;
							  
				clock_t t = clock();
				GMMDesc initial = initDescriptors[initIndex]->compute(input, k);
				double inittime = secondsSince(t);
				std::cout << "   computed in " << inittime << " seconds." << std::endl << std::endl;
					
				for (int algoIndex=0; algoIndex<algoDescriptors.size(); ++algoIndex)
				{
					std::stringstream outFileStream;
					outFileStream << initDescriptors[initIndex]->tag() << "__";
					AlgoDescriptor* ad = algoDescriptors[algoIndex];
					while (ad!=NULL)
					{
						outFileStream << ad->tag();
						ad = ad->getNext();
						if (ad!=NULL)
							outFileStream << "_and_";
					}
					
					try
					{
						double algotime = runAlgoSequence(algoDescriptors[algoIndex], initial, inittime,
							outDirStream.str(), outFileStream.str());
						totalAlgoruns++;
						totalAlgotime += algotime;
					}
					catch (std::exception& e)
					{
						std::cout << "runDefaultTests() - EXCEPTION CAUGHT!" << std::endl;
						std::cout << e.what() << std::endl;
						std::cout << "runDefaultTests() - Continuing with next algorithm sequence." << std::endl;
						ioutil::closeAllFiles();
					}
				}
			}			
		}
	}
}

void runDataGeneration()
{
	std::vector<DataDescriptor*> dataDescriptors = ioutil::createDataDescriptors(testlabSettings().datafile);

	std::stringstream pathStream;
	pathStream << testlabSettings().outputDir << "/";
	boost::filesystem::path outPath(pathStream.str());
	boost::filesystem::create_directories(outPath);
	
	// loop over all data sets
	for (unsigned int dataIndex=0; dataIndex<dataDescriptors.size(); ++dataIndex)
	{
		// clear input
		input.weights.resize(0);
		input.points.resize(0,0);
	
		std::cout << ">" << std::endl
				  << ">> retrieving new data set from descriptor: "
				  << dataDescriptors[dataIndex]->tag() << std::endl
				  << ">" << std::endl;

		// load or generate data from descriptor		
		dataDescriptors[dataIndex]->retrieve(input);
		if (input.weights.size()==0)
		{
			std::cout << "   skipping empty data set." << std::endl << std::endl;
			continue;
		}
		
		idx_type n = input.points.cols();
		idx_type d = input.points.rows();
		
		ioutil::openDataFile(pathStream.str() + dataDescriptors[dataIndex]->tag());
		for(idx_type i=0; i<n; i++)
			ioutil::storeDataEntry(input.points.col(i));
		ioutil::closeDataFile();
	}
}

int main(int argc, char* argv[])
{

#ifndef NDEBUG
	std::cout << std::endl
		<< " /-----------------------------------------\\" << std::endl
		<< "(  DEBUG mode is on: asserts are processed  )" << std::endl
		<< " \\-----------------------------------------/" << std::endl;
#endif

	init();
	parseTestLabConfiguration(argc, argv);

	switch(testlabSettings().mode)
	{

	case DEFAULT_TEST:
		std::cout << "--> Running Default Tests:" << std::endl << std::endl;
		try{
			runDefaultTests();
		}
		catch (std::exception& e)
		{
			std::cout << "main() - EXCEPTION CAUGHT!" << std::endl;
			std::cout << e.what() << std::endl;
		}
		break;

	case GENERATE_CSV:
		
		std::cout << "--> Running Data Generation:" << std::endl << std::endl;
		try{
			runDataGeneration();
		}
		catch (std::exception& e)
		{
			std::cout << "main() - EXCEPTION CAUGHT!" << std::endl;
			std::cout << e.what() << std::endl;
		}
		break;

	}

	return 0;
}
