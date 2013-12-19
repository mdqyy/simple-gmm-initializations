#include "configparser.h"
#include "settings.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>

namespace po = boost::program_options;

std::string getName(AlgorithmID algorithm)
{
    switch(algorithm){
		case EM: 
			return "EM";
		default: 
			gmmlab_throw("configparser::getName(algoid) - Unknown algorithm id.");
	}
}

std::string getName(MeasurementError me)
{
	switch(me){
		case NO_ME: 
			return "0";
		case UNIFORM_NOISE:
			return "UN";
		default:
			gmmlab_throw("configparser::getName(me) - Unknown measurement error.");
	}
}

void getID(std::string algorithm, AlgorithmID& target)
{ 
	boost::to_upper(algorithm);
	
	// gmm
	if (algorithm == "EM")
		target = EM;
	else if(algorithm == "")
		return;
	else
	{
		std::cout << "Failed to identify " << algorithm << std::endl;
		gmmlab_throw("configparser::getID() - Unknown algorithm.");
	}
}

void getInitMethod(std::string initmethod, InitMethod& target)
{
	boost::to_upper(initmethod);
    
	if (initmethod=="EMPTY")
		target = EMPTY;
	else if (initmethod=="UNIFORM_SPHERICAL")
		target = UNIFORM_SPHERICAL;
	else if(initmethod == "UNIFORM_SPHERICAL_WITH_PRUNING")
		target = UNIFORM_SPHERICAL_WITH_PRUNING;
	else if (initmethod=="EXPERIMENTAL_ADAPTIVE")
		target = EXPERIMENTAL_ADAPTIVE;
	else if (initmethod=="ADAPTIVE_SPHERICAL")
		target = ADAPTIVE_SPHERICAL;
	else if (initmethod=="UNIFORM_MEANS2GMM")
		target = UNIFORM_MEANS2GMM;
	else if(initmethod =="ADAPTIVE_MEANS2GMM")
		target = ADAPTIVE_MEANS2GMM;
	else if(initmethod == "ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM")
		target = ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM;
	else if(initmethod == "GONZALEZ2GMM")
		target = GONZALEZ2GMM;
	else if(initmethod == "RESTART_INIT")
		target = RESTART_INIT;
	else if(initmethod == "SAMPLE_AVGLINK_MEANS2GMM")
		target = SAMPLE_AVGLINK_MEANS2GMM;
	else if(initmethod == "GONZALEZ_FORGMM")
		target = GONZALEZ_FORGMM;
	else if(initmethod == "GONZALEZ_KWEDLO")
		target = GONZALEZ_KWEDLO;
	else if(initmethod == "")
		return;
	else
	{
		std::cout << "Failed to identify " << initmethod << std::endl;
		gmmlab_throw("configparser::getInitMethod() - Unknown method.");
	}
}

void getGenMethod(std::string genmethod, GenMethod& target)
{
    boost::to_upper(genmethod);
	if(genmethod == "UNIFORM_MEANS_AND_RESTRICTED_COVAR")
		target =  UNIFORM_MEANS_AND_RESTRICTED_COVAR;
	else if(genmethod == "")
		return;
	else
	{
		std::cout << "Failed to identify " << genmethod << std::endl;
		gmmlab_throw("configparser::genGenMethod() - Unknown method.");
	}
}

void getAntiPropMode(std::string apmode, AntiPropMode& target)
{
	boost::to_upper(apmode);
	if (apmode=="INVERSEDENSITY")
		target =  INVERSEDENSITY;
	else if (apmode=="MAHALANOBIS")
		target = MAHALANOBIS;
	else if(apmode=="")
		return;
	else
	{
		std::cout << "Failed to identify " << apmode << std::endl;
		gmmlab_throw("configparser::getAntiPropMode() - Unknown mode.")
	}
}

void getPropMode(std::string pmode, PropMode& target)
{	
	boost::to_upper(pmode);
	if (pmode=="DENSITY")
		target =  DENSITY;
	else if (pmode=="NEIGHBORHOOD")
		target =  NEIGHBORHOOD;
	else if(pmode == "")
		return;
	else
	{
		std::cout << "Failed to identify " << pmode << std::endl;
		gmmlab_throw("configparser::getPropMode() - Unknown mode.")
	}
}

void getCostMeasure(std::string costmeasure, CostMeasure& target)
{
    boost::to_upper(costmeasure);
	
	if(costmeasure=="NLL")
	  target = NLL;
	else
	{
		std::cout << "No cost measure defined. NLL will be used." << std::endl;
		target = NLL;
	}
}

void getMeasurementError(std::string me, MeasurementError& target)
{
	boost::to_upper(me);
	if(me=="NO_ME" || me=="")
		target = NO_ME;
	else if(me=="UN")
		target = UNIFORM_NOISE;
	else if(me== "")
		return;
	else
	{
		std::cout << "Failed to identify " << me << std::endl;
		gmmlab_throw("configparser::getMeasurementError() - Unknown mode");
	}
}

void parseTestLabConfiguration(int argc, char* argv[])
{
    std::cout << "command line: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    std::string config_file, mode;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
    ("version", "print version string")
    ("help,h", "produce help message")
    ("config", po::value<std::string>(&config_file)->default_value("testlab.cfg"),
     "name of configuration file.")
    ;
    
    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
    ("mode,m", po::value<std::string>(&mode), "mode of operation")
    ("data", po::value<std::string>(&testlabSettings().datafile), "data descriptor file")
    ("init", po::value<std::string>(&testlabSettings().initfile), "init descriptor file")
    ("algo", po::value<std::string>(&testlabSettings().algofile), "algo descriptor file")
    ("output,o", po::value<std::string>(&testlabSettings().outputDir), "output directory")
    ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("data", 1).add("init", 1).add("algo", 1).add("diff", 1);

    try
    {
        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        std::ifstream ifs(config_file.c_str());

        if (ifs)
        {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }

        if (vm.count("help"))
        {
            std::cout << visible << "\n";
            exit(0);
        }

        if (vm.count("version"))
        {
            std::cout << "Test Laboratory, version 0.1\n";
            exit(0);
        }

        boost::to_upper(mode);
        if (mode == "DEFAULT")
            testlabSettings().mode = DEFAULT_TEST;
        else if (mode == "GENERATE_CSV")
		testlabSettings().mode = GENERATE_CSV;

        if (commonSettings().verbose)
        {
            std::cout << "Verbose output is ON." << std::endl << std::endl;
        }
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}