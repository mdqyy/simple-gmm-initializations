#include "descriptors.h"

#include "../settings/configparser.h"

#include "data/filedatadescriptor.h"
#include "data/gen_uniform_means_and_restricted_covar_dd.h"

#include "init/gmm/uniform_means2gmm_id.h"
#include "init/gmm/uniform_spherical_id.h"
#include "init/gmm/uniform_spherical_with_pruning_id.h"
#include "init/gmm/adaptive_means2gmm_id.h"
#include "init/gmm/adaptive_spherical_id.h"
#include "init/gmm/gonzalez2gmm_id.h"
#include "init/gmm/alternating_adaptivemeans_means2gmm_id.h"
#include "init/gmm/restart_init_id.h"
#include "init/gmm/agglomerative_init_id.h"
#include "init/gmm/gonzalez_forgmm_id.h"
#include "init/gmm/gonzalez_kwedlo_id.h"

#include "init/gmm/empty_id.h"

#include "algo/gmm/emalgodescriptor.h"

#include <fstream>
#include <boost/algorithm/string.hpp>


void ioutil::readRangeSpec(std::string const& spec, unsigned int& from, unsigned int& to,
	unsigned int& step)
{
	std::istringstream subiss(spec);
	std::string entry;
	
	getline(subiss, entry, '-');
	from = atoi(entry.c_str());
	entry.clear();
	
	getline(subiss, entry, '.');
	if (entry.empty())
		to = from;
	else
		to = atoi(entry.c_str());
	entry.clear();
	
	getline(subiss, entry);
	if (entry.empty())
		step = 1;
	else
		step = atoi(entry.c_str());
//	entry.clear();
}


std::vector<unsigned int> ioutil::readIntList(std::istringstream& iss)
{
	std::string list;
	getline(iss, list, ':');
	
	std::istringstream subiss(list);
	std::string entry;
	std::vector<unsigned int> v;
	
	while (getline(subiss, entry, ','))
	{
		unsigned int from, to, step;
		readRangeSpec(entry, from, to, step);
		for (unsigned int i=from; i<=to; i+=step)
			v.push_back(i);
	}
	
	return v;
}

std::vector<fp_type> ioutil::readFpTypeList(std::istringstream& iss)
{
	std::string list;
	getline(iss, list, ':');
	
	std::istringstream subiss(list);
	std::string entry;
	std::vector<fp_type> v;
	
	while (getline(subiss, entry, ','))
		v.push_back(atof(entry.c_str()));
	
	return v;
}

fp_type ioutil::readFpType(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	return atof(entry.c_str());
}

unsigned int ioutil::readInt(std::istringstream& iss)
{
	std::string entry;
	getline(iss,entry, ':');
	return atoi(entry.c_str());
}

std::string ioutil::readUpperString(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	boost::to_upper(entry);
	return entry;
}

std::string ioutil::readString(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	return entry;
}

std::vector<DataDescriptor*> ioutil::createDataDescriptors(std::string const& filename)
{
	std::ifstream filestream;
	std::vector<DataDescriptor*> descriptors;
	std::cout << "Creating data descriptors from file: " << filename << std::endl;

	filestream.open(filename.c_str());
	std::string line;
	while (getline(filestream, line))
	{
		std::istringstream iss(line);
		
		std::string dataType;
		dataType = readUpperString(iss);
		if(dataType.empty())
			continue;
		
		if (dataType=="LOAD")
		{
			// read file descriptor from strings formatted as:
			//		load:<filename>:<first coordinate>-<last coordinate>.<block length>:<k range specs>
			
			std::string prefix = readString(iss);

			std::string entry = readString(iss);
			unsigned int from, to, block;
			readRangeSpec(entry, from, to, block);
			unsigned int num = to<from?0:to-from+1;
			
			std::vector<unsigned int> kList = readIntList(iss);
			
			entry = readUpperString(iss);
			bool normalize = (entry=="NORMALIZE");
			
			FileDataDescriptor* fdd = new FileDataDescriptor(prefix, from, num, block, normalize);
			for (std::size_t i=0; i<kList.size(); ++i)
			{
				fdd->add_k(kList[i]);
			}
			
			descriptors.push_back(fdd);
		}
		else if (dataType=="GEN")
		{
			
			// read data generation descriptor from strings formatted as:
			//		gen:genMode: ...  
				
			
			std::string genMethod = readUpperString(iss);
			
			if(genMethod == "UMRC")
			{

				// uniform_means_and_restricted_covar:
				//  	 gen:umrc:<size>:<dim>:<k>
				//               :<Wexp fp-range-specs>
				//               :<separation fp-range-specs>
				//               :<EWexp fp-range-specs>:<EWmin fp>:<EWmax fp>:<EWPropMin fp>:<EWPropMax fp>
				//               :<ME string>
				//               :<seed range-specs>
				
				unsigned int size = readInt(iss);
				unsigned int d = readInt(iss);
				unsigned int k = readInt(iss);
				std::vector<fp_type> wList = readFpTypeList(iss);
				std::vector<fp_type> sepList = readFpTypeList(iss);
				std::vector<fp_type> ewExpList = readFpTypeList(iss);
				fp_type minMinSqrtEW = readFpType(iss);
				fp_type maxMinSqrtEW = readFpType(iss);
				fp_type minSqrtEWProp = readFpType(iss);
				fp_type maxSqrtEWProp = readFpType(iss);
				MeasurementError me;
				getMeasurementError(readString(iss), me);
				std::vector<unsigned int> gList = readIntList(iss);
				
				for(std::size_t g=0; g<gList.size(); ++g)
					for(std::size_t w=0; w<wList.size(); ++w)
						for(std::size_t sep=0; sep<sepList.size(); ++sep)
							for(std::size_t ewExp=0; ewExp<ewExpList.size(); ++ewExp)
								descriptors.push_back(
									new GenUniformMeansAndRestrictedCovarDD(
											size, d, k, wList[w], 
											sepList[sep],
											ewExpList[ewExp], minMinSqrtEW, maxMinSqrtEW, minSqrtEWProp, maxSqrtEWProp, 
											me, gList[g])
								);
				
			}
			else
			{
					gmmlab_throw("descriptors::createDataDescriptors() - Unknown generation method.");
			}
		}
	}
	filestream.close();
	
	return descriptors;
}

std::vector<InitDescriptor*> ioutil::createInitDescriptors(std::string const& filename)
{
	std::ifstream filestream;
	std::vector<InitDescriptor*> descriptors;
	std::cout << "Creating init descriptors from file: " << filename << std::endl;

	filestream.open(filename.c_str());
	std::string line;
	while(getline(filestream, line))
	{
		std::istringstream iss(line);
		
		InitMethod algo;
		std::string algoString = readString(iss);
		if(algoString.empty())
			continue;
		getInitMethod(algoString, algo);
		
		// read init descriptor from strings formatted as:
		//                                   name:<seed range specs> or
		// e.g.:
		// - uniform_spherical_with_pruning:<oversamplingfactor range spec>:<seed range specs>
		// - alternating_adaptivemeans_means2gmm:<factor range spec>:<seed range specs>
		
		std::vector<unsigned int> oList;
		std::vector<unsigned int> initSubrunsList;
		std::vector<unsigned int> initSubstepsList;
		std::vector<fp_type> initSampleSizeFactorList;
		std::vector<fp_type> initUniformFactorList;
		std::string initSubInitMethod;
		std::vector<unsigned int> initUse2GMMList;
		
		if(algo == UNIFORM_SPHERICAL_WITH_PRUNING)
		{
			oList = readIntList(iss);
		}
		else if(algo == RESTART_INIT)
		{
			initSubInitMethod = readUpperString(iss);
			initSubrunsList = readIntList(iss);
			initSubstepsList = readIntList(iss);
			InitMethod tmp;
			getInitMethod(initSubInitMethod, tmp);
			if(tmp == ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM)
				initUniformFactorList = readFpTypeList(iss);
		}
		else if(algo == SAMPLE_AVGLINK_MEANS2GMM)
		{
			initSampleSizeFactorList = readFpTypeList(iss);
		}
		else if(algo == ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM)
		{
			initUniformFactorList = readFpTypeList(iss);
		}
		else if(algo == GONZALEZ_FORGMM)
		{
			initUse2GMMList = readIntList(iss);
			initSampleSizeFactorList = readFpTypeList(iss);
		}else if(algo == GONZALEZ_KWEDLO)
		{
			initSampleSizeFactorList = readFpTypeList(iss);
		}
			
		std::vector<unsigned int> sList = readIntList(iss);
				
		for (idx_type s=0; s<sList.size(); ++s)
		{
		  
			switch(algo){
				
				// gmm
				case UNIFORM_SPHERICAL: 
					descriptors.push_back(new UniformSphericalID(sList[s]));
					break;
				case UNIFORM_SPHERICAL_WITH_PRUNING:
					for(idx_type ok = 0; ok<oList.size(); ++ok)
						descriptors.push_back(new UniformSphericalWithPruningID(oList[ok], sList[s]));
					break;
				case ADAPTIVE_SPHERICAL:
					descriptors.push_back(new AdaptiveSphericalID(sList[s]));
					break;
				case UNIFORM_MEANS2GMM:
					descriptors.push_back(new UniformMeans2GMMID(sList[s]));
					break;
				case ADAPTIVE_MEANS2GMM:
					descriptors.push_back(new AdaptiveMeans2GMMID(sList[s]));
					break;
				case GONZALEZ2GMM:
					descriptors.push_back(new Gonzalez2GMMID(sList[s]));
					break;
				case ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM:
					for(idx_type f = 0; f<initUniformFactorList.size(); ++f)
						descriptors.push_back(new AlternatingAdaptiveMeansAndMeans2GMMID(sList[s], initUniformFactorList[f]));
					break;
				case RESTART_INIT:
					for(idx_type sr = 0; sr<initSubrunsList.size(); ++sr)
						for(idx_type sst = 0; sst<initSubstepsList.size(); ++sst)
						{
							InitMethod tmp;
							getInitMethod(initSubInitMethod, tmp);
							if(tmp == ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM)
								for(idx_type f = 0; f< initUniformFactorList.size(); ++f)
									descriptors.push_back(new RestartInitID(initSubInitMethod, initSubrunsList[sr], initSubstepsList[sst], sList[s], initUniformFactorList[f]));
							else
								descriptors.push_back(new RestartInitID(initSubInitMethod, initSubrunsList[sr], initSubstepsList[sst], sList[s]));
						}
					break;
				case SAMPLE_AVGLINK_MEANS2GMM:
					for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
						descriptors.push_back(new AgglomerativeInitID(initSampleSizeFactorList[sf],sList[s]));
					break;
				case GONZALEZ_FORGMM:
					for(idx_type u = 0; u<initUse2GMMList.size(); ++u)
						for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
							descriptors.push_back(new GonzalezForGMMID(sList[s],(initUse2GMMList[u]>0), initSampleSizeFactorList[sf]));
					break;
				case GONZALEZ_KWEDLO:
					for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
						descriptors.push_back(new GonzalezKwedlo(sList[s],initSampleSizeFactorList[sf]));
					break;
				
				// empty
				case EMPTY:
					descriptors.push_back(new EmptyID());
					break;
					
				default:
					std::cout << "Failed to create an initialization named " << algoString << std::endl;
					gmmlab_throw("ioutil::createInitDescriptors() - Unknown initialization.");
					
			}
		}	
	}
	filestream.close();
	
	return descriptors;
}

std::vector<AlgoDescriptor*> ioutil::createAlgoDescriptors(std::string const& filename)
{
	std::ifstream filestream;
	std::vector<AlgoDescriptor*> descriptors;
	std::cout << "Creating algo descriptors from file: " << filename << std::endl;

	filestream.open(filename.c_str());
	std::string line;
	std::vector<AlgoDescriptor*> last;
	while(getline(filestream, line))
	{
		bool first = false;
		std::istringstream iss(line);
		std::string entry;
		
		entry = readUpperString(iss);
		boost::erase_all(entry, " ");		
		if (entry=="FIRST")
		{
			first = true;
			last.clear();
		}
		else if (entry=="NEXT")
		{
			if (last.empty())
				gmmlab_throw("ioutil::createAlgoDescriptors() - next:... before first:... forbidden!");
		}
		else if (entry.empty())
		{
			last.clear();
			continue; // skip comments
		}

		AlgorithmID algo;
		getID(readUpperString(iss), algo);
		
		std::vector<unsigned int> stList;
		std::vector<unsigned int> sdList;
		
		switch(algo){
			
			// randomized algorithms needing a number of steps and a seed
			case EM:
				
				// read algo descriptor from strings formatted as:
				//		first:em:<steps range spec>:<seed range spec>
				// or:
				//		next:em:<steps>:<seed>
				
				stList = readIntList(iss);	
				sdList = readIntList(iss);
				
				if (first)	
					for (std::size_t st=0; st<stList.size(); ++st)
						for (std::size_t sd=0; sd<sdList.size(); ++sd)
						{
							AlgoDescriptor* ad = new EMAlgoDescriptor(stList[st], sdList[sd]);
							descriptors.push_back(ad);
							last.push_back(ad);
						}
				else
					for (std::size_t i=0; i<last.size(); ++i)
					{
						if (stList.size()>1)
							gmmlab_throw("ioutil::createAlgoDescriptors() - multiple step values forbidden for next:...!");
						if (sdList.size()>1)
							gmmlab_throw("ioutil::createAlgoDescriptors() - multiple seed values forbidden for next:...!");

						unsigned int steps = stList[0];
						// steps==0 stands for copy the number of steps from last
						if (steps==0)
							steps = last[i]->steps();
							
						AlgoDescriptor* ad = new EMAlgoDescriptor(steps, sdList[0]);
						last[i]->setNext(ad);
						last[i] = ad;
					}
					
					break;
					
			default:
				gmmlab_throw("descriptors::createAlgoDescriptors() - Unknown algorithm.")
		
		} // switch
	}
	filestream.close();
	
	return descriptors;
}
