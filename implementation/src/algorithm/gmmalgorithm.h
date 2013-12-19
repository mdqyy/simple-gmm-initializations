#ifndef GMMALGORITHM_H
#define GMMALGORITHM_H

#include "../base/gmmdesc.h"
#include "../base/commonutil.h"

/**
* @brief Abstract base class for GMM algorithms.
*/
class GMMAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	GMMAlgorithm(commonutil::DataSet const&, bool = false);

	GMMAlgorithm(const GMMAlgorithm&);
	GMMAlgorithm& operator=(const GMMAlgorithm&);

	virtual ~GMMAlgorithm()
	{
	}

	bool hasChanged() const;
	GMMDesc const& getGMMDesc() const;
	double getRuntime() const;

	virtual void setInput(commonutil::DataSet const&);
	virtual void init(GMMDesc const&);

	virtual void init(unsigned int) = 0;
	virtual void run(unsigned int numSteps = 1) = 0;

protected:
	commonutil::DataSet const* input;
	GMMDesc desc;
	double runtime;
	bool change;
	bool verbose;
};

#endif
