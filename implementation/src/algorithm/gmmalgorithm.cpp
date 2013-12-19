#include "gmmalgorithm.h"

#include <algorithm>
#include <ctime>
#include <random>


GMMAlgorithm::GMMAlgorithm(commonutil::DataSet const& data, bool v) : input(&data), runtime(0),
	change(false), verbose(v)
{
	assert(data.points.cols()==data.weights.size());
}

GMMAlgorithm::GMMAlgorithm(const GMMAlgorithm& rhs) : input(rhs.input), desc(rhs.desc),
	runtime(rhs.runtime), change(rhs.change), verbose(rhs.verbose)
{	
}

GMMAlgorithm& GMMAlgorithm::operator=(const GMMAlgorithm& rhs)
{
	this->input = rhs.input;
	this->desc = rhs.desc;
	this->runtime = rhs.runtime;
	this->change = rhs.change;
	return *this;
}

void GMMAlgorithm::setInput(commonutil::DataSet const& data)
{
	this->input = &data;
}

void GMMAlgorithm::init(GMMDesc const& dsc)
{
	this->desc = dsc;
	this->change = true;
}

bool GMMAlgorithm::hasChanged() const
{
	return this->change;
}

GMMDesc const& GMMAlgorithm::getGMMDesc() const
{
	return this->desc;
}

double GMMAlgorithm::getRuntime() const
{
	return this->runtime;
}
