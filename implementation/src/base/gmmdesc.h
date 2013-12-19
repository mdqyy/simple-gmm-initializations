#ifndef GMMDESC_H
#define GMMDESC_H

#include "base.h"

#include <vector>

/**
 * @brief Data structure for the parameters of a gaussian mixture.
 */
struct GMMDesc
{
	Vector weights;
	Matrix means;
	std::vector<Matrix> covariances;
	
	bool operator==(GMMDesc const& rhs) const;

	bool empty();
};

std::ostream& operator<<(std::ostream&, GMMDesc const&);

#endif
