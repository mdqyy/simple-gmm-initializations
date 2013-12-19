#ifndef SIMILARITIES_H
#define SIMILARITIES_H

#include "base.h"
#include "gmmdesc.h"


class AverageLinkage
{
public:
	fp_type operator()(Vector const& sum1, idx_type count1, Vector const& sum2, idx_type count2);
};

#endif