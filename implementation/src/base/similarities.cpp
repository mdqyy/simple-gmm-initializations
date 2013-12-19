#include "similarities.h"

#include "gmmutil.h"


fp_type AverageLinkage::operator()(const Vector& sum1, idx_type count1, const Vector& sum2, idx_type count2)
{
	if(count1 + count2 <= 0)
		gmmlab_throw("AverageLinkage::operator - count1+count2 <= 0")
	return 1./(count1 + count2) * sum1.dot(sum2);
}

