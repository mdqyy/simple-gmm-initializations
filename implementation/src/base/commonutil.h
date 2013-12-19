#ifndef COMMONUTIL_H
#define COMMONUTIL_H

#include "base.h"

#include <iostream>

namespace commonutil
{
	
	struct InitDescriptor
	{
		InitMethod method;
		uint32_t seed;
	};

	struct AlgoDescriptor
	{
		AlgorithmID algoId;
		unsigned int steps;
		uint32_t seed;
	};
		
	struct DataSet
	{
		Matrix points;
		Vector weights;
	};
	
	/**
	 * Erases the i-th element from a vector.
	 * @param v vector
	 * @param i the index
	 */
	void erase(Vector&, idx_type);
	
	/**
	 * Erases the i-th column from a matrix.
	 * @param m matrix
	 * @param i the index
	 */
	void eraseColumn(Matrix&, idx_type);
	
	/**
	 * fills the given vector by real values. The i-th entry is drawn uniformly at random from [min(i), max(i)].
	 */
	void fill(Vector& vector, std::mt19937& gen, Vector min, Vector max);
	
	/**
	 * fills the given vector by real values that are drawn uniformly at random from [min, max].
	 */
	void fill(Vector& vector, std::mt19937& gen, fp_type min, fp_type max);
	
	/**
	 * fills the given vector by real values that are drawn uniformly at random from [min, max].
	 */
	void fill(Matrix& matrix, std::mt19937& gen, fp_type min, fp_type max);
	
	/**
	 * returns true iff there does *not* exist an index i where p1(i)=p2(i). p1 and p2 must have same length.
	 */
	bool componentWiseDifferent(Vector, Vector);
	
	void appendData(DataSet& dataset, DataSet const& additionalData);


	template<typename T> const std::string toString(const T&);
	template<typename T> const T fromString(const std::string&);
	template<typename RndEngine> const idx_type randomIndex(Vector const&, RndEngine&);
	
}

std::ostream& operator<<(std::ostream&, commonutil::DataSet const&);

template<class T> const std::string commonutil::toString(const T& t)
{
     std::ostringstream stream;
     stream << t;
     return stream.str();
}

template<class T> const T commonutil::fromString(const std::string& s)
{
     std::istringstream stream (s);
     T t;
     stream >> t;
     return t;
}

template<typename RndEngine> const idx_type commonutil::randomIndex(Vector const& weights, RndEngine& re)
{
	idx_type size = weights.size();
	if (size==0)
		return 0;

	fp_type sum = weights.sum();
	assert(sum>=0);

	std::uniform_real_distribution<> urd(0, sum);
	fp_type r = urd(re);

	std::size_t index = 0;
	fp_type w = weights[0];
	while (r>w&&index<size-1)
	{
		r -= w;
		w = weights[++index];
	}
	return index;
}

#endif