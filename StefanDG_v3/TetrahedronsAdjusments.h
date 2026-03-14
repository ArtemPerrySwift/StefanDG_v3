#pragma once

struct ElementsAdjusmentsSet
{
	double* bilinear;
	double* linear;

	size_t* elementsTags;
	size_t nElements;

	static size_t calcAllTrtrahedronsCount(const ElementsAdjusmentsSet* tetrahedronsAdjusmentsSetIt, const size_t nSets);
};


