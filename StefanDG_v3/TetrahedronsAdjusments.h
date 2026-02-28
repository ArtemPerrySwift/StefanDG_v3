#pragma once

struct ElementsAdjusmentsSet
{
	double* bilinear;
	double* linear;

	size_t* tetrahedronsTags;
	size_t nTetrahedrons;

	static size_t calcAllTrtrahedronsCount(const ElementsAdjusmentsSet* tetrahedronsAdjusmentsSetIt, const size_t nSets);
};


