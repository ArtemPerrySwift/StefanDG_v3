#pragma once

struct CrossTetrahedronsAdjusments
{
	double* bilinear;

	size_t(*tetrahedronsIndexes)[2];
	size_t nCrossElements;

	static size_t calcAllCrossElementsCount(const CrossTetrahedronsAdjusments* crossTetrahedronsAdjusmentsSetIt, const size_t nSets);
};

