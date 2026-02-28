#pragma once

struct CrossElementsAdjusmentsSet
{
	const double* bilinear;

	const size_t* elementIndexes;
	size_t nCrossElements;

	static size_t calcAllCrossElementsCount(const CrossElementsAdjusmentsSet* crossTetrahedronsAdjusmentsSetIt, const size_t nSets);
};

