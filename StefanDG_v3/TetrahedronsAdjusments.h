#pragma once

struct TetrahedronsAdjusments
{
	double* bilinear;
	double* linear;

	size_t* tetrahedronsTags;
	size_t nTetrahedrons;

	static size_t calcAllTrtrahedronsCount(const TetrahedronsAdjusments* tetrahedronsAdjusmentsSetIt, const size_t nSets);
};


