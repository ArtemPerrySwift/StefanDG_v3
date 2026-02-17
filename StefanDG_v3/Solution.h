#pragma once
#include <unordered_map>
#include "LocalCoordinates3D.h"
class DGStefanTask;

class Solution
{
public:
	friend class DGStefanTask;
	//friend Solution&& solveEllipticTask(const Material materials[], uint8_t nMaterials);

	Solution();
	Solution(const Solution& solution);
	Solution(Solution&& solution);

	~Solution();

	double compute(const LocalCoordinates3D& localPoint, const size_t tetrahedronTag) const;
	void computeLocalGradient(const LocalCoordinates3D& localPoint, const size_t tetrahedronTag, LocalCoordinates3D& localGrad) const;
	void compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, const size_t* tetrahedronTagIt, double* valueIt) const;

private:
	std::unordered_map<size_t, const double*> _DOFsByTetrahedronTag;
	double* _DOFs;
	size_t _nDOFs;
};


