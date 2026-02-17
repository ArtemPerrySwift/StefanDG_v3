#include "Solution.h"
#include "LinearLagrangeBasis.h"


Solution::Solution() : _DOFs{ nullptr }, _nDOFs{ 0 }
{

}

Solution::Solution(const Solution& solution) : _DOFsByTetrahedronTag{ solution._DOFsByTetrahedronTag }, _nDOFs{ solution._nDOFs }
{
	_DOFs = new double[_nDOFs];
	std::copy_n(solution._DOFs, _nDOFs, _DOFs);
}

Solution::Solution(Solution&& solution)
{
	_DOFsByTetrahedronTag = std::move(solution._DOFsByTetrahedronTag);
	_DOFs = solution._DOFs;
	_nDOFs = solution._nDOFs;

	solution._DOFs = nullptr;
	solution._nDOFs = 0;
}

Solution::~Solution()
{
	_DOFsByTetrahedronTag.clear();
	delete[] _DOFs;
}

double Solution::compute(const LocalCoordinates3D& localPoint, const size_t tetrahedronTag) const
{
	return LinearLagrangeBasis::compute(localPoint, _DOFsByTetrahedronTag.at(tetrahedronTag));
}

void Solution::computeLocalGradient(const LocalCoordinates3D& localPoint, const size_t tetrahedronTag, LocalCoordinates3D& localGrad) const
{
	LinearLagrangeBasis::compute(localPoint, _DOFsByTetrahedronTag.at(tetrahedronTag), localGrad);
}

void Solution::compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, const size_t* tetrahedronTagIt, double* valueIt) const
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = LinearLagrangeBasis::compute(*localPointIt, _DOFsByTetrahedronTag.at(*tetrahedronTagIt));

		++localPointIt;
		++tetrahedronTagIt;
		++valueIt;
	}
}