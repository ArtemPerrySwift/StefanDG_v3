#include "LinearLagrangeBasis.h"


double LinearLagrangeBasis::compute(const LocalCoordinates3D& localPoint, const double* dofs)
{
	double startDOF = dofs[0];
	double sum = startDOF;
	sum += localPoint.u * (dofs[1] - startDOF);
	sum += localPoint.v * (dofs[2] - startDOF);
	sum += localPoint.w * (dofs[3] - startDOF);
	return sum;
}

void LinearLagrangeBasis::compute(const LocalCoordinates3D& localPoint, const double* dofs, LocalCoordinates3D& localGradient)
{
	localGradient.u = -dofs[0] + dofs[1];
	localGradient.v = -dofs[0] + dofs[2];
	localGradient.w = -dofs[0] + dofs[3];
}

void LinearLagrangeBasis::compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	valueIt = compute0Function(localPointIt, nPoints, valueIt);
	valueIt = compute1Function(localPointIt, nPoints, valueIt);
	valueIt = compute2Function(localPointIt, nPoints, valueIt);
	valueIt = compute3Function(localPointIt, nPoints, valueIt);
}

void LinearLagrangeBasis::compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	gradientIt = compute0Function(localPointIt, nPoints, gradientIt);
	gradientIt = compute1Function(localPointIt, nPoints, gradientIt);
	gradientIt = compute2Function(localPointIt, nPoints, gradientIt);
	gradientIt = compute3Function(localPointIt, nPoints, gradientIt);
}

double* LinearLagrangeBasis::compute0Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 1.0 - localPointIt->u - localPointIt->v - localPointIt->w;

		++valueIt;
		++localPointIt;
	}

	return valueIt;
}

double* LinearLagrangeBasis::compute1Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = localPointIt->u;

		++valueIt;
		++localPointIt;
	}

	return valueIt;
}

double* LinearLagrangeBasis::compute2Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = localPointIt->v;

		++valueIt;
		++localPointIt;
	}

	return valueIt;
}

double* LinearLagrangeBasis::compute3Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = localPointIt->w;

		++valueIt;
		++localPointIt;
	}

	return valueIt;
}

LocalCoordinates3D* LinearLagrangeBasis::compute0Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1;
		gradientIt->v = -1;
		gradientIt->w = -1;

		++gradientIt;
		++localPointIt;
	}

	return gradientIt;
}

LocalCoordinates3D* LinearLagrangeBasis::compute1Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = 1;
		gradientIt->v = 0;
		gradientIt->w = 0;

		++gradientIt;
		++localPointIt;
	}

	return gradientIt;
}

LocalCoordinates3D* LinearLagrangeBasis::compute2Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = 0;
		gradientIt->v = 1;
		gradientIt->w = 0;

		++gradientIt;
		++localPointIt;
	}

	return gradientIt;
}

LocalCoordinates3D* LinearLagrangeBasis::compute3Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = 0;
		gradientIt->v = 0;
		gradientIt->w = 1;

		++gradientIt;
		++localPointIt;
	}

	return gradientIt;
}


