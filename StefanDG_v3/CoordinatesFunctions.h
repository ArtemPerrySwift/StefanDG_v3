#pragma once
#include "LocalCoordinates2D.h"
#include "LocalCoordinates3D.h"
#include "Coordinates.h"
#include "GeometryConstants.h"

namespace CoordinatesFunctions
{
	void computeLocalJacobianMatrix(const double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
		const double determinant,
		double* localJacobianMatrixElementIt);

	void computeLocalJacobianMatrixies(const double(*transpJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
		const size_t nJacobianMatrixies,
		const double* determinnatIt,
		double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT]);

	void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D* localGradientIt, const uint16_t nPoints, Coordinates* gradientIt);

	void translate0FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate1FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate2FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate3FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);

	void translate(const LocalCoordinates3D localPoints[],
				   const uint8_t pointsCount,
				   const Coordinates& initPoint,
				   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
				   Coordinates points[]);

	void translate(const LocalCoordinates2D localPoints[],
				   const uint8_t pointsCount,
				   const Coordinates& initPoint,
				   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
				   Coordinates points[]);

	void coomputeDirectionalDerivative(const Coordinates* gradientIt, const uint16_t nGradients, const Coordinates direction, double* directionalDerivativeIt);

	void computeNormal(const Coordinates trianglePoints[constants::triangle::N_NODES], Coordinates& normal);
};

