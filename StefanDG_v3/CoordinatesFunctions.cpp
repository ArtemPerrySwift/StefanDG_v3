#include "CoordinatesFunctions.h"
#include <cmath>

namespace CoordinatesFunctions
{
    void computeLocalJacobianMatrix(const double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT], const double determinant, double* localJacobianMatrixElementIt)
    {
        *localJacobianMatrixElementIt = (transpJacobianMatrix[4] * transpJacobianMatrix[8] - transpJacobianMatrix[5] * transpJacobianMatrix[7]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[5] * transpJacobianMatrix[6] - transpJacobianMatrix[3] * transpJacobianMatrix[8]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[3] * transpJacobianMatrix[7] - transpJacobianMatrix[4] * transpJacobianMatrix[6]) / determinant;
        ++localJacobianMatrixElementIt;
        
        *localJacobianMatrixElementIt = (transpJacobianMatrix[2] * transpJacobianMatrix[7] - transpJacobianMatrix[1] * transpJacobianMatrix[8]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[0] * transpJacobianMatrix[8] - transpJacobianMatrix[2] * transpJacobianMatrix[6]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[1] * transpJacobianMatrix[6] - transpJacobianMatrix[0] * transpJacobianMatrix[7]) / determinant;
        ++localJacobianMatrixElementIt;
        
        *localJacobianMatrixElementIt = (transpJacobianMatrix[1] * transpJacobianMatrix[5] - transpJacobianMatrix[2] * transpJacobianMatrix[4]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[2] * transpJacobianMatrix[3] - transpJacobianMatrix[0] * transpJacobianMatrix[5]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[0] * transpJacobianMatrix[4] - transpJacobianMatrix[1] * transpJacobianMatrix[3]) / determinant;
        ++localJacobianMatrixElementIt;
    }

    void computeLocalJacobianMatrixies(const double(*transpJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                       const size_t nJacobianMatrixies,
                                       const double* determinnatIt,
                                       double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    { 
        for (size_t i = 0; i < nJacobianMatrixies; ++i)
        {
            computeLocalJacobianMatrix(*transpJacobianMatrixIt, *determinnatIt, *localJacobianMatrixIt);
            
            ++transpJacobianMatrixIt;
            ++determinnatIt;
            ++localJacobianMatrixIt;
        }
    }


    void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D* localGradientIt, const uint16_t nPoints, Coordinates* gradientIt)
	{
		for (uint16_t i = 0; i < nPoints; ++i)
		{
			gradientIt->x = localJacobian[0] * localGradientIt->u + localJacobian[1] * localGradientIt->v + localJacobian[2] * localGradientIt->w;
			gradientIt->y = localJacobian[3] * localGradientIt->u + localJacobian[4] * localGradientIt->v + localJacobian[5] * localGradientIt->w;
			gradientIt->z = localJacobian[6] * localGradientIt->u + localJacobian[7] * localGradientIt->v + localJacobian[8] * localGradientIt->w;

			++gradientIt;
			++localGradientIt;
		}
	}

    void translate0FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.v;
            localPoint3D.v = templateLocalPoint2D.u;
            localPoint3D.w = 0.0;
        }
    }

    void translate1FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.u;
            localPoint3D.v = 0.0;
            localPoint3D.w = templateLocalPoint2D.v;
        }
    }

    void translate2FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = 0.0;
            localPoint3D.v = templateLocalPoint2D.v;
            localPoint3D.w = templateLocalPoint2D.u;
        }
    }

    void translate3FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.u;
            localPoint3D.v = templateLocalPoint2D.v;
            localPoint3D.w = 1.0 - templateLocalPoint2D.u - templateLocalPoint2D.v;
        }
    }

    void translate(const LocalCoordinates3D localPoints[],
                   const uint8_t pointsCount,
                   const Coordinates& initPoint,
                   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                   Coordinates points[])
    {
        for (uint8_t i = 0; i < pointsCount; ++i)
        {
            Coordinates& point = points[i];
            point = initPoint;
            const LocalCoordinates3D& localPoint = localPoints[i];

            point.x += transpJacobianMatrix[0] * localPoint.u;
            point.y += transpJacobianMatrix[1] * localPoint.u;
            point.z += transpJacobianMatrix[2] * localPoint.u;

            point.x += transpJacobianMatrix[3] * localPoint.v;
            point.y += transpJacobianMatrix[4] * localPoint.v;
            point.z += transpJacobianMatrix[5] * localPoint.v;

            point.x += transpJacobianMatrix[6] * localPoint.w;
            point.y += transpJacobianMatrix[7] * localPoint.w;
            point.z += transpJacobianMatrix[8] * localPoint.w;
        }
    }

    void translate(const LocalCoordinates2D localPoints[],
                   const uint8_t pointsCount,
                   const Coordinates& initPoint,
                   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                   Coordinates points[])
    {
        for (uint8_t i = 0; i < pointsCount; ++i)
        {
            Coordinates& point = points[i];
            point = initPoint;
            const LocalCoordinates2D& localPoint = localPoints[i];

            point.x += transpJacobianMatrix[0] * localPoint.u;
            point.y += transpJacobianMatrix[1] * localPoint.u;
            point.z += transpJacobianMatrix[2] * localPoint.u;

            point.x += transpJacobianMatrix[3] * localPoint.v;
            point.y += transpJacobianMatrix[4] * localPoint.v;
            point.z += transpJacobianMatrix[5] * localPoint.v;
        }
    }

    void coomputeDirectionalDerivative(const Coordinates* gradientIt, const uint16_t nGradients, const Coordinates direction, double* directionalDerivativeIt)
    {
        for (uint8_t i = 0; i < nGradients; ++i)
        {
            *directionalDerivativeIt = gradientIt->x * direction.x + gradientIt->y * direction.y + gradientIt->z * direction.z;
            
            ++gradientIt;
            ++directionalDerivativeIt;
        }
    }

    static void computeDiffrence(const Coordinates& point1, const Coordinates& point2, Coordinates& vector)
    {
        vector.x = point1.x - point2.x;
        vector.y = point1.y - point2.y;
        vector.z = point1.z - point2.z;
    }

    void computeNormal(const Coordinates trianglePoints[constants::triangle::N_NODES], Coordinates &normal)
    {
        Coordinates vector1, vector2;

        computeDiffrence(trianglePoints[1], trianglePoints[0], vector1);
        computeDiffrence(trianglePoints[2], trianglePoints[0], vector2);

        normal.x = vector1.y * vector2.z - vector1.z * vector2.y;
        normal.y = vector1.z * vector2.x - vector1.x * vector2.z;
        normal.z = vector1.x * vector2.y - vector1.y * vector2.x;

        double norm = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x /= norm;
        normal.y /= norm;
        normal.z /= norm;
    }

}

