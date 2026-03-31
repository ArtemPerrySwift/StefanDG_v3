#include "CoordinatesFunctions.h"
#include <cmath>

namespace CoordinatesFunctions
{
    double computeTranspJacobianMatrix(const Coordinates* tetrahedronNodeBeginIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* directionNodesPtr = tetrahedronNodeBeginIt + 1;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        ++directionNodesPtr;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        ++directionNodesPtr;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        double det = transpJacobianMatrix[0] * transpJacobianMatrix[4] * transpJacobianMatrix[8];
        det += transpJacobianMatrix[2] * transpJacobianMatrix[3] * transpJacobianMatrix[7];
        det += transpJacobianMatrix[1] * transpJacobianMatrix[5] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[2] * transpJacobianMatrix[4] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[0] * transpJacobianMatrix[5] * transpJacobianMatrix[7];
        det -= transpJacobianMatrix[1] * transpJacobianMatrix[3] * transpJacobianMatrix[8];

        return det;
    }
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

    void copyTetrahedronsFacePoints(const uint8_t faceLocalIndex,
                                    const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES],
                                    Coordinates facePoints[constants::tetrahedron::N_NODES])
    {
        switch (faceLocalIndex)
        {
        case 0:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[2];
            *(++facePoints) = tetrahedronPoints[1];
            break;
        }
        case 1:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[1];
            *(++facePoints) = tetrahedronPoints[3];
            break;
        }
        case 2:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[3];
            *(++facePoints) = tetrahedronPoints[2];
            break;
        }
        case 3:
        {
            *facePoints = tetrahedronPoints[3];
            *(++facePoints) = tetrahedronPoints[1];
            *(++facePoints) = tetrahedronPoints[2];
            break;
        }
        }
    }

    void computeTetrahedronFaceDirections(const uint8_t faceLocalIndex, const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES], Coordinates directions[2])
    {
        switch (faceLocalIndex)
        {
        case 0:
        {
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 1:
        {
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[3], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 2:
        {
            computeDiffrence(tetrahedronPoints[3], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 3:
        {
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[3], directions[0]);
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[3], directions[1]);
            break;
        }
        }
    }

    double computeTriagnleJacobianDet(const Coordinates triangleDirections[2])
    {
        const Coordinates* triangleDirection2 = triangleDirections + 1;
        double e = triangleDirections->x * triangleDirections->x + triangleDirections->y * triangleDirections->y + triangleDirections->z * triangleDirections->z;
        double g = triangleDirection2->x * triangleDirection2->x + triangleDirection2->y * triangleDirection2->y + triangleDirection2->z * triangleDirection2->z;
        double f = triangleDirections->x * triangleDirection2->x + triangleDirections->y * triangleDirection2->y + triangleDirections->z * triangleDirection2->z;

        return sqrt(e * g - f * f);
    }

    void computeNormalViaTriangleDirections(const Coordinates triangleDirections[2], Coordinates& normal)
    {
        const Coordinates* triangleDirection2 = triangleDirections + 1;
        normal.x = triangleDirections->y * triangleDirection2->z - triangleDirections->z * triangleDirection2->y;
        normal.y = triangleDirections->z * triangleDirection2->x - triangleDirections->x * triangleDirection2->z;
        normal.z = triangleDirections->x * triangleDirection2->y - triangleDirections->y * triangleDirection2->x;

        double norm = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x /= norm;
        normal.y /= norm;
        normal.z /= norm;
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

