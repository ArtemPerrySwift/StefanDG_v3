#include "CoordinatesFunctions.h"

namespace CoordinatesFunctions
{
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
}

