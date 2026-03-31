#pragma once
#include "Coordinates.h"
#include "LocalCoordinates3D.h"

struct InterfaceSideElementsSet
{
	const size_t* indexes;
    const double (*localJacobians)[LocalCoordinates3D::COUNT * Coordinates::COUNT];
    const Coordinates* initPoints;
};
