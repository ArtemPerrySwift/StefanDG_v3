#pragma once
#include "Coordinates.h"
#include "LocalCoordinates3D.h"

struct InterfaceSideElements
{
	size_t* indexes;
    const uint8_t* facesLocalIndexes;
    const double (*localJacobians)[LocalCoordinates3D::COUNT * Coordinates::COUNT];
    const Coordinates* initPoints;
};
