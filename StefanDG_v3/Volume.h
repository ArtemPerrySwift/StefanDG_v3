#pragma once
#include "MaterialPhase.h"
#include "Boundary.h"

struct Volume
{
	int tag;
	const MaterialPhase* materialPhasePtr;
	Boundary* boundaries;
	//const int* boundariesTags;
	//const BoundaryCondition* boundariesConditions;
	unsigned int nBoundaries;
};
