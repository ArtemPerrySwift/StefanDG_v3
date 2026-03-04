#pragma once
#include "MaterialPhase.h"
#include "Boundary.h"

struct Volume
{
	int tag;
	const MaterialPhase* materialPhasePtr;
	const Boundary* boundaries;
	//const int* boundariesTags;
	//const BoundaryCondition* boundariesConditions;
	unsigned int nBoundaries;
};
