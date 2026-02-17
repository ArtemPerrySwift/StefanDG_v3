#pragma once
#include "MaterialPhase.h"
#include "BoundariesConditions.h"

struct Volume
{
	int tag;
	const MaterialPhase* materialPhasePtr;
	const int* boundariesTags;
	const BoundaryCondition* boundariesConditions;
	unsigned int nBoundaries;
};
