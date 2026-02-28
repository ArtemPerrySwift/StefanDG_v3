#pragma once

struct Interface
{
	int tag;
	double thermalConductivity;
	size_t volumesIndexes[2];
	size_t sidesIndexes[2];
};
