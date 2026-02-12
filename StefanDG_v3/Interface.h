#pragma once

struct Interface
{
	int tag;
	double thermalConductivity;
	size_t volumesIndexes[2];
	int sidesTags[2];
};
