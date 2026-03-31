#pragma once
#include "Boundary.h"

struct SharedBoundary
{
	int tag;
	bool isPlane;
	const Boundary::Condition* condition;

	size_t volumesIndexes[2];

	struct MeshBuffer
	{
		size_t* facesTags;
		size_t nFaces;

		double* det;
		Coordinates* normals;
	};

};
