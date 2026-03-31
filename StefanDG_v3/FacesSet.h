#pragma once
#include "cstdint"

struct FacesSet
{
	union ConditionRelated
	{
		int entityTag;
		double* firstElementNormal;
	};

	const size_t indexes;
	//const size_t* elementIndexes;
	//const uint8_t* localindexes;
	size_t count;

	//size_t* tagsBuffer;
};
