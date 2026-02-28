#pragma once
#include "cstdint"

struct FacesSet
{
	const size_t* nodesTags;
	const size_t* elementIndexes;
	const uint8_t* localindexes;
	size_t count;
};
