#pragma once
#include "CommonFragmentAdjusmentsBase.h"
#include <cstdint>

template<uint8_t N_DOFS1, uint8_t N_DOFS2>
struct CommonFragmentAdjusments : CommonFragmentAdjusmentsBase
{
	double bilinear[N_DOFS1 * N_DOFS2];
};
