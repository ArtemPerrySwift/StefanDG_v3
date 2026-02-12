#include "TetrahedronsAdjusments.h"


size_t TetrahedronsAdjusments::calcAllTrtrahedronsCount(const TetrahedronsAdjusments* tetrahedronsAdjusmentsSetIt, const size_t nSets)
{
    size_t nTetrahederons = 0;
    for (size_t i = 0; i < nSets; ++i)
    {
        nTetrahederons += tetrahedronsAdjusmentsSetIt->nTetrahedrons;
        ++tetrahedronsAdjusmentsSetIt;
    }
    
    return nTetrahederons;
}

