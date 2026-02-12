#include "CrossTetrahedronsAdjusments.h"

size_t CrossTetrahedronsAdjusments::calcAllCrossElementsCount(const CrossTetrahedronsAdjusments* crossTetrahedronsAdjusmentsSetIt, const size_t nSets)
{
    size_t nCrossElements = 0;
    for (size_t i = 0; i < nSets; ++i)
    {
        nCrossElements += crossTetrahedronsAdjusmentsSetIt->nCrossElements;
        ++crossTetrahedronsAdjusmentsSetIt;
    }
    
    return nCrossElements;
}

