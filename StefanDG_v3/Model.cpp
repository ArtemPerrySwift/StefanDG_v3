#include "Model.h"
#include "GMSHProxy.h"
#include <stdexcept>

static void setBoundariesZeroTags(Boundary* boundaryIt, const unsigned int nBoundaries)
{
    for (unsigned int i = 0; i < nBoundaries; ++i)
    {
        boundaryIt->tag = 0;
    }
}

static void initilizeBoundariesMap(const int* boundariesDimTagIt,
                                   const unsigned int nBoundaries,
                                   Boundary* boundariesMap,
                                   unsigned int& nSharedBoundaries,
                                   unsigned int& nFreeSharedBoundaries)
{
    Boundary* boundaryPtr;
    unsigned int j;
    unsigned int iSharedBoundary = 0;
    unsigned int iFreeSharedBoundary = 0;

    for (unsigned int i = 0; i < nBoundaries; ++i)
    {
        ++boundariesDimTagIt;
        j = *boundariesDimTagIt % nBoundaries;
        boundaryPtr = boundariesMap + j;

        if (boundaryPtr->tag != 0)
        {
            ++j;
            while (boundaryPtr->tag != 0 && j < nBoundaries)
            {
                ++boundaryPtr;
                ++j;
            }

            if (boundaryPtr->tag != 0)
            {
                boundaryPtr = boundariesMap;
                while (boundaryPtr->tag != 0)
                {
                    ++boundaryPtr;
                }
            }
        }

        boundaryPtr->tag = *boundariesDimTagIt;
        boundaryPtr->condition = nullptr;
        boundaryPtr->isPlane = GMSHProxy::model::isSurfacePlane(boundaryPtr->tag);
        boundaryPtr->isShared = GMSHProxy::model::getSurfaceUpEntitiesCount(boundaryPtr->tag) > 1;
        if (boundaryPtr->isShared)
        {
            boundaryPtr->index = iSharedBoundary;
            ++iSharedBoundary;

            if (GMSHProxy::model::isSurfaceFreeOfPhysicalGroups(boundaryPtr->tag))
            {
                ++iFreeSharedBoundary;
            }
        }

        ++boundariesDimTagIt;
    }

    nSharedBoundaries = iSharedBoundary;
    nFreeSharedBoundaries = iFreeSharedBoundary;
}

static void assignBoundaryConditionToGroupBoundaries(const int* groupBoundariesTagIt,
                                              const unsigned int nGroupBoundaries,
                                              Boundary* const boundariesMap,
                                              const Boundary::Condition* conditionPtr,
                                              const unsigned int nBoundaries)
{
    for (unsigned int j = 0; j < nGroupBoundaries; ++j)
    {
        Boundary* boundaryPtr = boundariesMap + *groupBoundariesTagIt % nBoundaries;

        if (boundaryPtr->tag != *groupBoundariesTagIt)
        {
            ++j;
            while (boundaryPtr->tag != *groupBoundariesTagIt && j < nBoundaries)
            {
                ++boundaryPtr;
                ++j;
            }

            if (boundaryPtr->tag != *groupBoundariesTagIt)
            {
                boundaryPtr = boundariesMap;
                while (boundaryPtr->tag != *groupBoundariesTagIt)
                {
                    ++boundaryPtr;
                }
            }
        }

        boundaryPtr->condition = conditionPtr;
    }
}

static Boundary::Condition* setExplicitBoundariesConditions(const int* physicalDimTagIt,
                                                            const unsigned int nPhysicalGroups,
                                                            Boundary* boundariesMap,
                                                            const unsigned int nBoundaries,
                                                            Boundary::Condition* conditionIt,
                                                            unsigned int& nNonconformInterfaces,
                                                            unsigned int& nProcessedSurfaces)
{
    char* physicalName;
    char* physicalNameCharIt;
    unsigned int nProcessedSurfacesSum = 0, iNonconformInterface = 0;

    for (unsigned int i = 0; i < nPhysicalGroups; ++i)
    {
        ++physicalDimTagIt;

        int* groupSurfacesTags;
        size_t nGroupSurfaces;
        GMSHProxy::model::getEntitiesForPhysicalGroup(constants::DIMENSION_2D, *physicalDimTagIt, groupSurfacesTags, nGroupSurfaces);

        GMSHProxy::model::getPhysicalName(constants::DIMENSION_2D, *physicalDimTagIt, physicalName);
        physicalNameCharIt = physicalName;

        switch (*physicalNameCharIt)
        {
        case 'D':
        {
            ++physicalNameCharIt;
            if (*physicalNameCharIt == ':')
            {
                ++physicalNameCharIt;

                char* endptr;
                double value = strtod(physicalNameCharIt, &endptr);
                if (endptr == physicalNameCharIt)
                {
                    throw std::runtime_error("Dirichlet value is not set in physical name");
                }
                conditionIt->type = Boundary::Condition::Type::DIRICHLET_VALUE;
                conditionIt->typeData.value = value;

            }
            else
            {
                conditionIt->type = Boundary::Condition::Type::DIRICHLET_FUNCTION;
            }
            break;
        }
        case 'N':
        {
            ++physicalNameCharIt;
            if (*physicalNameCharIt == ':')
            {
                ++physicalNameCharIt;
                char* endptr;
                double value = strtod(physicalNameCharIt, &endptr);
                if (endptr == physicalNameCharIt)
                {
                    throw std::runtime_error("Newman value is not set in physical name");
                }
                conditionIt->type = Boundary::Condition::Type::NEWMAN_VALUE;
                conditionIt->typeData.value = value;
            }
            else
            {
                conditionIt->type = Boundary::Condition::Type::NEWMAN_FUNCTION;
            }
            break;
        }
        case 'C':
        {
            ++iNonconformInterface;
            if (nGroupSurfaces != 2)
            {
                throw std::runtime_error("Each condition of nonconformity should include only 2 surfaces");
            }

            conditionIt->type = Boundary::Condition::Type::NONCONFORM_INTERFACE;
            conditionIt->typeData.index = iNonconformInterface;
            ++iNonconformInterface;
            break;
        }
        case 'S':
        {
            ++physicalNameCharIt;
            if (*physicalNameCharIt == ':')
            {
                ++physicalNameCharIt;

                char* endptr;
                double value = strtod(physicalNameCharIt, &endptr);
                if (endptr == physicalNameCharIt)
                {
                    throw std::runtime_error("Stefan value is not set in physical name");
                }

                conditionIt->type = Boundary::Condition::Type::STEFAN;
                conditionIt->typeData.value = value;
            }
            else
            {
                throw std::runtime_error("Stefan name should contain \':\' after \'S\'");
            }
            break;
        }

        default:
        {
            throw std::runtime_error("Unknown type of condition");
        }
        }

        assignBoundaryConditionToGroupBoundaries(groupSurfacesTags, nGroupSurfaces, boundariesMap, conditionIt, nBoundaries);

        nProcessedSurfacesSum += nGroupSurfaces;
        GMSHProxy::free(groupSurfacesTags);

        ++conditionIt;
        ++physicalDimTagIt;
    }

    nNonconformInterfaces = iNonconformInterface;
    nProcessedSurfaces = nProcessedSurfacesSum;
    return conditionIt;
}

static void setImpliticBoundaryConditions(Boundary* boundaryIt, const unsigned int nBoundaries, Boundary::Condition* conditionIt)
{
    conditionIt->type = Boundary::Condition::Type::HOMOGENEOUS_NEWMAN;
    Boundary::Condition* homogeneousNewmanCondition = conditionIt;
    ++conditionIt;
    unsigned int iConformInterface = 0;
    for (unsigned int i = 0; i < nBoundaries; ++i)
    {

        if (boundaryIt->condition == nullptr)
        {
            if (boundaryIt->isShared)
            {
                conditionIt->type = Boundary::Condition::Type::CONFORM_INTERFACE;
                conditionIt->typeData.index = iConformInterface;

                boundaryIt->condition = conditionIt;
                ++iConformInterface;
            }
            else
            {
                boundaryIt->condition = homogeneousNewmanCondition;
            }
        }

        ++boundaryIt;
    }
}

void Model::determineBoundaries(Boundary*& boundariesMap, unsigned int& nBoundaries)
{
    int* surfacesDimTags;
    GMSHProxy::model::getEntities(surfacesDimTags, nBoundaries, constants::DIMENSION_2D);
    boundariesMap = new Boundary[nBoundaries];

    setBoundariesZeroTags(boundariesMap, nBoundaries);

    initilizeBoundariesMap(surfacesDimTags, nBoundaries, boundariesMap, nSharedBoundaries, nConformInterfaces);

    int* physicalDimTags;
    unsigned int nPhysicalGroups;
    GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalGroups, constants::DIMENSION_2D);

    Boundary::Condition* conditionIt = new Boundary::Condition[nPhysicalGroups + nConformInterfaces + 1];
    conditions = conditionIt;

    unsigned int nProcessedBoundaries;
    conditionIt = setExplicitBoundariesConditions(physicalDimTags, nPhysicalGroups, boundariesMap, nBoundaries, conditionIt, nNonconformInterfaces, nProcessedBoundaries);

    if (nProcessedBoundaries != nBoundaries)
    {
        setImpliticBoundaryConditions(boundariesMap, nBoundaries, conditionIt);
    }

}

void Model::determineVolumesMaterials(std::pair<int, const MaterialPhase*>* const materialPhaseByVolumeTag, const unsigned int nVolumes)
{
    std::pair<int, const MaterialPhase*>* materialPhaseByVolumeTagIt = materialPhaseByVolumeTag;
    for (unsigned int i = 0; i < nVolumes; ++i)
    {
        materialPhaseByVolumeTagIt->first = 0.0;
        ++materialPhaseByVolumeTagIt;
    }

    int* physicalGroupsDimTags;
    GMSHProxy::model::getPhysicalGroups(physicalGroupsDimTags, nMaterialPhases, constants::DIMENSION_3D);

    MaterialPhase*& materialPhaseIt = new MaterialPhase[nMaterialPhases];
    materialPhases = materialPhaseIt;

    const int* physicalGroupsDimTagIt = physicalGroupsDimTags;
    for (unsigned int i = 0; i < nMaterialPhases; ++i)
    {
        ++physicalGroupsDimTagIt;
        char* physicalName;
        GMSHProxy::model::getPhysicalName3D(*physicalGroupsDimTagIt, physicalName);

        int nReadFields = sscanf(physicalName, "M: %lf; %lf; %lf;", &materialPhaseIt->thermalConductivity, &materialPhaseIt->heatCapacity, &materialPhaseIt->density);
        if (nReadFields != 3)
        {
            throw std::exception("Incorrect material name");
        }

        GMSHProxy::free(physicalName);

        int* groupVolumesTags;
        unsigned int nGroupVolumes;
        GMSHProxy::model::getEntitiesForPhysicalGroup3D(*physicalGroupsDimTagIt, groupVolumesTags, nGroupVolumes);

        const int* groupVolumeTagIt = groupVolumesTags;
        for (unsigned int j = 0; j < nGroupVolumes; ++j)
        {
            unsigned int index = *groupVolumeTagIt % nVolumes;
            materialPhaseByVolumeTagIt = materialPhaseByVolumeTag + index;
            if (materialPhaseByVolumeTagIt->first != 0)
            {
                ++index;
                while (materialPhaseByVolumeTagIt->first != *groupVolumeTagIt && index < nVolumes)
                {
                    ++materialPhaseByVolumeTagIt;
                    ++index;
                }

                if (materialPhaseByVolumeTagIt->first != *groupVolumeTagIt)
                {
                    materialPhaseByVolumeTagIt = materialPhaseByVolumeTag;
                    while (materialPhaseByVolumeTagIt->first != *groupVolumeTagIt)
                    {
                        ++materialPhaseByVolumeTagIt;
                    }
                }
            }

            materialPhaseByVolumeTagIt->second = materialPhaseIt;
        }

        ++physicalGroupsDimTagIt;
    }

    GMSHProxy::free(physicalGroupsDimTags);
}

void Model::initilizeByCurrentGMSHModel()
{
    Boundary* boundariesMap;
    unsigned int nBoundaries;

    determineBoundaries(boundariesMap, nBoundaries);

    NonconformInterface* nonconformInterfacesBegin = new NonconformInterface[nNonconformInterfaces];
    nonconformInterfaces = nonconformInterfacesBegin;
    bool* wasNonconfromInterfaceProcessed = new bool[nNonconformInterfaces] {};

    SharedBoundary* sharedBoundaryBegin = new SharedBoundary[nSharedBoundaries];
    sharedBoundaries = sharedBoundaryBegin;
    bool* wasSharedBoundaryProcessed = new bool[nSharedBoundaries] {};

    unsigned nVolumes;
    int* volumesDimTags;
    int** volumesBoundariesDimTags;
    unsigned int* volumesBoundariesCounts;
    unsigned int nVolumesBoundaries;
    GMSHProxy::model::getVolumes(volumesDimTags, nVolumes, volumesBoundariesDimTags, volumesBoundariesCounts, nVolumesBoundaries);

    Volume* volumeIt = new Volume[nVolumes];
    volumes = volumeIt;

    std::pair<int, const MaterialPhase*>* materialPhaseByVolumeTag = new std::pair<int, const MaterialPhase*>[nVolumes];
    determineVolumesMaterials(materialPhaseByVolumeTag, nVolumes);

    Boundary* volumeBoundaryIt = new Boundary[nVolumesBoundaries];
    volumesBoundaries = volumeBoundaryIt;

    const int* volumesDimTagsIt = volumesDimTags;
    const unsigned int* nVolumeBoundariesIt = volumesBoundariesCounts;
    int** volumeBoundariesDimTagIt = volumesBoundariesDimTags;
    for (unsigned int i = 0; i < nVolumes; ++i)
    {
        ++volumesDimTagsIt;

        int* boundariesDimTags;
        unsigned int nVolumeBoundaries = *nVolumeBoundariesIt;
        GMSHProxy::model::getBoundariesFor3DEntity(*volumesDimTagsIt, boundariesDimTags, nVolumeBoundaries);
        GMSHProxy::free(boundariesDimTags);

        volumeIt->tag = *volumesDimTagsIt;

        unsigned int index = volumeIt->tag % nVolumes;
        std::pair<int, const MaterialPhase*>* materialPhaseByVolumeTagIt = materialPhaseByVolumeTag + index;
        if (materialPhaseByVolumeTagIt->first != 0)
        {
            ++index;
            while (materialPhaseByVolumeTagIt->first != volumeIt->tag && index < nVolumes)
            {
                ++materialPhaseByVolumeTagIt;
                ++index;
            }

            if (materialPhaseByVolumeTagIt->first != volumeIt->tag)
            {
                materialPhaseByVolumeTagIt = materialPhaseByVolumeTag;
                while (materialPhaseByVolumeTagIt->first != volumeIt->tag)
                {
                    ++materialPhaseByVolumeTagIt;
                }
            }
        }

        volumeIt->materialPhasePtr = materialPhaseByVolumeTagIt->second;

        volumeIt->boundaries = volumeBoundaryIt;
        volumeIt->nBoundaries = nVolumeBoundaries;

        const int* boundaryDimTagIt = *volumeBoundariesDimTagIt;
        for (unsigned int j = 0; j < nVolumeBoundaries; ++j)
        {
            ++boundaryDimTagIt;
            int boundaryTag = *boundaryDimTagIt;
            int boundaryIndex = boundaryTag % nBoundaries;
            Boundary* boundaryPtr = boundariesMap + boundaryIndex;

            if (boundaryPtr->tag != boundaryTag)
            {
                ++boundaryIndex;
                while (boundaryPtr->tag != boundaryTag && boundaryIndex < nBoundaries)
                {
                    ++boundaryPtr;
                    ++boundaryIndex;
                }

                if (boundaryPtr->tag != boundaryTag)
                {
                    boundaryPtr = boundariesMap;
                    while (boundaryPtr->tag != boundaryTag)
                    {
                        ++boundaryPtr;
                    }
                }
            }

            *volumeBoundaryIt = *boundaryPtr;

            if (boundaryPtr->condition->type == Boundary::Condition::Type::NONCONFORM_INTERFACE)
            {
                unsigned int index = boundaryPtr->condition->typeData.index;
                NonconformInterface* nonconformInterface = nonconformInterfacesBegin + index;
                if (wasNonconfromInterfaceProcessed[index])
                {
                    nonconformInterface->sidesTags[1] = boundaryTag;
                    nonconformInterface->volumesIndexes[1] = i;
                }
                else
                {
                    nonconformInterface->thermalConductivity = volumeIt->materialPhasePtr->thermalConductivity;
                    nonconformInterface->sidesTags[0] = boundaryTag;
                    nonconformInterface->volumesIndexes[0] = i;
                }
            }

            if (boundaryPtr->isShared)
            {
                unsigned int index = boundaryPtr->condition->typeData.index;
                SharedBoundary* sharedBoundaryPtr = sharedBoundaryBegin + boundaryPtr->condition->typeData.index;
                if (wasSharedBoundaryProcessed[index])
                {
                    sharedBoundaryPtr->volumesIndexes[1] = j;
                }
                else
                {
                    sharedBoundaryPtr->tag = boundaryPtr->tag;
                    sharedBoundaryPtr->isPlane = boundaryPtr->isPlane;
                    sharedBoundaryPtr->condition = boundaryPtr->condition;
                    sharedBoundaryPtr->volumesIndexes[0] = j;
                }
            }

            ++boundaryDimTagIt;
            ++volumeBoundaryIt;
        }

        GMSHProxy::free(*volumeBoundariesDimTagIt);

        ++nVolumeBoundariesIt;
    }

    delete[] materialPhaseByVolumeTag;
    delete[] boundariesMap;
    delete[] wasNonconfromInterfaceProcessed;
    delete[] wasSharedBoundaryProcessed;

    GMSHProxy::free(volumesDimTags);
    GMSHProxy::free(volumesBoundariesDimTags);
    GMSHProxy::free(volumesBoundariesCounts);

}

void Model::clear()
{
    delete[] materialPhases;
    delete[] volumes;
    delete[] volumesBoundaries;
    delete[] conditions;
    delete[] conformInterfaces;
    delete[] nonconformInterfaces;
}
