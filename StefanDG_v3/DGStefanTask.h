#pragma once
#include "MaterialPhase.h"
#include "TetrahedronsAdjusments.h"
#include "CrossTetrahedronsAdjusments.h"
#include "ElementDGCalculator.h"
#include "FaceDGCalculator.h"
#include "LinearLagrangeBasis.h"
//#include "BoundariesConditions.h"
#include "DTGeometryKernel.h"
#include "GMSHProxy.h"
#include "Volume.h"
#include "Interface.h"
#include "NonconformInterface.h"
#include "ConformInterface.h"
#include "InterfaceSideElements.h"
#include "Solution.h"
#include "Map.h"
#include <cstdlib>

#include <stdexcept>
#include <cinttypes>

using Basis = LinearLagrangeBasis;

class DGStefanTask
{
public:
    static double dirichletCondition(const Coordinates& point);
    static double newmanCondition(const Coordinates& point);
    static double initialConditionLiquid(const Coordinates &point);
    static double initialConditionSolid(const Coordinates& point);

    static void computeDirichletCondition(const Coordinates* pointIt, const uint8_t nPoints, double* valueIt);
    static void computeNewmanCondition(const Coordinates* pointIt, const uint8_t nPoints, double* valueIt);
    static void computeLiquidTemperature(const Coordinates* pointIt, const uint8_t nPoints, double* valueIt);
    static void computeSolidTemperature(const Coordinates* pointIt, const uint8_t nPoints, double* valueIt);

    void solve(const double penalty,
               const double tMin,
               const double tMax,
               const size_t nTSteps,
               Solution* solutionIt)
    {
        Volume* volumes;
        unsigned int nVolumes;

        NonconformInterface* nonconformInterfaces;
        unsigned int nNonconformInterfaces;

        ConformInterface* conformInterfaces;
        unsigned int nConformInterfaces;
        
        Boundary* boundaries;
        Boundary::Condition* conditions;
        unsigned int nSharedBoundaries;

        MaterialPhase* materialPhases;
        unsigned int nMaterialPhases;

        getModel(materialPhases, nMaterialPhases, volumes, nVolumes, boundaries, conditions, nonconformInterfaces, nNonconformInterfaces, conformInterfaces, nConformInterfaces, nSharedBoundaries);

        ElementsAdjusmentsSet* volumesElementsAdjusmentsSets = new ElementsAdjusmentsSet[nVolumes];
        CrossElementsAdjusmentsSet* volumesElementsCrossAdjuesmentsSets = new CrossElementsAdjusmentsSet[nVolumes];
        CrossElementsAdjusmentsSet* nonconformInterfacesAdjuesmentsSets = new CrossElementsAdjusmentsSet[nNonconformInterfaces];

        double** firstApproximationIt = new double*[nVolumes];
        InterfaceSideElementsSet(*nonconformalInterfaceSideElementsSets)[2] = new InterfaceSideElementsSet[nNonconformInterfaces][2];
        InterfaceSideElementsSet(*conformalInterfaceSideElementsSets)[2] = new InterfaceSideElementsSet[nConformInterfaces][2];

        size_t** confromalInterfacesFacesTags;
        size_t* confromalInterfacesFacesTagsCounts;

        const size_t nTIntervals = nTSteps - 1;
        double dt = (tMax - tMin) / nTIntervals;

        const Volume* volumeIt = volumes;
        ElementsAdjusmentsSet* elementsAdjusmentsIt = volumesElementsAdjusments;
        CrossElementsAdjusmentsSet* elementsCrossAdjuesmentsIt = volumesInteriorFacesAdjusments;
        double** firstApproximationIt = volumesInitialX;
        InterfaceSideElementsSet* interfaceSideElementsIt = interfacesSidesElements;
        for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
        {

        }

        for (size_t i = 1; i < nTIntervals; ++i)
        {

        }

        delete[] volumes;
        delete[] nonconformInterfaces;
        delete[] conformInterfaces;
        delete[] volumesElementsAdjusmentsSets;
        delete[] volumesElementsCrossAdjuesmentsSets;
        delete[] nonconformInterfacesAdjuesmentsSets;
        delete[] firstApproximationIt;
        delete[] nonconformalInterfaceSideElementsSets;
        delete[] boundaries;
        delete[] conditions;
    }


private:

    void getMaterialPhasesTags(const MaterialPhase materialPhases[2], int materialPhaseTagIt[2])
    {
        int* physicalDimTags;
        size_t nPhysicalDimTags;
        GMSHProxy::model::getPhysicalGroups(constants::DIMENSION_3D, &physicalDimTags, &nPhysicalDimTags);

        materialPhaseTagIt[0] = 0;
        materialPhaseTagIt[1] = 0;

        size_t nPhysicalGroups = nPhysicalDimTags / 2;
        const int* physicalDimTagIt = physicalDimTags;
        for (size_t i = 0; i < nPhysicalGroups; ++i)
        {
            char* physicalName;
            ++physicalDimTagIt;
            GMSHProxy::model::getPhysicalName(constants::DIMENSION_3D, *physicalDimTagIt, physicalName);

            if (!strcmp(materialPhases[0].name, physicalName))
            {
                materialPhaseTagIt[0] = *physicalDimTagIt;
            }
            else
            {
                if (!strcmp(materialPhases[1].name, physicalName))
                {
                    materialPhaseTagIt[1] = *physicalDimTagIt;
                }
            }

            ++physicalDimTagIt;
            GMSHProxy::free(physicalName);
        }
    }

    void assignBoundaryConditionToGroupBoundaries(const int* groupBoundariesTagIt,
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

    /*
    void assignBoundaryConditionToGroupSurfaces(const int* surfacesTagIt,
                                                const unsigned int nSurfaces,
                                                const Boundary::Condition *boundaryCondition,
                                                Map::Element<int, const Boundary::Condition*>* conditionTypeByBoundaryTag,
                                                const unsigned int nMaxSurfaces)
    {
        for (unsigned int j = 0; j < nSurfaces; ++j)
        {
            const unsigned int hashIndex = *surfacesTagIt % nMaxSurfaces;
            Map::Element<int, const Boundary::Condition*>* mapElement = conditionTypeByBoundaryTag + hashIndex;

            if (mapElement->tag == 0)
            {
                mapElement->tag = *surfacesTagIt;
                mapElement->value = boundaryCondition;

                return;
            }
            else
            {
                while (mapElement->next != nullptr)
                {
                    mapElement = mapElement->next;
                }

                mapElement->next = new Map::Element<int, const Boundary::Condition*>(*surfacesTagIt, boundaryCondition);
            }

            ++surfacesTagIt;
        }
    }
    unsigned int getConditionIndex(char strIndex[])
    {
        if (*strIndex == '\0')
        {
            return 0;
        }

        unsigned int index;
        sscanf(strIndex, "%u", &index);
        return index;
    }
    */
    void setBoundariesZeroTags(Boundary* boundaryIt, const unsigned int nBoundaries)
    {
        for (unsigned int i = 0; i < nBoundaries; ++i)
        {
            boundaryIt->tag = 0;
        }
    }

    void initilizeBoundariesMap(const int* boundariesDimTagIt, const unsigned int nBoundaries, Boundary* boundariesMap, unsigned int &nSharedBoundaries, unsigned int &nFreeSharedBoundaries)
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

    Boundary::Condition* setExplicitBoundariesConditions(const int* physicalDimTagIt,
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

    void setImpliticBoundaryConditions(Boundary* boundaryIt, const unsigned int nBoundaries, Boundary::Condition* conditionIt)
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

    void determineBoundaries(Boundary* &boundariesMap,
                             unsigned int &nBoundaries,
                             Boundary::Condition*& conditions,
                             unsigned int& nNonconformalInterfaces,
                             unsigned int& nConformalInterfaces,
                             unsigned int& nSharedBoundaries)
    {
        int* surfacesDimTags;
        GMSHProxy::model::getEntities(surfacesDimTags, nBoundaries, constants::DIMENSION_2D);
        boundariesMap = new Boundary[nBoundaries];

        setBoundariesZeroTags(boundariesMap, nBoundaries);

        initilizeBoundariesMap(surfacesDimTags, nBoundaries, boundariesMap, nSharedBoundaries, nConformalInterfaces);

        int* physicalDimTags;
        unsigned int nPhysicalGroups;
        GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalGroups, constants::DIMENSION_2D);

        Boundary::Condition* conditionIt = new Boundary::Condition[nPhysicalGroups + nConformalInterfaces + 1];
        conditions = conditionIt;

        unsigned int nProcessedBoundaries;
        conditionIt = setExplicitBoundariesConditions(physicalDimTags, nPhysicalGroups, boundariesMap, nBoundaries, conditionIt, nNonconformalInterfaces, nProcessedBoundaries);

        if (nProcessedBoundaries != nBoundaries)
        {
            setImpliticBoundaryConditions(boundariesMap, nBoundaries, conditionIt);
        }

    }
/*
    unsigned int setExplicitBoundariesConditions(Map::Element<int, const Boundary::Condition*>* &conditionTypeByBoundaryTag,
                                               unsigned int &nSurfaces,
                                               Boundary::Condition* &conditions,
                                               unsigned int &nNonconformalInterfaces)
    {
        int* surfacesDimTags;
        GMSHProxy::model::getEntities(surfacesDimTags, nSurfaces, constants::DIMENSION_2D);
        conditionTypeByBoundaryTag = new Map::Element<int, const Boundary::Condition*>[nSurfaces];
        
        int* physicalDimTags;
        unsigned int nPhysicalGroups;
        GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalGroups, constants::DIMENSION_2D);

        Boundary::Condition* conditionIt = new Boundary::Condition[nPhysicalGroups + 1];
        conditions = conditionIt;

        const int* physicalDimTagIt = physicalDimTags;
        char* physicalName;
        char* physicalNameCharIt;
        unsigned int nProcessedSurfaces = 0, iNonconformInterface = 0;
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

            assignBoundaryConditionToGroupSurfaces(groupSurfacesTags, nGroupSurfaces, conditionIt, conditionTypeByBoundaryTag, nSurfaces);

            nProcessedSurfaces += nGroupSurfaces;
            GMSHProxy::free(groupSurfacesTags);

            ++conditionIt;
            ++physicalDimTagIt;
        }

        GMSHProxy::free(physicalDimTags);
        nNonconformalInterfaces = iNonconformInterface;
        conditionIt->type = Boundary::Condition::Type::HOMOGENEOUS_NEWMAN;

        if (nProcessedSurfaces != nSurfaces)
        {
            const int* surfaceDimTagIt = surfacesDimTags;
            for (unsigned int i = 0; i < nSurfaces; ++i)
            {
                ++surfaceDimTagIt;
                Map::Element<int, const Boundary::Condition*>* mapElement = conditionTypeByBoundaryTag + (*surfaceDimTagIt) % nSurfaces;

                if (mapElement->tag != *surfaceDimTagIt)
                {
                    if (mapElement->tag == 0)
                    {
                        mapElement->tag = *surfaceDimTagIt;
                        mapElement->value = conditionIt;
                        return;
                    }
                    else
                    {
                        while (mapElement->next != nullptr && mapElement->tag != *surfaceDimTagIt)
                        {
                            mapElement = mapElement->next;
                        }

                        if (mapElement->tag != *surfaceDimTagIt)
                        {
                            mapElement->tag = *surfaceDimTagIt;
                            mapElement->value = conditionIt;

                        }

                    }
                }

                ++surfaceDimTagIt;
            }
        }

        GMSHProxy::free(surfacesDimTags);

        return nProcessedSurfaces;
    }
    */
    /*
    void determineImplicitBoundariesConditions(const int* surfaceDimTagIt,
                                               const unsigned int nSurfaces,
                                               const unsigned int nUnprocessedSurfaces,
                                               Map::Element<int, const Boundary::ICondition*>* conditionTypeByBoundaryTag,
                                               const Boundary::ICondition**& implicitConditions,
                                               unsigned int &nImplicitConditions)
    {
        const Boundary::ICondition* homogeneousNewmanCondition = new Boundary::Condition<Boundary::ICondition::Type::HOMOGENEOUS_NEWMAN>();
        const Boundary::ICondition** implicitConditionsBuffer = new const Boundary::ICondition*[nUnprocessedSurfaces];

        const Boundary::ICondition** implicitConditionIt = implicitConditionsBuffer;
        unsigned int nConformInterafces = 0;
        for (unsigned int i = 0; i < nSurfaces; ++i)
        {
            ++surfaceDimTagIt;
            Map::Element<int, const Boundary::ICondition*>* mapElement = conditionTypeByBoundaryTag + (*surfaceDimTagIt) % nSurfaces;

            if (mapElement->tag != *surfaceDimTagIt)
            {
                if (mapElement->tag == 0)
                {
                    mapElement->tag = *surfaceDimTagIt;

                    int* upEntitiesTags;
                    size_t nUpEntities;
                    GMSHProxy::model::getSurfaceUpEntities(*surfaceDimTagIt, upEntitiesTags, nUpEntities);
                    GMSHProxy::free(upEntitiesTags);

                    if (nUpEntities == 2)
                    {
                        *implicitConditionIt = new Boundary::Condition<Boundary::ICondition::Type::CONFORM_INTERFACE>(nConformInterafces);
                        mapElement->value = *implicitConditionIt;
                        ++implicitConditionIt;
                        ++nConformInterafces;
                    }
                    else
                    {
                        mapElement->value = homogeneousNewmanCondition;
                    }
                    return;
                }
                else
                {
                    while (mapElement->next != nullptr && mapElement->tag != *surfaceDimTagIt)
                    {
                        mapElement = mapElement->next;
                    }

                    if (mapElement->tag != *surfaceDimTagIt)
                    {
                        int* upEntitiesTags;
                        size_t nUpEntities;
                        GMSHProxy::model::getSurfaceUpEntities(*surfaceDimTagIt, upEntitiesTags, nUpEntities);
                        GMSHProxy::free(upEntitiesTags);

                        if (nUpEntities == 2)
                        {
                            *implicitConditionIt = new Boundary::Condition<Boundary::ICondition::Type::CONFORM_INTERFACE>(nConformInterafces);
                            mapElement->value = *implicitConditionIt;
                            ++implicitConditionIt;
                            ++nConformInterafces;
                        }
                        else
                        {
                            mapElement->value = homogeneousNewmanCondition;
                        }
                        
                    }

                }
            }

            ++surfaceDimTagIt;
        }

        unsigned int nImplicitConditions = nConformInterafces ;
        implicitConditions = new const Boundary::ICondition*[nConformInterafces];

        const Boundary::ICondition** implicitConditionsIt = implicitConditions;
        *implicitConditionsIt = homogeneousNewmanCondition;
        ++implicitConditionsIt;
        implicitConditionIt = implicitConditionsBuffer;

        implicitConditionIt = 

        for (unsigned int i = 0; i < nConformInterafces; ++i)
        {

        }
    }
    */

    /*
    void determineBoundariesConditions(Map::Element<int, Boundary::Condition>* &conditionTypeByBoundaryTag, unsigned int &nSurfaces, unsigned int &nNonconformalGroups)
    {
        int* surfacesDimTags;
        GMSHProxy::model::getEntities(surfacesDimTags, nSurfaces, constants::DIMENSION_2D);
        unsigned int nProcessedSurfaces = determineExplicitBoundariesConditions(conditionTypeByBoundaryTag, nSurfaces, nNonconformalGroups);

        if (nProcessedSurfaces != nSurfaces)
        {
            determineImplicitBoundariesConditions(surfacesDimTags, nSurfaces, conditionTypeByBoundaryTag);
        }

        GMSHProxy::free(surfacesDimTags);
    }
    */

    unsigned int determineVolumesMaterialPhase(const MaterialPhase* materialPhasePtr, Map::Element<int, const MaterialPhase*>* materialPhasePtrByVolumeTag, const unsigned int nVolumes)
    {
        int* dimTags;
        unsigned int nPhaseVolumes;
        GMSHProxy::model::getEntitiesForPhysicalName(materialPhasePtr->name, dimTags, nPhaseVolumes);

        const int* dimTagIt = dimTags;
        for (unsigned int i = 0; i < nPhaseVolumes; ++i)
        {
            ++dimTagIt;
            Map::Element<int, const MaterialPhase*>* materialPhasePtrByVolumeTagElement = materialPhasePtrByVolumeTag + (*dimTagIt) % nVolumes;
            if (materialPhasePtrByVolumeTagElement->tag == 0)
            {
                materialPhasePtrByVolumeTagElement->tag = *dimTagIt;
                materialPhasePtrByVolumeTagElement->value = materialPhasePtr;
            }
            else
            {
                while (materialPhasePtrByVolumeTag->next != nullptr)
                {
                    materialPhasePtrByVolumeTag = materialPhasePtrByVolumeTag->next;
                }

                materialPhasePtrByVolumeTag->next = new Map::Element<int, const MaterialPhase*>(*dimTagIt, materialPhasePtr);
            }
        }

        GMSHProxy::free(dimTags);

        return nPhaseVolumes;
    }

    void determineVolumesMaterialPhases(const MaterialPhase materialPhases[2], Map::Element<int, const MaterialPhase*>* materialPhasePtrByVolumeTag, const unsigned int nVolumes)
    {
        unsigned int nProcessedVolumes = determineVolumesMaterialPhase(materialPhases, materialPhasePtrByVolumeTag, nVolumes);
        nProcessedVolumes += determineVolumesMaterialPhase(materialPhases + 1, materialPhasePtrByVolumeTag, nVolumes);

        if (nProcessedVolumes != nVolumes)
        {
            throw std::runtime_error("For not all volumes material phase is set");
        }
    }

    void determineVolumesMaterials(std::pair<int, const MaterialPhase*>* const materialPhaseByVolumeTag, const unsigned int nVolumes, MaterialPhase* &materialPhases, unsigned int &nMaterialPhases)
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


    void getModel(MaterialPhase* &materialPhases,
                  unsigned int &nMaterialPhases,
                  Volume*& volumes,
                  unsigned int& nVolumes,
                  Boundary* &volumesBoundaries,
                  Boundary::Condition*& conditions,
                  NonconformInterface*& nonconformInterfaces,
                  unsigned int& nNonconformInterfaces,
                  ConformInterface* conformInterfaces,
                  unsigned int nConformInterfaces,
                  unsigned int nSharedBoundaries)
    {
        Boundary* boundariesMap;
        unsigned int nBoundaries;

        determineBoundaries(boundariesMap, nBoundaries, conditions, nNonconformInterfaces, nConformInterfaces, nSharedBoundaries);

        nonconformInterfaces = new NonconformInterface[nNonconformInterfaces];
        bool* wasNonconfromInterfaceProcessed = new bool[nNonconformInterfaces]{};

        ConformInterface* conformInterfaces = new ConformInterface[nConformInterfaces];
        bool* wasConfromInterfaceProcessed = new bool[nNonconformInterfaces]{};

        int* volumesDimTags;
        int** volumesBoundariesDimTags;
        unsigned int* volumesBoundariesCounts;
        unsigned int nVolumesBoundaries;
        GMSHProxy::model::getVolumes(volumesDimTags, nVolumes, volumesBoundariesDimTags, volumesBoundariesCounts, nVolumesBoundaries);

        Volume* volumeIt = new Volume[nVolumes];
        volumes = volumeIt;

        std::pair<int, const MaterialPhase*>* materialPhaseByVolumeTag = new std::pair<int, const MaterialPhase*>[nVolumes];
        determineVolumesMaterials(materialPhaseByVolumeTag, nVolumes, materialPhases, nMaterialPhases);

        Boundary* volumeBoundaryIt = new Boundary[nVolumesBoundaries];
        volumesBoundaries = volumeBoundaryIt;

        const int* volumesDimTagsIt = volumesDimTags;
        const unsigned int* nVolumeBoundariesIt = volumesBoundariesCounts;
        int** volumeBoundariesDimTagIt  = volumesBoundariesDimTags;
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
                    NonconformInterface* nonconformInterface = nonconformInterfaces + index;
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

                if (boundaryPtr->condition->type == Boundary::Condition::Type::CONFORM_INTERFACE)
                {
                    unsigned int index = boundaryPtr->condition->typeData.index;
                    ConformInterface* nonconformInterface = conformInterfaces + index;
                    if (wasNonconfromInterfaceProcessed[index])
                    {
                        nonconformInterface->volumesIndexes[1] = i;
                    }
                    else
                    {
                        nonconformInterface->thermalConductivity = volumeIt->materialPhasePtr->thermalConductivity;
                        nonconformInterface->volumesIndexes[0] = i;
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
        delete[] wasConfromInterfaceProcessed;

        GMSHProxy::free(volumesDimTags);
        GMSHProxy::free(volumesBoundariesDimTags);
        GMSHProxy::free(volumesBoundariesCounts);
    }

    double* dblBuffer;
    int* intBuffer;
    Coordinates* coordinatesBuffer;
    void** ptrBuffer;
    size_t* sizesBuffer;

    double _penalty;

	MaterialPhase* materialPhases;

    Volume* volumes;
	size_t nVolumes;

	ElementsAdjusmentsSet* volumesElementsAdjusments;
	CrossElementsAdjusmentsSet* volumesInteriorFacesAdjusments;
	double** volumesInitialX;

    Interface* nonconformalInterfaces;
	size_t nNonconformalInterfaces;

    InterfaceSideElementsSet* interfacesSidesElements;
	CrossElementsAdjusmentsSet* nonconformalInterfacesFragmentsAdjusments;


    void computeTetrahedronsAdjusments(const double (*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                       const double (*transpMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                       const double* determinantIt,
                                       const Coordinates* initPointIt,
                                       const size_t nTetrahedrons,
                                       const double lambda,
                                       const double gamma,
                                       const double dt,
                                       double (*initialCondition)(const Coordinates& point),
                                       double (*initialXIt)[Basis::N_FUNCTIONS],
                                       double (*bilinearTetrahedronsAdjusmentsIt)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                       double (*linearTetrahedronsAdjusmentsIt)[Basis::N_FUNCTIONS]);

    void computeInteriorFacesAdjusments(const int interiorFacesEntityTag,
                                        const size_t(*interiorFaceTetrahedronsIndexesIt)[2],
                                        const uint8_t(*interiorFaceLocalIndexesIt)[2],
                                        const double lambda,
                                        const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                        double (bilinearTetrahedronsAdjusments)[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                        double (*bilinearInterfaceAdjusmentsIt)[FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS]);

    void computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                         const size_t* boundaryFaceTetrahedronIndexIt,
                                         const uint8_t* boundaryFacesLocalIndexesIt,
                                         const Coordinates outwardNormal,
                                         const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         double(* const dirichletCondition)(const Coordinates&),
                                         const double lambda,
                                         double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]);

    void computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                         const size_t* boundaryFaceTetrahedronIndexIt,
                                         const uint8_t* boundaryFacesLocalIndexesIt,
                                         const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         double(* const dirichletCondition)(const Coordinates&),
                                         const double lambda,
                                         double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]);

    void computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                         const size_t* boundaryFaceTetrahedronIndexIt,
                                         const uint8_t* boundaryFacesLocalIndexesIt,
                                         const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         const double condition,
                                         const double lambda,
                                         double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]);

    void computeNewmanFacesAdjusments(const int boundaryFacesEntityTag,
                                      const size_t* boundaryFacesTetrahedronIndexesIt,
                                      const uint8_t* boundaryFacesLocalIndexesIt,
                                      double(* const newmanCondition)(const Coordinates&),
                                      double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]);


    void processVolumes(const double dt);


};

