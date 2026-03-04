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

    void solve(const MaterialPhase* materialPhases,
               const double stefanTemperature,
               const double penalty,
               const double tMin,
               const double tMax,
               const size_t nTSteps,
               Solution* solutionIt)
    {
        Volume* volumes;
        NonconformInterface* nonconformInterfaces;
        unsigned int nVolumes;
        
        Boundary* boundaries;
        Boundary::ConditionsSet conditionsSet;
        getModel(materialPhases, volumes, nVolumes, boundaries, conditionsSet, nonconformInterfaces);
        unsigned int nNonconformInterafaces = conditionsSet.nNonconformInterfaces;

        ElementsAdjusmentsSet* volumesElementsAdjusmentsSets = new ElementsAdjusmentsSet[nVolumes];
        CrossElementsAdjusmentsSet* volumesElementsCrossAdjuesmentsSets = new CrossElementsAdjusmentsSet[nVolumes];
        CrossElementsAdjusmentsSet* nonconformInterfacesAdjuesmentsSets = new CrossElementsAdjusmentsSet[nNonconformInterafaces];

        double** firstApproximationIt = new double*[nVolumes];
        InterfaceSideElementsSet(*nonconformalInterfaceSideElementsSets)[2] = new InterfaceSideElementsSet[nNonconformInterafaces][2];

        const size_t nTIntervals = nTSteps - 1;
        double dt = (tMax - tMin) / nTIntervals;

        const Volume* volumeIt = volumes;
        ElementsAdjusmentsSet* elementsAdjusmentsIt = volumesElementsAdjusments;
        CrossElementsAdjusmentsSet* elementsCrossAdjuesmentsIt = volumesInteriorFacesAdjusments;
        double** firstApproximationIt = volumesInitialX;
        InterfaceSideElementsSet* interfaceSideElementsIt = interfacesSidesElements;
        for (unsigned int i = 0; i < nNonconformInterafaces; ++i)
        {

        }

        for (size_t i = 1; i < nTIntervals; ++i)
        {

        }

        delete[] volumes;
        delete[] nonconformInterfaces;
        delete[] volumesElementsAdjusmentsSets;
        delete[] volumesElementsCrossAdjuesmentsSets;
        delete[] nonconformInterfacesAdjuesmentsSets;
        delete[] firstApproximationIt;
        delete[] nonconformalInterfaceSideElementsSets;
        delete[] conditionsSet.valuesDirichlet;
        delete[] conditionsSet.valuesNewman;
        delete[] conditionsSet.nonconformInterfaces;
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


    void assignBoundaryConditionToGroupSurfaces(const int* surfacesTagIt,
                                                const unsigned int nSurfaces,
                                                const Boundary::ICondition *boundaryCondition,
                                                Map::Element<int, const Boundary::ICondition*>* conditionTypeByBoundaryTag,
                                                const unsigned int nMaxSurfaces)
    {
        for (unsigned int j = 0; j < nSurfaces; ++j)
        {
            const unsigned int hashIndex = *surfacesTagIt % nMaxSurfaces;
            Map::Element<int, const Boundary::ICondition*>* mapElement = conditionTypeByBoundaryTag + hashIndex;

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

                mapElement->next = new Map::Element<int, const Boundary::ICondition*>(*surfacesTagIt, boundaryCondition);
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

    unsigned int determineBoundariesConditions(Map::Element<int, const Boundary::ICondition*>* &conditionTypeByBoundaryTag,
                                               unsigned int &nSurfaces,
                                               Boundary::ConditionsSet& conditionsSet)
    {
        int* physicalDimTags;
        unsigned int nPhysicalGroups;
        GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalGroups, constants::DIMENSION_2D);

        char** physicalNames = new char*[nPhysicalGroups];

        unsigned int iValueDirichlet = 0;
        unsigned int iValueNewman = 0;
        unsigned int iNonconformInterface = 0;

        unsigned int nProcessedSurfaces = 0;

        const int* physicalDimTagIt = physicalDimTags;
        char** physicalNameIt = physicalNames;
        for (unsigned int i = 0; i < nPhysicalGroups; ++i)
        {
            ++physicalDimTagIt;
            GMSHProxy::model::getPhysicalName(constants::DIMENSION_2D, *physicalDimTagIt, *physicalNameIt);
            char* physicalNameCharIt = *physicalNameIt;
            
            switch (*physicalNameCharIt)
            {
            case 'D':
            {
                ++physicalNameCharIt;
                if (*physicalNameCharIt == 'V')
                {
                    ++iValueDirichlet;
                }
                break;
            }
            case 'N':
            {
                ++physicalNameCharIt;
                if (*physicalNameCharIt == 'V')
                {
                    ++iValueNewman;
                }
                break;
            }
            case 'C':
            {
                ++physicalNameCharIt;
                ++iNonconformInterface;
                break;
            }
            case 'S':
            {
                break;
            }
            default:
            {
                throw std::runtime_error("Unknown type of condition");
            }
            }
            ++physicalDimTagIt;
            ++physicalNameIt;
        }

        GMSHProxy::free(physicalDimTags);

        Boundary::Condition<Boundary::ICondition::Type::DIRICHLET_VALUE>* valueDirichletIt = new Boundary::Condition<Boundary::ICondition::Type::DIRICHLET_VALUE>[iValueDirichlet];
        Boundary::Condition<Boundary::ICondition::Type::NEWMAN_VALUE>* valueNewmanIt = new Boundary::Condition<Boundary::ICondition::Type::NEWMAN_VALUE>[iValueNewman];
        Boundary::Condition<Boundary::ICondition::Type::NONCONFORM_INTERFACE>* nonconformInterfaceIt = new Boundary::Condition<Boundary::ICondition::Type::NONCONFORM_INTERFACE>[iNonconformInterface];
        
        conditionsSet.valuesDirichlet = valueDirichletIt;
        conditionsSet.nValuesDirichlet = iValueDirichlet;
        conditionsSet.valuesNewman = valueNewmanIt;
        conditionsSet.nValuesNewman = iValueNewman;
        conditionsSet.nonconformInterfaces = nonconformInterfaceIt;
        conditionsSet.nNonconformInterfaces = iNonconformInterface;

        int* surfacesDimTags;
        GMSHProxy::model::getEntities(surfacesDimTags, nSurfaces, constants::DIMENSION_2D);
        conditionTypeByBoundaryTag = new Map::Element<int, const Boundary::ICondition*>[nSurfaces];

        physicalDimTagIt = physicalDimTags;
        physicalNameIt = physicalNames;

        const Boundary::ICondition* condition;
        for (unsigned int i = 0; i < nPhysicalGroups; ++i)
        {
            ++physicalDimTagIt;

            int* groupSurfacesTags;
            size_t nGroupSurfaces;
            GMSHProxy::model::getEntitiesForPhysicalGroup(constants::DIMENSION_2D, *physicalDimTagIt, groupSurfacesTags, nGroupSurfaces);

            char* physicalNameCharIt = *physicalNameIt;
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
                    
                    valueDirichletIt->value = value;
                    condition = valueDirichletIt;
                    ++valueDirichletIt;
                }
                else
                {
                    condition = &Boundary::ConditionsSet::dirichletFunctionCondition;
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
                    valueNewmanIt->value = value;
                    condition = valueNewmanIt;
                    ++valueNewmanIt;
                }
                else
                {
                    condition = &Boundary::ConditionsSet::newmanFunctionCondition;
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
                nonconformInterfaceIt->index = iNonconformInterface;
                condition = nonconformInterfaceIt;
                ++iNonconformInterface;
                ++nonconformInterfaceIt;
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

                    conditionsSet.stefan.value = value;
                    condition = &conditionsSet.stefan;
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

            nProcessedSurfaces += nGroupSurfaces;
            assignBoundaryConditionToGroupSurfaces(groupSurfacesTags, nGroupSurfaces, condition, conditionTypeByBoundaryTag, nSurfaces);

            GMSHProxy::free(groupSurfacesTags);
            ++physicalNameIt;
        }
        condition = &Boundary::ConditionsSet::homogeneousNewmanCondition;
        if (nProcessedSurfaces != nSurfaces)
        {
            const int* surfaceDimTagIt;
            for (unsigned int i = 0; i < nSurfaces; ++i)
            {
                ++surfaceDimTagIt;
                Map::Element<int, const Boundary::ICondition*>* mapElement = conditionTypeByBoundaryTag + (*surfaceDimTagIt) % nSurfaces;

                if (mapElement->tag != *surfaceDimTagIt)
                {
                    if (mapElement->tag == 0)
                    {
                        mapElement->tag = *surfaceDimTagIt;
                        mapElement->value = condition;
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
                            mapElement->value = condition;

                        }

                    }
                }

                ++surfaceDimTagIt;
            }
        }

        GMSHProxy::free(surfacesDimTags);

        return nProcessedSurfaces;
    }
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

    void getModel(const MaterialPhase materialPhases[2],
                  Volume*& volumes,
                  unsigned int& nVolumes,
                  Boundary* &boundaries,
                  Boundary::ConditionsSet& conditionsSet,
                  NonconformInterface*& nonconformInterfaces)
    {
        Map::Element<int, const Boundary::ICondition*>* conditionTypeByBoundaryTag;
        unsigned int nSurfaces;

        determineBoundariesConditions(conditionTypeByBoundaryTag, nSurfaces, conditionsSet);

        nonconformInterfaces = new NonconformInterface[conditionsSet.nNonconformInterfaces];
        bool* wasNonconfromInterfaceProcessed = new bool[conditionsSet.nNonconformInterfaces]{};

        int* volumesDimTags;
        int* boundariesDimTags;
        unsigned int* volumesBoundariesCounts;
        unsigned int nBoundaries;
        GMSHProxy::model::getVolumes(volumesDimTags, nVolumes, boundariesDimTags, volumesBoundariesCounts, nBoundaries);
        //GMSHProxy::model::getEntities(volumesDimTags, nVolumes, constants::DIMENSION_3D);

        Volume* volumeIt = new Volume[nVolumes];
        volumes = volumeIt;

        Map::Element<int, const MaterialPhase*>* materialPhaseByVolumeTag = new Map::Element<int, const MaterialPhase*>[nVolumes];
        determineVolumesMaterialPhases(materialPhases, materialPhaseByVolumeTag, nVolumes);

        Boundary* boundaryIt = new Boundary[nBoundaries];
        boundaries = boundaryIt;

        const int* boundaryDimTagIt = boundariesDimTags;
        const int* volumesDimTagsIt = volumesDimTags;
        const unsigned int* nVolumeBoundariesIt = volumesBoundariesCounts;

        for (unsigned int i = 0; i < nVolumes; ++i)
        {
            ++volumesDimTagsIt;
            volumeIt->tag = *volumesDimTagsIt;
            volumeIt->materialPhasePtr = Map::getExistedValue(*volumesDimTagsIt, materialPhaseByVolumeTag + volumeIt->tag % nVolumes);
            volumeIt->boundaries = boundaryIt;

            const unsigned int nVolumeBoundaries = *nVolumeBoundariesIt;
            for (unsigned int j = 0; j < nVolumeBoundaries; ++j)
            {
                ++boundaryDimTagIt;
                int boundaryTag = *boundaryDimTagIt;
                boundaryIt->tag = boundaryTag;
                boundaryIt->condition = Map::getExistedValue(boundaryTag, conditionTypeByBoundaryTag + boundaryTag % nSurfaces);
                boundaryIt->isPlane = GMSHProxy::model::isSurfacePlane(boundaryTag);

                if (boundaryIt->condition->type == Boundary::ICondition::Type::NONCONFORM_INTERFACE)
                {
                    Boundary::Condition<Boundary::ICondition::Type::NONCONFORM_INTERFACE>* nonconformCondition = (Boundary::Condition<Boundary::ICondition::Type::NONCONFORM_INTERFACE>*)boundaryIt->condition;
                    NonconformInterface* nonconformInterface = nonconformInterfaces + nonconformCondition->index;
                    if (wasNonconfromInterfaceProcessed[nonconformCondition->index])
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
                ++boundaryDimTagIt;
                ++boundaryIt;
            }
            ++nVolumeBoundariesIt;
        }

        delete[] materialPhaseByVolumeTag;
        delete[] conditionTypeByBoundaryTag;
        delete[] wasNonconfromInterfaceProcessed;

        GMSHProxy::free(volumesDimTags);
        GMSHProxy::free(boundariesDimTags);
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

