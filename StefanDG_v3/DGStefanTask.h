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
        unsigned int nVolumes, nNonconformInterafaces;

        getModel(materialPhases, volumes, nVolumes, nonconformInterfaces, nNonconformInterafaces);

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
                                                const Boundary::Condition boundaryCondition,
                                                Map::Element<int, Boundary::Condition>* conditionTypeByBoundaryTag,
                                                const unsigned int nMaxSurfaces)
    {
        for (unsigned int j = 0; j < nSurfaces; ++j)
        {
            const unsigned int hashIndex = *surfacesTagIt % nMaxSurfaces;
            Map::Element<int, Boundary::Condition>* mapElement = conditionTypeByBoundaryTag + hashIndex;

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

                mapElement->next = new Map::Element<int, Boundary::Condition>(*surfacesTagIt, boundaryCondition);
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

    unsigned int determineExplicitBoundariesConditions(Map::Element<int, Boundary::Condition>* conditionTypeByBoundaryTag,
                                                       const unsigned int nSurfaces,
                                                       unsigned int &nNonconformalGroups)
    {
        int* physicalDimTags;
        unsigned int nPhysicalGroups;
        GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalGroups, constants::DIMENSION_2D);

        nNonconformalGroups = 0;
        unsigned int nProcessedSurfaces = 0;
        const int* physicalDimTagIt = physicalDimTags;
        for (unsigned int i = 0; i < nPhysicalGroups; ++i)
        {
            ++physicalDimTagIt;

            int* groupSurfacesTags;
            size_t nGroupSurfaces;
            GMSHProxy::model::getEntitiesForPhysicalGroup(constants::DIMENSION_2D, *physicalDimTagIt, groupSurfacesTags, nGroupSurfaces);

            Boundary::Condition boundaryCondition;

            char* physicalName;
            GMSHProxy::model::getPhysicalName(constants::DIMENSION_2D, *physicalDimTagIt, physicalName);
            
            switch (*physicalName)
            {
            case 'D':
            {
                boundaryCondition.type = Boundary::Condition::Type::DIRICHLET;
                boundaryCondition.index = getConditionIndex(physicalName + 1);
                break;
            }
            case 'N':
            {
                boundaryCondition.type = Boundary::Condition::Type::NEWMAN;
                boundaryCondition.index = getConditionIndex(physicalName + 1);
                break;
            }
            case 'C':
            {
                if (nGroupSurfaces != 2)
                {
                    throw std::runtime_error("Each condition of nonconformity should include only 2 surfaces");
                }

                boundaryCondition.type = Boundary::Condition::Type::NONCONFORM_INTERFACE;
                boundaryCondition.index = nNonconformalGroups;
                ++nNonconformalGroups;
                break;
            }
            default:
            {
                throw std::runtime_error("Unknown type of condition");
            }
            }

            GMSHProxy::free(physicalName);

            nProcessedSurfaces += nGroupSurfaces;
            assignBoundaryConditionToGroupSurfaces(groupSurfacesTags, nGroupSurfaces, boundaryCondition, conditionTypeByBoundaryTag, nSurfaces);

            GMSHProxy::free(groupSurfacesTags);
            ++physicalDimTagIt;
        }

        GMSHProxy::free(physicalDimTags);

        return nProcessedSurfaces;
    }

    void determineImplicitBoundariesConditions(const int* surfaceDimTagIt,
                                               const unsigned int nSurfaces,
                                               Map::Element<int, Boundary::Condition>* conditionTypeByBoundaryTag)
    {
        unsigned int nConformInterafces = 0;
        for (unsigned int i = 0; i < nSurfaces; ++i)
        {
            ++surfaceDimTagIt;
            Map::Element<int, Boundary::Condition>* mapElement = conditionTypeByBoundaryTag + (*surfaceDimTagIt) % nSurfaces;

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
                        mapElement->value = { Boundary::Condition::Type::HOMOGENEOUS_NEWMAN, 0 };
                    }
                    else
                    {
                        mapElement->value = { Boundary::Condition::Type::HOMOGENEOUS_NEWMAN, nConformInterafces };
                        ++nConformInterafces;
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
                            mapElement->next = new Map::Element<int, Boundary::Condition>(*surfaceDimTagIt, { Boundary::Condition::Type::HOMOGENEOUS_NEWMAN, 0 });
                        }
                        else
                        {
                            mapElement->next = new Map::Element<int, Boundary::Condition>(*surfaceDimTagIt, { Boundary::Condition::Type::HOMOGENEOUS_NEWMAN, nConformInterafces });
                            ++nConformInterafces;
                        }
                        
                    }

                }
            }

            ++surfaceDimTagIt;
        }
    }

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

    void getModel(const MaterialPhase materialPhases[2], Volume*& volumes, unsigned int& nVolumes, NonconformInterface*& nonconformInterfaces, unsigned int& nNonconformInterafaces)
    {
        Map::Element<int, Boundary::Condition>* conditionTypeByBoundaryTag;
        unsigned int nSurfaces;

        determineBoundariesConditions(conditionTypeByBoundaryTag, nSurfaces, nNonconformInterafaces);

        nonconformInterfaces = new NonconformInterface[nNonconformInterafaces];
        bool* wasNonconfromInterfaceProcessed = new bool[nNonconformInterafaces]{};

        int* volumesDimTags;
        GMSHProxy::model::getEntities(volumesDimTags, nVolumes, constants::DIMENSION_3D);

        volumes = new Volume[nVolumes];
        Map::Element<int, const MaterialPhase*>* materialPhaseByVolumeTag = new Map::Element<int, const MaterialPhase*>[nVolumes];
        determineVolumesMaterialPhases(materialPhases, materialPhaseByVolumeTag, nVolumes);

        const int* volumesDimTagsIt = volumesDimTags;
        Volume* volumeIt = volumes;
        for (unsigned int i = 0; i < nVolumes; ++i)
        {
            ++volumesDimTagsIt;
            volumeIt->tag = *volumesDimTagsIt;
            volumeIt->materialPhasePtr = Map::getExistedValue(*volumesDimTagsIt, materialPhaseByVolumeTag + volumeIt->tag % nVolumes);

            int* boundariesDimTags;
            unsigned int nBoundaries;
            GMSHProxy::model::getBoundaries(constants::DIMENSION_3D, volumeIt->tag, boundariesDimTags, nBoundaries);
            volumeIt->boundaries = new Boundary[nBoundaries];

            const int* boundariesDimTagIt = boundariesDimTags;
            Boundary* boundaryIt = volumeIt->boundaries;
            for (unsigned int j = 0; j < nBoundaries; ++j)
            {
                ++boundariesDimTagIt;
                int boundaryTag = *boundariesDimTagIt;
                boundaryIt->tag = boundaryTag;
                boundaryIt->condition = Map::getExistedValue(boundaryTag, conditionTypeByBoundaryTag + boundaryTag % nSurfaces);
                boundaryIt->isPlane = GMSHProxy::model::isSurfacePlane(boundaryTag);

                if (boundaryIt->condition.type == Boundary::Condition::Type::NONCONFORM_INTERFACE)
                {
                    NonconformInterface* nonconformInterface = nonconformInterfaces + boundaryIt->condition.index;
                    if (wasNonconfromInterfaceProcessed[boundaryIt->condition.index])
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
                ++boundariesDimTagIt;
                ++boundaryIt;
            }

            GMSHProxy::free(boundariesDimTags);
        }

        delete[] materialPhaseByVolumeTag;
        delete[] conditionTypeByBoundaryTag;
        delete[] wasNonconfromInterfaceProcessed;
        GMSHProxy::free(volumesDimTags);
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

