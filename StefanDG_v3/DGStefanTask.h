#pragma once
#include "MaterialPhase.h"
#include "TetrahedronsAdjusments.h"
#include "CrossTetrahedronsAdjusments.h"
#include "ElementDGCalculator.h"
#include "FaceDGCalculator.h"
#include "LinearLagrangeBasis.h"
#include "BoundariesConditions.h"
#include "DTGeometryKernel.h"
#include "GMSHProxy.h"
#include "Volume.h"
#include "Interface.h"
#include "InterfaceSideElements.h"
#include "Solution.h"
#include "Map.h"

#include <stdexcept>


class DGStefanTask
{
public:

    using Basis = LinearLagrangeBasis;

    static double dirichletCondition(const Coordinates& point);
    static double newmanCondition(const Coordinates& point);
    static double initialConditionLiquid(const Coordinates &point);
    static double initialConditionSolid(const Coordinates& point);

    void solve(const MaterialPhase* materialPhases,
               const double stefanTemperature,
               const double penalty,
               const double tMin,
               const double tMax,
               const size_t tSteps,
               Solution* solutionIt)
    {
        Volume* volumes;
        unsigned int nVolumes;


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

    const MaterialPhase* getMaterialPhaseForVolume(const int tag, const MaterialPhase materialPhases[2], int materialPhasesTags[2])
    {
        int* entityPGsTags;
        size_t nEntityPGsTags;
        GMSHProxy::model::getPhysicalGroupsForEntity(constants::DIMENSION_3D, tag, entityPGsTags, nEntityPGsTags);

        const int* entityPGsTagIt = entityPGsTags;
        for (size_t j = 0; j < nEntityPGsTags; ++j)
        {
            if (materialPhasesTags[0] == *entityPGsTagIt)
            {
                return materialPhases;
            }
            else
            {
                if (materialPhasesTags[1] == *entityPGsTagIt)
                {
                    return materialPhases + 1;
                }
            }

            ++entityPGsTagIt;
        }

        return nullptr;
    }

    void getVolumeBoundaries(const int tag, int* &boundariesTags, unsigned int &nBoundaries)
    {
        int* boundariesDimTags;
        size_t nBoundariesDimTags;
        GMSHProxy::model::getBoundaries(constants::DIMENSION_3D, tag, boundariesDimTags, nBoundariesDimTags);

        nBoundaries = nBoundariesDimTags / 2;
        boundariesTags = new int[nBoundaries];

        int* boundaryDimTagIt = boundariesDimTags;
        int* boundaryTagIt = boundariesTags;
        for (unsigned int i = 0; i < nBoundaries; ++i)
        {
            ++boundaryDimTagIt;
            *boundaryTagIt = *boundaryDimTagIt;

            ++boundaryDimTagIt;
            ++boundaryTagIt;
        }
    }

    void getVolumes(const MaterialPhase materialPhases[2], Volume*& volumes, unsigned int& nVolumes)
    {
        int* volumesDimTags;
        size_t nVolumeDimTags;
        GMSHProxy::model::getEntities(volumesDimTags, nVolumeDimTags, constants::DIMENSION_3D);

        nVolumes = nVolumeDimTags / 2;
        volumes = new Volume[nVolumes];

        int materialPhasesTags[2];
        getMaterialPhasesTags(materialPhases, materialPhasesTags);

        const int* volumesDimTagsIt = volumesDimTags;
        Volume* volumeIt = volumes;
        for (unsigned int i = 0; i < nVolumes; ++i)
        {
            ++volumesDimTagsIt;
            volumeIt->tag = *volumesDimTagsIt;
            volumeIt->materialPhasePtr = getMaterialPhaseForVolume(volumeIt->tag, materialPhases, materialPhasesTags);

            int* boundariesTags;
            getVolumeBoundaries(volumeIt->tag, boundariesTags, volumeIt->nBoundaries);
            volumeIt->boundariesTags = boundariesTags;

            BoundaryCondition* boundariesConditions = new BoundaryCondition[volumeIt->nBoundaries];
            DTGeometryKernel::determineBoundaryConditions(volumeIt->boundariesTags, volumeIt->nBoundaries, boundariesConditions);
            volumeIt->boundariesConditions = boundariesConditions;
        }
    }

    unsigned int assignBoundariesTagsAndConditionTypes(const int* physicalDimTagIt,
                                                       const unsigned int nPhysicalGroups,
                                                       Map::MapElement<int, BoundaryCondition>* conditionTypeByBoundaryTag,
                                                       const unsigned int nMaxSurfaces)
    {
        unsigned int nSurfaces = 0;
        for (unsigned int i = 0; i < nPhysicalGroups; ++i)
        {
            ++physicalDimTagIt;

            char* physicalName;
            GMSHProxy::model::getPhysicalName(constants::DIMENSION_2D, *physicalDimTagIt, physicalName);
            BoundaryCondition boundaryCondition;
            switch (*physicalName)
            {
            case 'D':
            {
                boundaryCondition = BoundaryCondition::DIRICHLET;
                break;
            }
            case 'N':
            {
                boundaryCondition = BoundaryCondition::NEWMAN;
                break;
            }
            case 'C':
            {
                boundaryCondition = BoundaryCondition::NONCONFORM_INTERAFACE;
                break;
            }
            default:
            {
                throw std::runtime_error("Unknown type of condition");
            }
            }
            GMSHProxy::free(physicalName);

            int* groupSurfacesTags;
            size_t nGroupSurfacesTags;
            GMSHProxy::model::getEntitiesForPhysicalGroup(constants::DIMENSION_2D, *physicalDimTagIt, groupSurfacesTags, nGroupSurfacesTags);

            nSurfaces += nGroupSurfacesTags;

            const int* groupSurfacesTagIt = groupSurfacesTags;
            for (unsigned int j = 0; j < nGroupSurfacesTags; ++j)
            {
                Map::addElement<int, BoundaryCondition>(*groupSurfacesTagIt, boundaryCondition, conditionTypeByBoundaryTag, nMaxSurfaces);
                ++groupSurfacesTagIt;
            }

            GMSHProxy::free(groupSurfacesTags);
            ++physicalDimTagIt;
        }

        return nSurfaces;
    }

    void getInterafces(const Volume* volumes, const unsigned int nVolumes, Interface*& nonconformInterfaces, unsigned int& nInterafaces)
    {
        int* surfacesDimTags;
        size_t nSurfacesDimTags;
        GMSHProxy::model::getEntities(surfacesDimTags, nSurfacesDimTags, constants::DIMENSION_2D);

        unsigned int nSurfaces = nSurfacesDimTags / 2;
        int* physicalDimTags;
        size_t nPhysicalDimTags;
        GMSHProxy::model::getPhysicalGroups(physicalDimTags, nPhysicalDimTags, constants::DIMENSION_2D);
        unsigned int nPhysicalGroups2D = nPhysicalDimTags / 2;

        Map::MapElement<int, BoundaryCondition>* conditionTypeByBoundaryTag = new Map::MapElement<int, BoundaryCondition>[nSurfaces];

        unsigned int nProcessedSurfaces = assignBoundariesTagsAndConditionTypes(physicalDimTags, nPhysicalGroups2D, conditionTypeByBoundaryTag, nSurfaces);

        if (nProcessedSurfaces != nSurfaces)
        {
            const int* surfaceDimTagIt = surfacesDimTags;
            for (unsigned int i = 0; i < nSurfaces; ++i)
            {
                ++surfaceDimTagIt;
                const Map::MapElement<int, BoundaryCondition>* conditionTypeByBoundaryTagIt = Map::tryElementAdding<int, BoundaryCondition>(*surfaceDimTagIt, BoundaryCondition::HOMOGENEOUS_NEWMAN, conditionTypeByBoundaryTag, nSurfaces);
            }
        }
        //std::unordered_map<int, BoundaryCondition> conditionTypeByBoundaryTag;
        //conditionTypeByBoundaryTag.reserve(nSurfaces);

        //std::pair<int, BoundaryCondition>* pairsTagCondition = new std::pair<int, BoundaryCondition>[nSurfaces];
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

	TetrahedronsAdjusments* volumesElementsAdjusments;
	CrossTetrahedronsAdjusments* volumesInteriorFacesAdjusments;
	double** volumesInitialX;

    Interface* nonconformalInterfaces;
	size_t nNonconformalInterfaces;

    InterfaceSideElements* interfacesSidesElements;
	CrossTetrahedronsAdjusments* nonconformalInterfacesFragmentsAdjusments;


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

