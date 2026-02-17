#pragma once
#include "BoundariesConditions.h"
#include "MaterialPhase.h"
#include "GeometryConstants.h"

namespace DTGeometryKernel
{
    const uint8_t N_CROSS_ELEMENTS = 2;

    void assignVolumesMaterialsPhases(const int volumesTags[], const size_t nVolumes, const MaterialPhase* materialPhaseIt, const uint8_t nMaterialsPhases, const MaterialPhase* volumesMaterialsPhasesPtrs[]);

    size_t countNonconformalVolumesInterfaces();

    void extractFaceNodesTags(const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES], const uint8_t faceLocalIndex, size_t* faceNodesTagsIt);

    void extractBoundaryTetrahedrons(const double tetrahedronsNativeJacobians[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
        const Coordinates tetrahedronsInitPoints[],
        const size_t* boundaryTetrahedronsIndexIt,
        const size_t nBoundaryTetrahedrons,
        double (*boundaryTetrahedronsNativeJacobiansIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
        Coordinates* boundaryTetrahedronInitPointsIt);

    void buildFragmentsModel(const int surafcesTags[N_CROSS_ELEMENTS],
        const char newModelName[],
        int*& fragmentTags,
        size_t& nFacesFragments,
        size_t(*&surfacesFragmentsTrianglesIndexes)[N_CROSS_ELEMENTS]);


    void determineFragmentsTetrahderonsIndexes(const size_t nFragments,
                                               const size_t(*fragmentsFacesIndexesIt)[N_CROSS_ELEMENTS],
                                               const size_t* const intefaceTetrahedronsIndexes[N_CROSS_ELEMENTS],
                                               size_t(*fragmentsTetrahedronsIndexesIt)[N_CROSS_ELEMENTS]);

    size_t createVolumeFacesEntities(const int volumeTag,
        const int boundariesTags[],
        const size_t nBoundaries,
        const BoundaryCondition boundariesConditions[],
        int& interiorFacesEntity,
        size_t& nInteriorFaces,
        int boundariesFacesEntytiesTags[],
        size_t nBoundariesFaces[],
        size_t(*&interiorFacesTetrahedronsIndexes)[2],
        uint8_t(*&interiorFacesLocalIndexes)[2],
        size_t* boundariesFacesTetrahedronsIndexes[],
        uint8_t* boundariesFacesLocalIndexes[]);

    void removeCreatedVolumeFacesEntities(const int* entitiesTagIt, const size_t nEntities);

    void determineBoundaryConditions(const int* boundaryTagIt, const size_t nBoundaries, BoundaryCondition* conditionTypeIt);

    void determineBoundaryConditions(const int* boundaryTagIt, const size_t nBoundaries, int* conditionTagIt, char** conditionNameIt, BoundaryCondition* conditionTypeIt);


    void computeInterfacesThermalConductivities(const MaterialPhase* volumeMaterialPhasePtrsIt[],
                                                const size_t(*intefacesVolumesIndexesIt)[2],
                                                const size_t nInterfaces,
                                                double* thermalConductivityIt);

    void get3DElementsByPoints(const Coordinates* pointIt, const uint8_t nPoints, size_t* elementTagIt, LocalCoordinates3D* localCoordinatesIt);
    
    void getModelNodesCoordinates(const int* nodeTagIt, size_t nNodes, Coordinates* coordinatesIt);

    void computeLocalJacobianMatrix(const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                    const double determinant,
                                    double nativeJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT]);

    uint8_t getClosestNodeIndex(const Coordinates& node, const size_t* nodeTagIt, const uint8_t nNodes);
};

