#pragma once
#include "MaterialPhase.h"
#include "TetrahedronsAdjusments.h"
#include "CrossTetrahedronsAdjusments.h"
#include "ElementDGCalculator.h"
#include "FaceDGCalculator.h"
#include "LinearLagrangeBasis.h"
#include "BoundariesConditions.h"
#include "DTGeometryKernel.h"
#include "Volume.h"
#include "Interface.h"
#include "InterfaceSideElements.h"

class DGStefanTask
{
public:

    using Basis = LinearLagrangeBasis;

    static double dirichletCondition(const Coordinates& point);
    static double newmanCondition(const Coordinates& point);
    static double initialConditionLiquid(const Coordinates &point);
    static double initialConditionSolid(const Coordinates& point);
private:
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

    InterfaceSideElements* (*nonconformalInterfacesSidesElements)[2];
	CrossTetrahedronsAdjusments* nonconformalInterfacesFragmentsAdjusments;

    void processVolume(const Volume &volume,
                       const double dt,
                       TetrahedronsAdjusments& elementsAdjusments,
                       CrossTetrahedronsAdjusments& elementsCrossAdjuesments,
                       double* &firstApproximatuion)
    {
        size_t* tetrahedronsTags;
        size_t nTetrahedrons;

        size_t(*tetrahedronsNodesTags)[constants::tetrahedron::N_NODES];
        GMSHProxy::model::mesh::getTetrahedrons(tetrahedronsTags,
            nTetrahedrons,
            tetrahedronsNodesTags,
            volume.tag);

        elementsAdjusments.tetrahedronsTags = tetrahedronsTags;
        elementsAdjusments.nTetrahedrons = nTetrahedrons;

        double (*bilinearElementsAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[nTetrahedrons][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        double (*linearElementsAdjusments)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];
        double (*volumeInitialX)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];

        elementsAdjusments.bilinear = (double*)bilinearElementsAdjusments;
        elementsAdjusments.linear = (double*)linearElementsAdjusments;
        firstApproximatuion = (double*)volumeInitialX;

        double (*transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
        double* determinants;
        Coordinates* elementsInitPoints;

        GMSHProxy::model::mesh::getTetrahedronJacobians(transpJacobianMatrixes, determinants, elementsInitPoints, nTetrahedrons, volume.tag);

        double (*localJacobianMatrixes)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nTetrahedrons][LocalCoordinates3D::COUNT * Coordinates::COUNT];
        CoordinatesFunctions::computeLocalJacobianMatrixies(transpJacobianMatrixes, nTetrahedrons, determinants, localJacobianMatrixes);

        computeTetrahedronsAdjusments(localJacobianMatrixes,
            transpJacobianMatrixes,
            determinants,
            elementsInitPoints,
            nTetrahedrons,
            volume.materialPhasePtr->thermalConductivity,
            volume.materialPhasePtr->heatCapacity * volume.materialPhasePtr->density,
            dt,
            volume.materialPhasePtr->state == MaterialPhase::LIQUID ? initialConditionLiquid : initialConditionSolid,
            volumeInitialX,
            bilinearElementsAdjusments,
            linearElementsAdjusments);


        int interiorFacesEntityTag;
        int* boundariesFacesEntitiesTags = intBuffer;

        size_t(*interiorFacesTetrahedronsIndexes)[2];
        uint8_t(*interiorFacesLocalIndexes)[2];
        size_t interiorFaces;

        size_t** boundariesFacesTetrahedronsIndexes = (size_t**)ptrBuffer;
        uint8_t** boundariesFacesLocalIndexes = (uint8_t**)(ptrBuffer + volume.nBoundaries);
        size_t* nBoundariesFaces = sizesBuffer;

        DTGeometryKernel::createVolumeFacesEntities(volume.tag,
            volume.boundariesTags,
            volume.nBoundaries,
            volume.boundariesConditions,
            interiorFacesEntityTag,
            interiorFaces,
            boundariesFacesEntitiesTags,
            nBoundariesFaces,
            interiorFacesTetrahedronsIndexes,
            interiorFacesLocalIndexes,
            boundariesFacesTetrahedronsIndexes,
            boundariesFacesLocalIndexes);

        elementsCrossAdjuesments.tetrahedronsIndexes = interiorFacesTetrahedronsIndexes;
        elementsCrossAdjuesments.nCrossElements = interiorFaces;

        double (*interiorFacesAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[interiorFaces][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        elementsCrossAdjuesments.bilinear = (double*)interiorFacesAdjusments;

        computeInteriorFacesAdjusments(interiorFacesEntityTag,
                                       interiorFacesTetrahedronsIndexes,
                                       interiorFacesLocalIndexes,
                                       volume.materialPhasePtr->thermalConductivity,
                                       localJacobianMatrixes,
                                       bilinearElementsAdjusments,
                                       interiorFacesAdjusments);

        delete[] interiorFacesLocalIndexes;

        const int* boundariesTagIt = volume.boundariesTags;
        const int* boundaryFacesEntityTagIt = boundariesFacesEntitiesTags;

        const size_t* const* boundaryFacesTetrahedronsIndexesIt = boundariesFacesTetrahedronsIndexes;
        const uint8_t* const* boundaryFacesLocalIndexesIt = boundariesFacesLocalIndexes;

        const BoundaryCondition* boundaryConditionIt = volume.boundariesConditions;

        for (size_t j = 0; j < volume.nBoundaries; ++j)
        {
            switch (*boundaryConditionIt)
            {
            case BoundaryCondition::DIRICHLET:
            {
                char* surfaceType;
                GMSHProxy::model::getSurafceType(*boundariesTagIt, surfaceType);

                if (strcmp(surfaceType, "Plane"))
                {
                    computeDirichletFacesAdjusments(*boundaryFacesEntityTagIt,
                                                    *boundaryFacesTetrahedronsIndexesIt,
                                                    *boundaryFacesLocalIndexesIt,
                                                    localJacobianMatrixes,
                                                    dirichletCondition,
                                                    volume.materialPhasePtr->thermalConductivity,
                                                    bilinearElementsAdjusments,
                                                    linearElementsAdjusments);
                }
                else
                {
                    size_t faceNodesIndexes[constants::tetrahedron::N_FACES * constants::triangle::N_NODES];
                    DTGeometryKernel::extractFaceNodesTags(tetrahedronsNodesTags[**boundaryFacesTetrahedronsIndexesIt], **boundaryFacesLocalIndexesIt, faceNodesIndexes);

                    Coordinates faceNodes[constants::triangle::N_NODES];
                    GMSHProxy::model::mesh::get3Nodes(faceNodesIndexes, faceNodes);

                    Coordinates normal;
                    CoordinatesFunctions::computeNormal(faceNodes, normal);

                    computeDirichletFacesAdjusments(*boundaryFacesEntityTagIt,
                                                    *boundaryFacesTetrahedronsIndexesIt,
                                                    *boundaryFacesLocalIndexesIt,
                                                    normal,
                                                    localJacobianMatrixes,
                                                    dirichletCondition,
                                                    volume.materialPhasePtr->thermalConductivity,
                                                    bilinearElementsAdjusments,
                                                    linearElementsAdjusments);

                }

                GMSHProxy::free(surfaceType);

                break;

            }

            case BoundaryCondition::NEWMAN:
            {
                computeNewmanFacesAdjusments(*boundaryFacesEntityTagIt,
                    *boundaryFacesTetrahedronsIndexesIt,
                    *boundaryFacesLocalIndexesIt,
                    newmanCondition,
                    linearElementsAdjusments);
                break;
            }

            case BoundaryCondition::NONCONFORM:
            {
                double (*boundaryTetrahedronsLocalJacobians)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nBoundariesFaces[j]][LocalCoordinates3D::COUNT * Coordinates::COUNT];
                Coordinates* boundaryTetrahedronsInitPoints = new Coordinates[nBoundariesFaces[j]];
                DTGeometryKernel::extractBoundaryTetrahedrons(localJacobianMatrixes, elementsInitPoints, *boundaryFacesTetrahedronsIndexesIt, nBoundariesFaces[j], boundaryTetrahedronsLocalJacobians, boundaryTetrahedronsInitPoints);

                InterfaceSideBuffer* interfaceSideBuffer;

                auto intefaceIndexIt = interfaceIndexByPhysicalTag.find(boundariesConditionsTags[j]);
                if (intefaceIndexIt == interfaceIndexByPhysicalTag.end())
                {
                    interfaceIndexByPhysicalTag.emplace(boundariesConditionsTags[j], interfaceIndex);

                    interfaceSideBuffer = *volumesNonconformalInterfaceBoundaryIt;

                    size_t faceNodesTags[constants::tetrahedron::N_FACES * constants::triangle::N_NODES];
                    DTGeometryKernel::extractFaceNodesTags(tetrahedronsNodesTags[**boundaryFacesTetrahedronsIndexesIt], **boundaryFacesLocalIndexesIt, faceNodesTags);
                    Coordinates faceNodes[constants::triangle::N_NODES];
                    GMSHProxy::model::mesh::get3Nodes(faceNodesTags, faceNodes);
                    CoordinatesFunctions::computeNormal(faceNodes, interfaceSideBuffer->normal);


                    ++volumesNonconformalInterfaceBoundaryIt;
                    ++nonconformalIntefaceSidesTagsIt;
                    ++nonconformalIntefaceTetrahedronsIndexesIt;
                    ++interfaceIndex;
                }
                else
                {
                    InterfaceSideBuffer* volumeInterfaceBoundaryOtherSide = volumesNonconformalInterfaceBoundaries[intefaceIndexIt->second];
                    interfaceSideBuffer = volumeInterfaceBoundaryOtherSide;
                    ++interfaceSideBuffer;

                    nonconformalIntefacesVolumesIndexes[intefaceIndexIt->second][1] = i;
                    nonconformalIntefacesSidesTags[intefaceIndexIt->second][1] = *boundariesTagIt;
                    nonconformalIntefacesTetrahedronsIndexes[intefaceIndexIt->second][1] = *boundaryFacesTetrahedronsIndexesIt;

                    interfaceSideBuffer->normal.x = -volumeInterfaceBoundaryOtherSide->normal.x;
                    interfaceSideBuffer->normal.y = -volumeInterfaceBoundaryOtherSide->normal.y;
                    interfaceSideBuffer->normal.z = -volumeInterfaceBoundaryOtherSide->normal.z;
                }

                interfaceSideBuffer->tetrahedronsLocalJacobians = boundaryTetrahedronsLocalJacobians;
                interfaceSideBuffer->tetrahedronsInitPoints = boundaryTetrahedronsInitPoints;
                interfaceSideBuffer->facesLocalIndexes = *boundaryFacesLocalIndexesIt;


                computeVolumesInterafaceAdjusments(*boundariesTagIt,
                    *boundaryFacesTetrahedronsIndexesIt,
                    *boundaryFacesLocalIndexesIt,
                    interfaceSideBuffer->normal,
                    penaltyCoeff,
                    volumeMaterialPhasePtr->thermalConductivity,
                    localJacobianMatrixes,
                    bilinearElementsAdjusments);

                *boundaryFacesTetrahedronsIndexesIt = nullptr;
                *boundaryFacesLocalIndexesIt = nullptr;

                break;
            }
            }

            delete[] * boundaryFacesTetrahedronsIndexesIt;
            delete[] * boundaryFacesLocalIndexesIt;

            ++boundaryConditionIt;
            ++boundariesTagIt;
            ++boundaryFacesEntityTagIt;
            ++boundaryFacesTetrahedronsIndexesIt;
            ++boundaryFacesLocalIndexesIt;

        }

        DTGeometryKernel::removeCreatedVolumeFacesEntities(boundariesFacesEntitiesTags, nBoundaries);
    }


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


};

