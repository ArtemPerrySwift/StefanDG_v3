#pragma once
#include "TetrahedronsAdjusments.h"
#include "CrossTetrahedronsAdjusments.h"
#include "ElementDGCalculator.h"
#include "FaceDGCalculator.h"
#include "LinearLagrangeBasis.h"
//#include "BoundariesConditions.h"
#include "DTGeometryKernel.h"
#include "GMSHProxy.h"
#include "Model.h"
#include "InterfaceSideElements.h"
#include "Solution.h"
#include "Container.h"
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
        Model model;
        model.initilizeByCurrentGMSHModel();

        ElementsAdjusmentsSet* volumesElementsAdjusmentsSets = new ElementsAdjusmentsSet[model.nVolumes];
        CrossElementsAdjusmentsSet* volumesElementsCrossAdjuesmentsSets = new CrossElementsAdjusmentsSet[model.nVolumes];
        CrossElementsAdjusmentsSet* nonconformInterfacesAdjuesmentsSets = new CrossElementsAdjusmentsSet[model.nNonconformInterfaces];

        double** firstApproximationIt = new double*[model.nVolumes];
        InterfaceSideElementsSet(*nonconformalInterfaceSideElementsSets)[2] = new InterfaceSideElementsSet[model.nNonconformInterfaces][2];
        InterfaceSideElementsSet(*conformalInterfaceSideElementsSets)[2] = new InterfaceSideElementsSet[model.nConformInterfaces][2];

        Container<size_t*, size_t>* sharedBoundariesFacesTagsSets = new Container<size_t*, size_t>[model.nSharedBoundaries];

        initilize(sharedBoundariesFacesTagsSets, model.nSharedBoundaries);

        const size_t nTIntervals = nTSteps - 1;
        double dt = (tMax - tMin) / nTIntervals;

        /*
        const Volume* volumeIt = volumes;
        ElementsAdjusmentsSet* elementsAdjusmentsIt = volumesElementsAdjusments;
        CrossElementsAdjusmentsSet* elementsCrossAdjuesmentsIt = volumesInteriorFacesAdjusments;
        double** firstApproximationIt = volumesInitialX;
        InterfaceSideElementsSet* interfaceSideElementsIt = interfacesSidesElements;
        */
        for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
        {

        }

        for (size_t i = 1; i < nTIntervals; ++i)
        {

        }


        model.clear();

        delete[] volumesElementsAdjusmentsSets;
        delete[] volumesElementsCrossAdjuesmentsSets;
        delete[] nonconformInterfacesAdjuesmentsSets;
        delete[] firstApproximationIt;
        delete[] nonconformalInterfaceSideElementsSets;
    }


private:

    void initilize(Container<size_t*, size_t>* sharedBoundariesFacesTagsSetIt, const unsigned int nPtrs)
    {
        for (unsigned int i = 0; i < nPtrs; ++i)
        {
            sharedBoundariesFacesTagsSetIt->elements = nullptr;
            ++sharedBoundariesFacesTagsSetIt;
        }
    }

    void preprocess(const Volume* volumeIt, const unsigned int nVolumes, size_t** volumeElementsTagsSetIt, size_t* nVolumeElementsIt, size_t& nElementsAdjusments, size_t& nDOFs)
    {
        size_t nElements = 0;
        for (unsigned int i = 0; i < nVolumes; ++i)
        {
            GMSHProxy::model::mesh::getTetrahedrons(*volumeElementsTagsSetIt, *nVolumeElementsIt, volumeIt->tag);
            nElements += *nVolumeElementsIt;
        }

        nElementsAdjusments = nElements * ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
        nDOFs = nElements * Basis::N_FUNCTIONS;
    }
    /*
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
    */


    void computeTetrahedronsAdjusments(const double (*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                       const double (*transpMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                       const double* determinantIt,
                                       const Coordinates* initPointIt,
                                       const size_t nElements,
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

