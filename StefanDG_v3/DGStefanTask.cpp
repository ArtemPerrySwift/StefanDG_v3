#include "DGStefanTask.h"
#include "CoordinatesFunctions.h"
#include "GMSHProxy.h"
#include "FacesSet.h"
#include <algorithm>

template<size_t N_FUNCTIONS>
static void computeTetrahedronAdjusments(const double* massMatrixElementIt, const double gamma, const double* stiffnessMatrixElementIt, const double lambda, const double* trackPowerVectorElementIt, const double determinant, double* bilinearAdjusmentIt, double* linearAdjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *linearAdjusmentIt = -determinant * gamma * (*trackPowerVectorElementIt);
        for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
        {
            *bilinearAdjusmentIt = determinant * (*massMatrixElementIt - *stiffnessMatrixElementIt);

            ++massMatrixElementIt;
            ++stiffnessMatrixElementIt;
            ++bilinearAdjusmentIt;
        }

        ++trackPowerVectorElementIt;
        ++linearAdjusmentIt;
    }
}

void computeTetrahedronsAdjusments(const double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
    const double(*transpMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
    const double* determinantIt,
    const Coordinates* initPointIt,
    const size_t nTetrahedrons,
    const double lambda,
    const double gamma,
    const double dt,
    double(*initialCondition)(const Coordinates& point),
    void* memoryBuffer,
    double(*initialXIt)[Basis::N_FUNCTIONS],
    double(*bilinearTetrahedronsAdjusmentsIt)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
    double(*linearTetrahedronsAdjusmentsIt)[Basis::N_FUNCTIONS])
{
    const double* templateMassMatrix = ElementDGCalculator<Basis>::getMassMatrix();
    double* stiffnessMatrix = (double*)memoryBuffer;
    double* localVector = stiffnessMatrix + ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double* initlialConditionValues = localVector + Basis::N_FUNCTIONS;

    const LocalCoordinates3D* localGradients = ElementDGCalculator<Basis>::getIntegrationLocalGradients();
    Coordinates* gradients = (Coordinates*)(initlialConditionValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);
    Coordinates* integrationPoints = gradients + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

    double massCoeff = gamma / dt;

    const LLt* massSLAESolverPtr = ElementDGCalculator<Basis>::getMassSLAESolverPtr();

    for (size_t tetrahedronIndex = 0; tetrahedronIndex < nTetrahedrons; tetrahedronIndex++)
    {
        CoordinatesFunctions::translate(*localJacobianMatrixIt, localGradients, NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
        ElementDGCalculator<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

        CoordinatesFunctions::translate(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints,
            NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps,
            *initPointIt,
            *transpMatrixIt,
            integrationPoints);

        std::transform(integrationPoints, integrationPointsEndIt, initlialConditionValues, initialCondition);
        ElementDGCalculator<Basis>::computePowerVector(initlialConditionValues, localVector);
        massSLAESolverPtr->solve(localVector, *initialXIt);

        computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                                                         massCoeff,
                                                         stiffnessMatrix,
                                                         lambda,
                                                         localVector,
                                                         *determinantIt,
                                                         *bilinearTetrahedronsAdjusmentsIt,
                                                         *linearTetrahedronsAdjusmentsIt);

        ++initialXIt;
        ++initPointIt;
        ++determinantIt;
        ++transpMatrixIt;
        ++localJacobianMatrixIt;

        ++bilinearTetrahedronsAdjusmentsIt;
        ++linearTetrahedronsAdjusmentsIt;
    }
}

template<size_t N_FUNCTIONS>
static void addInteriorFaceAdjusmets(const double* flowMatrixElementIt, const double* massMatrixElementIt, const double penaltyCoeff, const double detLambdaCoeff, double* adjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *adjusmentIt = detLambdaCoeff * (*flowMatrixElementIt - penaltyCoeff * (*massMatrixElementIt));
        const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
        const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
        for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
        {
            ++massMatrixElementIt;
            ++ijFlowMatrixElementPtr;
            jiFlowMatrixElementPtr += N_FUNCTIONS;

            *adjusmentIt = detLambdaCoeff * ((*ijFlowMatrixElementPtr + *jiFlowMatrixElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
        }

        flowMatrixElementIt += N_FUNCTIONS + 1;
    }
}

template<size_t N_FUNCTIONS>
void computeInteriorFaceCrossAdjusmets(const double* crossFlowMatrix1ElementIt,
                                       const double* crossFlowMatrix2ElementIt,
                                       const double* massMatrixElementIt,
                                       const double penaltyCoeff,
                                       const double detLambdaCoeff,
                                       double* adjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *adjusmentIt = detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *crossFlowMatrix2ElementIt) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
        const double* jiFlowMatrix2ElementPtr = crossFlowMatrix2ElementIt;
        for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
        {
            ++massMatrixElementIt;
            ++crossFlowMatrix1ElementIt;
            jiFlowMatrix2ElementPtr += N_FUNCTIONS;

            *adjusmentIt = detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *jiFlowMatrix2ElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
        }

        ++crossFlowMatrix2ElementIt;
    }
}


void computeInteriorFacesAdjusments(const int interiorFacesEntityTag,
                                    const FacesSet &interiorFacesSet,
                                    const double lambda,
                                    const double penalty,
                                    const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                    void* memoryBuffer,
                                    double(bilinearTetrahedronsAdjusments)[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                    double(*bilinearInterfaceAdjusmentsIt)[FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS])
{
    const double* const* massMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceDGCalculator<Basis>::getCrossMassMatrixies();

    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();

    double* flowMatrix = (double*)memoryBuffer;
    double* crossFlowMatrix1 = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* crossFlowMatrix2 = crossFlowMatrix1 + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* normalDerivatives = crossFlowMatrix2 + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS);

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, interiorFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    size_t(*facesNodesTags)[constants::triangle::N_NODES];
    Coordinates(*facesNodes)[constants::triangle::N_NODES];
    GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, interiorFacesEntityTag);

    Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
    GMSHProxy::free(facesNodesTags);

    const size_t* faceElementIndexIt = interiorFacesSet.elementIndexes;
    const uint8_t* faceLocalIndexIt = interiorFacesSet.localindexes;

    size_t nFaces = interiorFacesSet.count;
    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;
        double commonMultiplier = *determinantIt * lambda;

        Coordinates normal;
        CoordinatesFunctions::computeNormal(*faceNodesIt, normal);

        CoordinatesFunctions::translate(localJacobianMatrix[faceElementIndexIt[0]],
            gradientsByFace[faceLocalIndexIt[0]],
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[0], normalDerivatives, flowMatrix);
        addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[faceLocalIndexIt[0]], penalty, commonMultiplier, bilinearTetrahedronsAdjusments[faceElementIndexIt[0]]);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[1], normalDerivatives, crossFlowMatrix1);

        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;

        CoordinatesFunctions::translate(localJacobianMatrix[faceElementIndexIt[1]],
            gradientsByFace[faceLocalIndexIt[1]],
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[1], normalDerivatives, flowMatrix);
        addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[faceLocalIndexIt[1]], penalty, commonMultiplier, bilinearTetrahedronsAdjusments[faceElementIndexIt[1]]);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[0], normalDerivatives, crossFlowMatrix2);

        computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                                                              crossFlowMatrix1,
                                                              crossMassMatrixies[faceLocalIndexIt[0]][faceLocalIndexIt[1]],
                                                              facePenalty,
                                                              commonMultiplier,
                                                              *bilinearInterfaceAdjusmentsIt);

        ++faceElementIndexIt;
        ++faceElementIndexIt;
        ++faceLocalIndexIt;
        ++faceLocalIndexIt;
        ++bilinearInterfaceAdjusmentsIt;
        ++determinantIt;
        ++faceNodesIt;
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);
    GMSHProxy::free(facesNodes);

}

template<size_t N_FUNCTIONS>
static void addDirichletFaceBilinearAdjusmets(const double* flowMatrixElementIt, const double* massMatrixElementIt, const double penaltyCoeff, const double detLambdaCoeff, double* adjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *adjusmentIt = detLambdaCoeff * (*flowMatrixElementIt + *flowMatrixElementIt - penaltyCoeff * (*massMatrixElementIt));
        const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
        const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
        for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
        {
            ++massMatrixElementIt;
            ++ijFlowMatrixElementPtr;
            jiFlowMatrixElementPtr += N_FUNCTIONS;

            *adjusmentIt += detLambdaCoeff * (*ijFlowMatrixElementPtr + *jiFlowMatrixElementPtr - penaltyCoeff * (*massMatrixElementIt));
        }

        flowMatrixElementIt += N_FUNCTIONS + 1;
    }
}

template<size_t N_FUNCTIONS>
static void addDirichletFaceLinearAdjusments(const double* powerVectorIt, const double* flowVectorElementIt, const double penaltyCoeff, const double commonCoeff, double* adjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *adjusmentIt += commonCoeff * (penaltyCoeff * (*powerVectorIt) - *flowVectorElementIt);

        ++powerVectorIt;
        ++flowVectorElementIt;
        ++adjusmentIt;
    }
}


void computeDirichletFacesAdjusments(const FacesSet& boundaryFacesSet,
                                     const Coordinates outwardNormal,
                                     const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                     double(*const dirichletCondition)(const Coordinates&),
                                     const double lambda,
                                     const double penalty,
                                     void* memoryBuffer,
                                     double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                     double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    double* flowMatrix = (double*)memoryBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* powerVector = flowVector + Basis::N_FUNCTIONS;

    double* dirichletValues = powerVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = (Coordinates*)(normalDerivatives + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, boundaryFacesSet.tag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    const size_t* boundaryFaceTetrahedronIndexIt = boundaryFacesSet.elementIndexes;
    const uint8_t* boundaryFacesLocalIndexesIt = boundaryFacesSet.localindexes;

    size_t nFaces = boundaryFacesSet.count;
    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double determinant = *determinantIt;
        double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
        double commonMultiplier = determinant * lambda;


        CoordinatesFunctions::translate(localJacobianMatrix[*boundaryFaceTetrahedronIndexIt],
                                        gradientsByFace[*boundaryFacesLocalIndexesIt],
                                        NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                                        gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, outwardNormal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(*boundaryFacesLocalIndexesIt, normalDerivatives, flowMatrix);
        addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix,
                                                              templatePenaltyMatrixByFace[*boundaryFacesLocalIndexesIt],
                                                              facePenalty,
                                                              commonMultiplier,
                                                              bilinearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

        CoordinatesFunctions::translate(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                                        NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps,
                                        *initPointIt,
                                        *transpJacobianMatrixIt,
                                        integrationPoints);

        std::transform(integrationPoints, integrationPointsEndIt, dirichletValues, dirichletCondition);

        FaceDGCalculator<Basis>::computeFlowVector(normalDerivatives, dirichletValues, flowVector);
        FaceDGCalculator<Basis>::computePowerVector(*boundaryFacesLocalIndexesIt, dirichletValues, powerVector);

        addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(powerVector, flowVector, facePenalty, commonMultiplier, linearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

        ++boundaryFaceTetrahedronIndexIt;
        ++boundaryFacesLocalIndexesIt;
        ++determinantIt;
        ++initPointIt;
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);
}


void computeDirichletFacesAdjusments(const int boundaryTag,
                                     const FacesSet& boundaryFacesSet,
                                     const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                     double(* const dirichletCondition)(const Coordinates&),
                                     const double lambda,
                                     const double penalty,
                                     void* memoryBuffer,
                                     double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                     double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    double* flowMatrix = (double*)memoryBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* powerVector = flowVector + Basis::N_FUNCTIONS;

    double* dirichletValues = powerVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = (Coordinates*)(normalDerivatives + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, boundaryFacesSet.tag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;


    const size_t* boundaryFaceTetrahedronIndexIt = boundaryFacesSet.elementIndexes;
    const uint8_t* boundaryFacesLocalIndexesIt = boundaryFacesSet.localindexes;

    const size_t nFaces = boundaryFacesSet.count;

    if (GMSHProxy::model::isSurfacePlane(boundaryTag))
    {

    }
    else
    {
        size_t(*facesNodesTags)[constants::triangle::N_NODES];
        Coordinates(*facesNodes)[constants::triangle::N_NODES];
        GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, boundaryFacesSet.tag);

        Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
        GMSHProxy::free(facesNodesTags);

        for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
        {
            double determinant = *determinantIt;
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
            double commonMultiplier = determinant * lambda;

            CoordinatesFunctions::translate(localJacobianMatrix[*boundaryFaceTetrahedronIndexIt],
                gradientsByFace[*boundaryFacesLocalIndexesIt],
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                gradients);

            Coordinates normal;
            CoordinatesFunctions::computeNormal(*faceNodesIt, normal);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

            FaceDGCalculator<Basis>::computeFlowMatrix(*boundaryFacesLocalIndexesIt, normalDerivatives, flowMatrix);
            addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix,
                templatePenaltyMatrixByFace[*boundaryFacesLocalIndexesIt],
                facePenalty,
                commonMultiplier,
                bilinearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

            CoordinatesFunctions::translate(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                *initPointIt,
                *transpJacobianMatrixIt,
                integrationPoints);

            std::transform(integrationPoints, integrationPointsEndIt, dirichletValues, dirichletCondition);

            FaceDGCalculator<Basis>::computeFlowVector(normalDerivatives, dirichletValues, flowVector);
            FaceDGCalculator<Basis>::computePowerVector(*boundaryFacesLocalIndexesIt, dirichletValues, powerVector);

            addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(powerVector, flowVector, facePenalty, commonMultiplier, linearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

            ++boundaryFaceTetrahedronIndexIt;
            ++boundaryFacesLocalIndexesIt;
            ++determinantIt;
            ++initPointIt;
            ++faceNodesIt;
        }

        GMSHProxy::free(facesNodes);
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);
}

void computeDirichletFacesAdjusments(const FacesSet& boundaryFacesSet,
                                     const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                     const double condition,
                                     const double lambda,
                                     const double penalty,
                                     void* memoryBuffer,
                                     double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                     double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    const double (*massVectorByFace)[Basis::N_FUNCTIONS] = FaceDGCalculator<Basis>::getMassVectors();
    double* flowMatrix = (double*)memoryBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double* dirichletValues = flowVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = (Coordinates*)(normalDerivatives + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, boundaryFacesSet.tag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    size_t(*facesNodesTags)[constants::triangle::N_NODES];
    Coordinates(*facesNodes)[constants::triangle::N_NODES];
    GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, boundaryFacesSet.tag);

    Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
    GMSHProxy::free(facesNodesTags);

    const size_t* boundaryFaceTetrahedronIndexIt = boundaryFacesSet.elementIndexes;
    const uint8_t* boundaryFacesLocalIndexesIt = boundaryFacesSet.localindexes;

    const size_t nFaces = boundaryFacesSet.count;
    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double determinant = *determinantIt;
        double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
        double commonMultiplier = determinant * lambda * condition;

        size_t tetrahedronIndex = *boundaryFaceTetrahedronIndexIt;
        uint8_t faceLocalIndex = *boundaryFacesLocalIndexesIt;

        CoordinatesFunctions::translate(localJacobianMatrix[*boundaryFaceTetrahedronIndexIt],
            gradientsByFace[*boundaryFacesLocalIndexesIt],
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            gradients);

        Coordinates normal;
        CoordinatesFunctions::computeNormal(*faceNodesIt, normal);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(*boundaryFacesLocalIndexesIt, normalDerivatives, flowMatrix);
        addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix,
            templatePenaltyMatrixByFace[*boundaryFacesLocalIndexesIt],
            facePenalty,
            commonMultiplier,
            bilinearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

        CoordinatesFunctions::translate(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            *initPointIt,
            *transpJacobianMatrixIt,
            integrationPoints);

        FaceDGCalculator<Basis>::computeFlowVector(normalDerivatives, flowVector);

        addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(massVectorByFace[*boundaryFacesLocalIndexesIt], flowVector, facePenalty, commonMultiplier, linearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

        ++boundaryFaceTetrahedronIndexIt;
        ++boundaryFacesLocalIndexesIt;
        ++determinantIt;
        ++initPointIt;
        ++faceNodesIt;
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);
    GMSHProxy::free(facesNodes);
}

template<size_t N_FUNCTIONS>
static void addNewmanFaceLinearAdjusments(const double* newmanVectorElementIt, const double determinant, double* adjusmentIt)
{
    for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
    {
        *adjusmentIt += determinant * (*newmanVectorElementIt);

        ++adjusmentIt;
        ++newmanVectorElementIt;
    }
}


void computeNewmanFacesAdjusments(const FacesSet& boundaryFacesSet,
                                  double(* const newmanCondition)(const Coordinates&),
                                  void* memoryBuffer,
                                  double linearTetrahedronsAdjusments[][LinearLagrangeBasis::N_FUNCTIONS])
{
    double* newmanValues = (double*)memoryBuffer;
    double* newmanVector = newmanValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = (Coordinates*)(newmanVector + LinearLagrangeBasis::N_FUNCTIONS);
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, boundaryFacesSet.tag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    const size_t* boundaryFaceTetrahedronIndexIt = boundaryFacesSet.elementIndexes;
    const uint8_t* boundaryFacesLocalIndexesIt = boundaryFacesSet.localindexes;

    const size_t nFaces = boundaryFacesSet.count;
    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {

        CoordinatesFunctions::translate(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, *initPointIt, *transpJacobianMatrixIt, integrationPoints);
        std::transform(integrationPoints, integrationPointsEndIt, newmanValues, newmanCondition);

        FaceDGCalculator<Basis>::computePowerVector(*boundaryFacesLocalIndexesIt, newmanValues, newmanVector);
        addNewmanFaceLinearAdjusments<Basis::N_FUNCTIONS>(newmanVector, *determinantIt, linearTetrahedronsAdjusments[*boundaryFaceTetrahedronIndexIt]);

        ++initPointIt;
        ++transpJacobianMatrixIt;
        ++determinantIt;

        ++boundaryFacesLocalIndexesIt;
        ++boundaryFaceTetrahedronIndexIt;
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);

    delete[] newmanValues;
    delete[] integrationPoints;
}


void DGStefanTask::processVolumes(const double dt)
{
    const Volume* volumeIt = volumes;
    ElementsAdjusmentsSet* elementsAdjusmentsIt = volumesElementsAdjusments;
    CrossElementsAdjusmentsSet* elementsCrossAdjuesmentsIt = volumesInteriorFacesAdjusments;
    double** firstApproximationIt = volumesInitialX;
    InterfaceSideElementsSet* interfaceSideElementsIt = interfacesSidesElements;

    for (size_t i = 0; i < nVolumes; ++i)
    {
        size_t* tetrahedronsTags;
        size_t nTetrahedrons;

        size_t(*tetrahedronsNodesTags)[constants::tetrahedron::N_NODES];
        GMSHProxy::model::mesh::getTetrahedrons(tetrahedronsTags,
            nTetrahedrons,
            tetrahedronsNodesTags,
            volumeIt->tag);

        elementsAdjusmentsIt->tetrahedronsTags = tetrahedronsTags;
        elementsAdjusmentsIt->nTetrahedrons = nTetrahedrons;

        double (*bilinearElementsAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[nTetrahedrons][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        double (*linearElementsAdjusments)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];
        double (*volumeInitialX)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];

        elementsAdjusmentsIt->bilinear = (double*)bilinearElementsAdjusments;
        elementsAdjusmentsIt->linear = (double*)linearElementsAdjusments;
        *firstApproximationIt = (double*)volumeInitialX;

        double (*transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
        double* determinants;
        Coordinates* elementsInitPoints;

        GMSHProxy::model::mesh::getTetrahedronsJacobians(transpJacobianMatrixes, determinants, elementsInitPoints, nTetrahedrons, volumeIt->tag);

        double (*localJacobianMatrixes)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nTetrahedrons][LocalCoordinates3D::COUNT * Coordinates::COUNT];
        CoordinatesFunctions::computeLocalJacobianMatrixies(transpJacobianMatrixes, nTetrahedrons, determinants, localJacobianMatrixes);

        computeTetrahedronsAdjusments(localJacobianMatrixes,
            transpJacobianMatrixes,
            determinants,
            elementsInitPoints,
            nTetrahedrons,
            volumeIt->materialPhasePtr->thermalConductivity,
            volumeIt->materialPhasePtr->heatCapacity * volumeIt->materialPhasePtr->density,
            dt,
            volumeIt->materialPhasePtr->state == MaterialPhase::LIQUID ? initialConditionLiquid : initialConditionSolid,
            volumeInitialX,
            bilinearElementsAdjusments,
            linearElementsAdjusments);


        int interiorFacesEntityTag;
        int* boundariesFacesEntitiesTags = intBuffer;

        size_t(*interiorFacesTetrahedronsIndexes)[2];
        uint8_t(*interiorFacesLocalIndexes)[2];
        size_t interiorFaces;

        size_t** boundariesFacesTetrahedronsIndexes = (size_t**)ptrBuffer;
        uint8_t** boundariesFacesLocalIndexes = (uint8_t**)(ptrBuffer + volumeIt->nBoundaries);
        size_t* nBoundariesFaces = sizesBuffer;

        DTGeometryKernel::createVolumeFacesEntities(volumeIt->tag,
            volumeIt->boundariesTags,
            volumeIt->nBoundaries,
            volumeIt->boundariesConditions,
            interiorFacesEntityTag,
            interiorFaces,
            boundariesFacesEntitiesTags,
            nBoundariesFaces,
            interiorFacesTetrahedronsIndexes,
            interiorFacesLocalIndexes,
            boundariesFacesTetrahedronsIndexes,
            boundariesFacesLocalIndexes);

        elementsCrossAdjuesmentsIt->elementIndexes = interiorFacesTetrahedronsIndexes;
        elementsCrossAdjuesmentsIt->nCrossElements = interiorFaces;

        double (*interiorFacesAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[interiorFaces][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        elementsCrossAdjuesmentsIt->bilinear = (double*)interiorFacesAdjusments;

        computeInteriorFacesAdjusments(interiorFacesEntityTag,
            interiorFacesTetrahedronsIndexes,
            interiorFacesLocalIndexes,
            volumeIt->materialPhasePtr->thermalConductivity,
            localJacobianMatrixes,
            bilinearElementsAdjusments,
            interiorFacesAdjusments);

        delete[] interiorFacesLocalIndexes;

        const int* boundariesTagIt = volumeIt->boundariesTags;
        const int* boundaryFacesEntityTagIt = boundariesFacesEntitiesTags;

        const size_t** boundaryFacesTetrahedronsIndexesIt = boundariesFacesTetrahedronsIndexes;
        const uint8_t** boundaryFacesLocalIndexesIt = boundariesFacesLocalIndexes;

        const BoundaryCondition* boundaryConditionIt = volumeIt->boundariesConditions;

        for (size_t j = 0; j < volumeIt->nBoundaries; ++j)
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
                        volumeIt->materialPhasePtr->thermalConductivity,
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
                        volumeIt->materialPhasePtr->thermalConductivity,
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

            case BoundaryCondition::NONCONFORM_INTERFACE:
            {
                double (*boundaryTetrahedronsLocalJacobians)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nBoundariesFaces[j]][LocalCoordinates3D::COUNT * Coordinates::COUNT];
                Coordinates* boundaryTetrahedronsInitPoints = new Coordinates[nBoundariesFaces[j]];
                DTGeometryKernel::extractBoundaryTetrahedrons(localJacobianMatrixes, elementsInitPoints, *boundaryFacesTetrahedronsIndexesIt, nBoundariesFaces[j], boundaryTetrahedronsLocalJacobians, boundaryTetrahedronsInitPoints);

                interfaceSideElementsIt->indexes = *boundaryFacesTetrahedronsIndexesIt;
                interfaceSideElementsIt->localJacobians = boundaryTetrahedronsLocalJacobians;
                interfaceSideElementsIt->initPoints = boundaryTetrahedronsInitPoints;
                interfaceSideElementsIt->facesLocalIndexes = *boundaryFacesLocalIndexesIt;

                *boundaryFacesTetrahedronsIndexesIt = nullptr;
                *boundaryFacesLocalIndexesIt = nullptr;

                ++interfaceSideElementsIt;

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

        DTGeometryKernel::removeCreatedVolumeFacesEntities(boundariesFacesEntitiesTags, volumeIt->nBoundaries);

        ++volumeIt;
        ++elementsAdjusmentsIt;
        ++elementsCrossAdjuesmentsIt;
        ++firstApproximationIt;
    }
}

void computeTetrahedronsAdjusments(const double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                   const double(*transpMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                   const double* determinantIt,
                                   const Coordinates* initPointIt,
                                   const size_t nTetrahedrons,
                                   const MaterialPhase& materialPhase,
                                   const double dt,
                                   void* memoryBuffer,
                                   double(*initialXIt)[Basis::N_FUNCTIONS],
                                   double(*bilinearTetrahedronsAdjusmentsIt)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                   double(*linearTetrahedronsAdjusmentsIt)[Basis::N_FUNCTIONS])
{
    const double* templateMassMatrix = ElementDGCalculator<Basis>::getMassMatrix();
    double* stiffnessMatrix = (double*)memoryBuffer;
    double* localVector = stiffnessMatrix + ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double* initlialConditionValues = localVector + Basis::N_FUNCTIONS;

    const LocalCoordinates3D* localGradients = ElementDGCalculator<Basis>::getIntegrationLocalGradients();
    Coordinates* gradients = (Coordinates*)(initlialConditionValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);
    Coordinates* integrationPoints = gradients + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

    double lambda = materialPhase.thermalConductivity;
    double massCoeff = materialPhase.heatCapacity * materialPhase.density / dt;

    const LLt* massSLAESolverPtr = ElementDGCalculator<Basis>::getMassSLAESolverPtr();
    if (materialPhase.state == MaterialPhase::LIQUID)
    {
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < nTetrahedrons; tetrahedronIndex++)
        {
            CoordinatesFunctions::translate(*localJacobianMatrixIt, localGradients, NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
            ElementDGCalculator<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

            CoordinatesFunctions::translate(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints,
                NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps,
                *initPointIt,
                *transpMatrixIt,
                integrationPoints);

            DGStefanTask::computeLiquidTemperature(integrationPoints, NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, initlialConditionValues);
            ElementDGCalculator<Basis>::computePowerVector(initlialConditionValues, localVector);
            massSLAESolverPtr->solve(localVector, *initialXIt);

            computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                massCoeff,
                stiffnessMatrix,
                lambda,
                localVector,
                *determinantIt,
                *bilinearTetrahedronsAdjusmentsIt,
                *linearTetrahedronsAdjusmentsIt);

            ++initialXIt;
            ++initPointIt;
            ++determinantIt;
            ++transpMatrixIt;
            ++localJacobianMatrixIt;

            ++bilinearTetrahedronsAdjusmentsIt;
            ++linearTetrahedronsAdjusmentsIt;
        }
    }
    else
    {
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < nTetrahedrons; tetrahedronIndex++)
        {
            CoordinatesFunctions::translate(*localJacobianMatrixIt, localGradients, NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
            ElementDGCalculator<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

            CoordinatesFunctions::translate(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints,
                NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps,
                *initPointIt,
                *transpMatrixIt,
                integrationPoints);

            DGStefanTask::computeSolidTemperature(integrationPoints, NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, initlialConditionValues);
            ElementDGCalculator<Basis>::computePowerVector(initlialConditionValues, localVector);
            massSLAESolverPtr->solve(localVector, *initialXIt);

            computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                massCoeff,
                stiffnessMatrix,
                lambda,
                localVector,
                *determinantIt,
                *bilinearTetrahedronsAdjusmentsIt,
                *linearTetrahedronsAdjusmentsIt);

            ++initialXIt;
            ++initPointIt;
            ++determinantIt;
            ++transpMatrixIt;
            ++localJacobianMatrixIt;

            ++bilinearTetrahedronsAdjusmentsIt;
            ++linearTetrahedronsAdjusmentsIt;
        }
    }

}


void processVolumes(const Volume* volumeIt,
                    const unsigned int nVolumes,
                    const double penalty,
                    const double dt,
                    void* memoryBuffer,
                    ElementsAdjusmentsSet* elementsAdjusmentsIt,
                    CrossElementsAdjusmentsSet* elementsCrossAdjuesmentsIt,
                    double** firstApproximationIt,
                    NonconformInterface* nonconformInterfaces,
                    InterfaceSideElementsSet(*nonconformalInterfaceSideElementsSetIt)[2])
{
    for (size_t i = 0; i < nVolumes; ++i)
    {
        size_t(*tetrahedronsNodesTags)[constants::tetrahedron::N_NODES];
        GMSHProxy::model::mesh::getTetrahedrons(elementsAdjusmentsIt->tetrahedronsTags,
            elementsAdjusmentsIt->nTetrahedrons,
            tetrahedronsNodesTags,
            volumeIt->tag);

        const size_t nTetrahedrons = elementsAdjusmentsIt->nTetrahedrons;

        double (*bilinearElementsAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[nTetrahedrons][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        double (*linearElementsAdjusments)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];
        double (*volumeInitialX)[Basis::N_FUNCTIONS] = new double[nTetrahedrons][Basis::N_FUNCTIONS];

        elementsAdjusmentsIt->bilinear = (double*)bilinearElementsAdjusments;
        elementsAdjusmentsIt->linear = (double*)linearElementsAdjusments;
        *firstApproximationIt = (double*)volumeInitialX;

        double (*transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
        double* determinants;
        Coordinates* elementsInitPoints;

        GMSHProxy::model::mesh::getTetrahedronsJacobians(transpJacobianMatrixes, determinants, elementsInitPoints, volumeIt->tag);

        double (*localJacobianMatrixes)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nTetrahedrons][LocalCoordinates3D::COUNT * Coordinates::COUNT];
        CoordinatesFunctions::computeLocalJacobianMatrixies(transpJacobianMatrixes, nTetrahedrons, determinants, localJacobianMatrixes);

        computeTetrahedronsAdjusments(localJacobianMatrixes,
                                      transpJacobianMatrixes,
                                      determinants,
                                      elementsInitPoints,
                                      nTetrahedrons,
                                      *volumeIt->materialPhasePtr,
                                      dt,
                                      memoryBuffer,
                                      volumeInitialX,
                                      bilinearElementsAdjusments,
                                      linearElementsAdjusments);

        FacesSet interiorFacesSet;
        FacesSet* boundariesFacesSets = nullptr; //!

        DTGeometryKernel::getVolumeFacesSets(volumeIt->tag, volumeIt->boundaries, volumeIt->nBoundaries, interiorFacesSet, boundariesFacesSets, memoryBuffer);

        elementsCrossAdjuesmentsIt->elementIndexes = interiorFacesSet.elementIndexes;
        elementsCrossAdjuesmentsIt->nCrossElements = interiorFacesSet.count;

        double (*interiorFacesAdjusments)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = new double[interiorFacesSet.count][ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS];
        elementsCrossAdjuesmentsIt->bilinear = (double*)interiorFacesAdjusments;

        computeInteriorFacesAdjusments(interiorFacesSet,
                                       volumeIt->materialPhasePtr->thermalConductivity,
                                       penalty,
                                       localJacobianMatrixes,
                                       memoryBuffer,
                                       bilinearElementsAdjusments,
                                       interiorFacesAdjusments);


        delete[] interiorFacesLocalIndexes;

        const int* boundariesTagIt = volumeIt->boundariesTags;
        const int* boundaryFacesEntityTagIt = boundariesFacesEntitiesTags;

        const size_t** boundaryFacesTetrahedronsIndexesIt = boundariesFacesTetrahedronsIndexes;
        const uint8_t** boundaryFacesLocalIndexesIt = boundariesFacesLocalIndexes;

        const BoundaryCondition* boundaryConditionIt = volumeIt->boundariesConditions;

        for (size_t j = 0; j < volumeIt->nBoundaries; ++j)
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
                        volumeIt->materialPhasePtr->thermalConductivity,
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
                        volumeIt->materialPhasePtr->thermalConductivity,
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

            case BoundaryCondition::NONCONFORM_INTERFACE:
            {
                double (*boundaryTetrahedronsLocalJacobians)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = new double[nBoundariesFaces[j]][LocalCoordinates3D::COUNT * Coordinates::COUNT];
                Coordinates* boundaryTetrahedronsInitPoints = new Coordinates[nBoundariesFaces[j]];
                DTGeometryKernel::extractBoundaryTetrahedrons(localJacobianMatrixes, elementsInitPoints, *boundaryFacesTetrahedronsIndexesIt, nBoundariesFaces[j], boundaryTetrahedronsLocalJacobians, boundaryTetrahedronsInitPoints);

                interfaceSideElementsIt->indexes = *boundaryFacesTetrahedronsIndexesIt;
                interfaceSideElementsIt->localJacobians = boundaryTetrahedronsLocalJacobians;
                interfaceSideElementsIt->initPoints = boundaryTetrahedronsInitPoints;
                interfaceSideElementsIt->facesLocalIndexes = *boundaryFacesLocalIndexesIt;

                *boundaryFacesTetrahedronsIndexesIt = nullptr;
                *boundaryFacesLocalIndexesIt = nullptr;

                ++interfaceSideElementsIt;

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

        DTGeometryKernel::removeCreatedVolumeFacesEntities(boundariesFacesEntitiesTags, volumeIt->nBoundaries);

        ++volumeIt;
        ++elementsAdjusmentsIt;
        ++elementsCrossAdjuesmentsIt;
        ++firstApproximationIt;
    }
}