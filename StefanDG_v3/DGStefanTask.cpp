#include "DGStefanTask.h"
#include "CoordinatesFunctions.h"
#include "GMSHProxy.h"
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

void DGStefanTask::computeTetrahedronsAdjusments(const double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT], const double(*transpMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT], const double* determinantIt, const Coordinates* initPointIt, const size_t nTetrahedrons, const double lambda, const double gamma, const double dt, double(*initialCondition)(const Coordinates& point), double(*initialXIt)[Basis::N_FUNCTIONS], double(*bilinearTetrahedronsAdjusmentsIt)[ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS], double(*linearTetrahedronsAdjusmentsIt)[Basis::N_FUNCTIONS])
{
    const double* templateMassMatrix = ElementDGCalculator<Basis>::getMassMatrix();
    double* stiffnessMatrix = dblBuffer;
    double* localVector = stiffnessMatrix + ElementDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double* initlialConditionValues = localVector + Basis::N_FUNCTIONS;

    const LocalCoordinates3D* localGradients = ElementDGCalculator<Basis>::getIntegrationLocalGradients();
    Coordinates* gradients = coordinatesBuffer;
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


void DGStefanTask::computeInteriorFacesAdjusments(const int interiorFacesEntityTag,
                                                  const size_t(*interiorFaceTetrahedronsIndexesIt)[2],
                                                  const uint8_t(*interiorFaceLocalIndexesIt)[2],
                                                  const double lambda,
                                                  const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                  double(bilinearTetrahedronsAdjusments)[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                  double(*bilinearInterfaceAdjusmentsIt)[FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS])
{
    const double* const* massMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceDGCalculator<Basis>::getCrossMassMatrixies();

    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();

    double* flowMatrix = dblBuffer;
    double* crossFlowMatrix1 = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* crossFlowMatrix2 = crossFlowMatrix1 + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    Coordinates* gradients = coordinatesBuffer;
    double* normalDerivatives = crossFlowMatrix2 + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    size_t nFaces;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, nFaces, interiorFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    size_t(*facesNodesTags)[constants::triangle::N_NODES];
    Coordinates(*facesNodes)[constants::triangle::N_NODES];
    GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, nFaces, interiorFacesEntityTag);

    Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
    GMSHProxy::free(facesNodesTags);


    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double facePenalty = _penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;
        double commonMultiplier = *determinantIt * lambda;

        Coordinates normal;
        CoordinatesFunctions::computeNormal(*faceNodesIt, normal);

        const size_t* faceTetrahedronIndexeIt = *interiorFaceTetrahedronsIndexesIt;
        const uint8_t* faceLocalIndexIt = *interiorFaceLocalIndexesIt;

        CoordinatesFunctions::translate(localJacobianMatrix[faceTetrahedronIndexeIt[0]],
            gradientsByFace[faceLocalIndexIt[0]],
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[0], normalDerivatives, flowMatrix);
        addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[faceLocalIndexIt[0]], _penalty, commonMultiplier, bilinearTetrahedronsAdjusments[*faceTetrahedronIndexeIt]);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[1], normalDerivatives, crossFlowMatrix1);

        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;

        CoordinatesFunctions::translate(localJacobianMatrix[faceTetrahedronIndexeIt[1]],
            gradientsByFace[faceLocalIndexIt[1]],
            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
            gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceDGCalculator<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[1], normalDerivatives, flowMatrix);
        addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[faceLocalIndexIt[0]], _penalty, commonMultiplier, bilinearTetrahedronsAdjusments[*faceTetrahedronIndexeIt]);

        FaceDGCalculator<Basis>::computeFlowMatrix(faceLocalIndexIt[0], normalDerivatives, crossFlowMatrix2);

        computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                                                              crossFlowMatrix1,
                                                              crossMassMatrixies[faceLocalIndexIt[0]][faceLocalIndexIt[1]],
                                                              facePenalty,
                                                              commonMultiplier,
                                                              *bilinearInterfaceAdjusmentsIt);

        ++interiorFaceTetrahedronsIndexesIt;
        ++interiorFaceLocalIndexesIt;
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


void DGStefanTask::computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                              const size_t* boundaryFaceTetrahedronIndexIt,
                                              const uint8_t* boundaryFacesLocalIndexesIt,
                                              const Coordinates outwardNormal,
                                              const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                              double(* const dirichletCondition)(const Coordinates&),
                                              const double lambda,
                                              double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                              double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    double* flowMatrix = dblBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* powerVector = flowVector + Basis::N_FUNCTIONS;

    double* dirichletValues = powerVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = coordinatesBuffer;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;


    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    size_t nFaces;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, nFaces, boundaryFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double determinant = *determinantIt;
        double facePenalty = _penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
        double commonMultiplier = determinant * lambda;

        size_t tetrahedronIndex = *boundaryFaceTetrahedronIndexIt;
        uint8_t faceLocalIndex = *boundaryFacesLocalIndexesIt;

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


void DGStefanTask::computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                                   const size_t* boundaryFaceTetrahedronIndexIt,
                                                   const uint8_t* boundaryFacesLocalIndexesIt,
                                                   const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                   double(* const dirichletCondition)(const Coordinates&),
                                                   const double lambda,
                                                   double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                   double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    double* flowMatrix = dblBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;
    double* powerVector = flowVector + Basis::N_FUNCTIONS;

    double* dirichletValues = powerVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = coordinatesBuffer;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    size_t nFaces;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, nFaces, boundaryFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    size_t(*facesNodesTags)[constants::triangle::N_NODES];
    Coordinates(*facesNodes)[constants::triangle::N_NODES];
    GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, nFaces, boundaryFacesEntityTag);

    Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
    GMSHProxy::free(facesNodesTags);

    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double determinant = *determinantIt;
        double facePenalty = _penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
        double commonMultiplier = determinant * lambda;

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

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);
    GMSHProxy::free(facesNodes);
}

void DGStefanTask::computeDirichletFacesAdjusments(const int boundaryFacesEntityTag,
                                                   const size_t* boundaryFaceTetrahedronIndexIt,
                                                   const uint8_t* boundaryFacesLocalIndexesIt,
                                                   const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                   const double condition,
                                                   const double lambda,
                                                   double bilinearTetrahedronsAdjusments[][FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                   double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
{
    const double* const* templatePenaltyMatrixByFace = FaceDGCalculator<Basis>::getMassMatrixies();
    const LocalCoordinates3D(*gradientsByFace)[FaceDGCalculator<Basis>::N_BASIS_VALUES] = FaceDGCalculator<Basis>::getLocalGradientsByFace();
    const double (*massVectorByFace)[Basis::N_FUNCTIONS] = FaceDGCalculator<Basis>::getMassVectors();
    double* flowMatrix = dblBuffer;

    double* flowVector = flowMatrix + FaceDGCalculator<Basis>::N_LOCAL_MATRIX_ELEMENTS;

    double* dirichletValues = flowVector + Basis::N_FUNCTIONS;
    double* normalDerivatives = dirichletValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    Coordinates* integrationPoints = coordinatesBuffer;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
    Coordinates* gradients = integrationPointsEndIt;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    size_t nFaces;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, nFaces, boundaryFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    size_t(*facesNodesTags)[constants::triangle::N_NODES];
    Coordinates(*facesNodes)[constants::triangle::N_NODES];
    GMSHProxy::model::mesh::getTrianglesNodes(facesNodesTags, facesNodes, nFaces, boundaryFacesEntityTag);

    Coordinates(*faceNodesIt)[constants::triangle::N_NODES] = facesNodes;
    GMSHProxy::free(facesNodesTags);

    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        double determinant = *determinantIt;
        double facePenalty = _penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;
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


void DGStefanTask::computeNewmanFacesAdjusments(const int boundaryFacesEntityTag, const size_t* boundaryFacesTetrahedronIndexesIt, const uint8_t* boundaryFacesLocalIndexesIt, double(* const newmanCondition)(const Coordinates&), double linearTetrahedronsAdjusments[][LinearLagrangeBasis::N_FUNCTIONS])
{
    Coordinates* integrationPoints = coordinatesBuffer;
    Coordinates* integrationPointsEndIt = integrationPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
    double* newmanValues = dblBuffer;
    double* newmanVector = newmanValues + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

    double (*transpJacobianMatrixies)[Coordinates::COUNT * LocalCoordinates3D::COUNT];
    double* determinants;
    Coordinates* initPoints;
    size_t nFaces;
    GMSHProxy::model::mesh::getTriangleJacobians(transpJacobianMatrixies, determinants, initPoints, nFaces, boundaryFacesEntityTag);

    const double (*transpJacobianMatrixIt)[Coordinates::COUNT * LocalCoordinates3D::COUNT] = transpJacobianMatrixies;
    const double* determinantIt = determinants;
    const Coordinates* initPointIt = initPoints;

    for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
    {
        size_t tetrahedronIndex = *boundaryFacesTetrahedronIndexesIt;
        uint8_t faceLocalIndex = *boundaryFacesLocalIndexesIt;

        CoordinatesFunctions::translate(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, *initPointIt, *transpJacobianMatrixIt, integrationPoints);
        std::transform(integrationPoints, integrationPointsEndIt, newmanValues, newmanCondition);

        FaceDGCalculator<Basis>::computePowerVector(*boundaryFacesLocalIndexesIt, newmanValues, newmanVector);
        addNewmanFaceLinearAdjusments<Basis::N_FUNCTIONS>(newmanVector, *determinantIt, linearTetrahedronsAdjusments[*boundaryFacesTetrahedronIndexesIt]);

        ++initPointIt;
        ++transpJacobianMatrixIt;
        ++determinantIt;

        ++boundaryFacesLocalIndexesIt;
        ++boundaryFacesTetrahedronIndexesIt;
    }

    GMSHProxy::free(transpJacobianMatrixies);
    GMSHProxy::free(determinants);
    GMSHProxy::free(initPoints);

    delete[] newmanValues;
    delete[] integrationPoints;
}
