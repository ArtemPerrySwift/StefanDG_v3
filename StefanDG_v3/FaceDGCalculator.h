#pragma once
#include "LinearLagrangeBasis.h"
#include "NumericalIntegration.h"
#include "GeometryConstants.h"

template <class Basis>
class FaceDGCalculator
{
	using NumericalIntegrationMethod = NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>;

public:
	static const uint16_t N_BASIS_VALUES = Basis::N_FUNCTIONS * NumericalIntegrationMethod::nSteps;
	static const uint8_t N_LOCAL_MATRIX_ELEMENTS = Basis::N_FUNCTIONS * Basis::N_FUNCTIONS;

	static void computeFlowMatrix(const uint8_t faceIndex, const double normalDerivatives[N_BASIS_VALUES], double flowMatrix[N_LOCAL_MATRIX_ELEMENTS]);
	static void computeFlowVector(const double normalDerivatives[N_BASIS_VALUES],
								  const double targetFunctionValues[NumericalIntegrationMethod::nSteps],
								  double flowVector[Basis::N_FUNCTIONS]);

	static void computeFlowVector(const double normalDerivatives[N_BASIS_VALUES],
								  double flowVector[Basis::N_FUNCTIONS]);

	static void computePowerVector(const uint8_t faceIndex, const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS]);

	static const double (*getValuesByFace())[N_BASIS_VALUES];
	static const LocalCoordinates3D (*getLocalGradientsByFace())[N_BASIS_VALUES];
	static const double* const* getMassMatrixies();
	static const double (*getCrossMassMatrixies())[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];
	static const double (*getMassVectors())[Basis::N_FUNCTIONS];

	static void init();
	static void finalize();
	static bool isInitialized();
private:
	FaceDGCalculator();

	static bool _initializationState;

	static void computeBasisDataByFace();

	static void computeMassMatrix(const double values[N_BASIS_VALUES], double massMatrix[N_LOCAL_MATRIX_ELEMENTS]);
	static void computeCrossMassMatrix(const double values1[N_BASIS_VALUES],
									   const double values2[N_BASIS_VALUES],
									   double crossMassMatrix1[N_LOCAL_MATRIX_ELEMENTS],
									   double crossMassMatrix2[N_LOCAL_MATRIX_ELEMENTS]);
	static void computeMassMatrixies();

	static void computeMassVector(const double* faceValueIt, double* massVectorIElementt);
	static void computeMassVectors();

	static double integrateFlowFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionNormalDerivativeIt);
	static double integrateFlowFunction(const double* iBasisFunctionNormalDerivativeIt);
	static double integrateMassFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionValueIt);
	static double integrateMassFunction(const double* iBasisFunctionValueIt);
	static double integratePowerFunction(const double* iBasisFunctionValueIt, const double* targetFunctionValueIt);

	static double (*_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];
	static double** _massMatrixByFace;
	static double (*_massVectorByFace)[Basis::N_FUNCTIONS];
	
	static double (*_valuesByFace)[N_BASIS_VALUES];
	static LocalCoordinates3D (*_localGradientsByFace)[N_BASIS_VALUES];
	

};


