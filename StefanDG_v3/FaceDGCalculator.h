#pragma once
#include "LinearLagrangeBasis.h"
#include "NumericalIntegration.h"
#include "GeometryConstants.h"


class FaceDGCalculator
{
	using Basis = LinearLagrangeBasis;
	using NumericalIntegrationMethod = NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>;
	static const uint16_t N_BASIS_VALUES = Basis::N_FUNCTIONS * NumericalIntegrationMethod::nSteps;
	static const uint8_t N_LOCAL_MATRIX_ELEMENTS = Basis::N_FUNCTIONS * Basis::N_FUNCTIONS;

public:
	static void computeFlowMatrix(const uint8_t faceIndex, const double normalDerivatives[N_BASIS_VALUES], double flowMatrix[N_LOCAL_MATRIX_ELEMENTS]);
	static const double (*getValuesByFace())[N_BASIS_VALUES];
	static const LocalCoordinates3D (*getLocalgradientsByFace())[N_BASIS_VALUES];
	static const double* const* getMassMatrixies();
	static const double (*getCrossMassMatrixies())[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];

	static void init();
	static void finalize();
private:
	FaceDGCalculator();
	static void computeBasisDataByFace();
	static void computeMassMatrix(const double values[N_BASIS_VALUES], double massMatrix[N_LOCAL_MATRIX_ELEMENTS]);
	static void computeCrossMassMatrix(const double values1[N_BASIS_VALUES],
								const double values2[N_BASIS_VALUES],
								double crossMassMatrix1[N_LOCAL_MATRIX_ELEMENTS],
								double crossMassMatrix2[N_LOCAL_MATRIX_ELEMENTS]);
	static void computeMassMatrixies();

	static double integrateFlowFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionNormalDerivativeIt);
	static double integrateMassFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionValueIt);

	static double (*_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];
	static double** _massMatrixByFace;
	
	static double (*_valuesByFace)[N_BASIS_VALUES];
	static LocalCoordinates3D (*_localGradientsByFace)[N_BASIS_VALUES];
	

};


