#pragma once
#include "LocalCoordinates3D.h"
#include "Coordinates.h"
#include "NumericalIntegration.h"

template<class Basis>
class GlobalCalculator
{
	//using Basis = LinearLagrangeBasis;
	using NumericalIntegrationMethod = NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>;
public:
	static const uint16_t N_BASIS_VALUES = Basis::N_FUNCTIONS * NumericalIntegrationMethod::nSteps;
	static const uint8_t N_LOCAL_MATRIX_ELEMENTS = Basis::N_FUNCTIONS * Basis::N_FUNCTIONS;

	static void init();
	static void finalize();
	static const double* getMassMatrix();
	static const double* getMassVector();
	static const double* getIntegrationValues();
	static const LocalCoordinates3D* getIntegrationLocalGradients();
	
	static void computeStiffnessMatrix(const Coordinates gradients[N_BASIS_VALUES], double stiffnessMatrix[N_LOCAL_MATRIX_ELEMENTS]);
	static void computePowerVector(const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS]);

private:
	GlobalCalculator() = default;
	static void computeMassMatrix();
	static void computeMassVector();

	static double integrateMassFunction(const double* iBasisFunctionValueIt);
	static double integrateMassFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionValueIt);
	static double integrateStiffnessFunction(const Coordinates* iBasisFunctionGradientIt, const Coordinates* jBasisFunctionGradientIt);
	static double integratePowerFunction(const double* iBasisFunctionValueIt, const double* targetFunctionValueIt);

	static double* _memoryBuffer;
	static double* _massMatrix;
	static double* _massVector;

	static double* _values;
	static LocalCoordinates3D* _localGradients;
};


