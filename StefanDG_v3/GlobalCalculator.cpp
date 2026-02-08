#include "GlobalCalculator.h"
#include "Bases.h"

template<class Basis>
void GlobalCalculator<Basis>::computeMassMatrix()
{
	double* di = _massMatrix;
	const double* iFunctionValues = _values;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		*di = integrateMassFunction(iFunctionValues, iFunctionValues);
		const double* jFunctionValues = iFunctionValues + NumericalIntegrationMethod::nSteps;
		double* ijElement = di + 1;
		double* jiElement = di + NumericalIntegrationMethod::nSteps;
		for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
		{
			*ijElement = *jiElement = integrateMassFunction(iFunctionValues, jFunctionValues);
			++ijElement;

			jiElement += Basis::N_FUNCTIONS;
			jFunctionValues += NumericalIntegrationMethod::nSteps;
		}

		iFunctionValues += NumericalIntegrationMethod::nSteps;
		di += Basis::N_FUNCTIONS + 1;
	}
}

template<class Basis>
void GlobalCalculator<Basis>::computeMassVector()
{
	double* massVectorElementIt = _massVector;
	const double* iBasisFunctionValues = _values;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		*massVectorElementIt = integrateMassFunction(iBasisFunctionValues);

		iBasisFunctionValues += NumericalIntegrationMethod::nSteps;
		++massVectorElementIt;
	}
}

template<class Basis>
double GlobalCalculator<Basis>::integrateMassFunction(const double* iBasisFunctionValueIt)
{
	const double* weightsIt = NumericalIntegrationMethod::weights;
	double sum = 0.0;
	for (uint8_t k = 0; k < NumericalIntegrationMethod::nSteps; ++k)
	{
		sum += (*weightsIt) * (*iBasisFunctionValueIt);
		++weightsIt;
		++iBasisFunctionValueIt;
	}

	return sum;
}

template<class Basis>
double GlobalCalculator<Basis>::integrateMassFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionValueIt)
{
	const double* weightsIt = NumericalIntegrationMethod::weights;
	double sum = 0.0;
	for (uint8_t k = 0; k < NumericalIntegrationMethod::nSteps; ++k)
	{
		sum += (*weightsIt) * (*iBasisFunctionValueIt) * (*jBasisFunctionValueIt);
		++weightsIt;
		++iBasisFunctionValueIt;
		++jBasisFunctionValueIt;
	}
	
	return sum;
}

template<class Basis>
double GlobalCalculator<Basis>::integrateStiffnessFunction(const Coordinates* iBasisFunctionGradientIt, const Coordinates* jBasisFunctionGradientIt)
{
	const double* weightsIt = NumericalIntegrationMethod::weights;
	double sum = 0.0;
	for (uint8_t k = 0; k < NumericalIntegrationMethod::nSteps; ++k)
	{
		sum += (*weightsIt) * (iBasisFunctionGradientIt->x * jBasisFunctionGradientIt->x + iBasisFunctionGradientIt->y * jBasisFunctionGradientIt->y + iBasisFunctionGradientIt->z * jBasisFunctionGradientIt->z);
		++weightsIt;
		++iBasisFunctionGradientIt;
		++jBasisFunctionGradientIt;
	}

	return sum;
}

template<class Basis>
double GlobalCalculator<Basis>::integratePowerFunction(const double* iBasisFunctionValueIt, const double* targetFunctionValueIt)
{
	const double* weightsIt = NumericalIntegrationMethod::weights;
	double sum = 0.0;
	for (uint8_t k = 0; k < NumericalIntegrationMethod::nSteps; ++k)
	{
		sum += (*weightsIt) * (*iBasisFunctionValueIt) * (*targetFunctionValueIt);
		++weightsIt;
		++iBasisFunctionValueIt;
		++targetFunctionValueIt;
	}

	return sum;
}

template<class Basis>
void GlobalCalculator<Basis>::init()
{
	_memoryBuffer = new double[N_LOCAL_MATRIX_ELEMENTS + Basis::N_FUNCTIONS + N_BASIS_VALUES * (1 + LocalCoordinates3D::COUNT)];
	_massMatrix = _memoryBuffer;
	_massVector = _massMatrix + N_LOCAL_MATRIX_ELEMENTS;
	_values = _massVector + Basis::N_FUNCTIONS;
	_localGradients = (LocalCoordinates3D*)(_values + N_BASIS_VALUES);

	Basis::compute(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, _values);
	Basis::compute(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, _localGradients);

	computeMassMatrix();
	computeMassVector();
}

template<class Basis>
void GlobalCalculator<Basis>::finalize()
{
	delete[] _memoryBuffer;
	_memoryBuffer = nullptr;
	_massMatrix = nullptr;
	_massVector = nullptr;

	_values = nullptr;
	_localGradients = nullptr;
}

template<class Basis>
const double* GlobalCalculator<Basis>::getMassMatrix()
{
	return _massMatrix;
}

template<class Basis>
const double* GlobalCalculator<Basis>::getMassVector()
{
	return _massVector;
}

template<class Basis>
const double* GlobalCalculator<Basis>::getIntegrationValues()
{
	return _values;
}

template<class Basis>
const LocalCoordinates3D* GlobalCalculator<Basis>::getIntegrationLocalGradients()
{
	return _localGradients;
}

template<class Basis>
void GlobalCalculator<Basis>::computeStiffnessMatrix(const Coordinates gradients[N_BASIS_VALUES], double stiffnessMatrix[N_LOCAL_MATRIX_ELEMENTS])
{
	double* di = stiffnessMatrix;
	const Coordinates* iFunctionGradients = gradients;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		*di = integrateStiffnessFunction(iFunctionGradients, iFunctionGradients);
		const Coordinates* jFunctionGradients = iFunctionGradients + NumericalIntegrationMethod::nSteps;
		double* ijElement = di + 1;
		double* jiElement = di + NumericalIntegrationMethod::nSteps;
		for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
		{
			*ijElement = *jiElement = integrateStiffnessFunction(iFunctionGradients, jFunctionGradients);
			++ijElement;

			++ijElement;
			jiElement += Basis::N_FUNCTIONS;
			jFunctionGradients += NumericalIntegrationMethod::nSteps;
		}

		iFunctionGradients += NumericalIntegrationMethod::nSteps;
		di += Basis::N_FUNCTIONS + 1;
	}
}

template<class Basis>
void GlobalCalculator<Basis>::computePowerVector(const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS])
{
	const double* iBasisFunctionValues = _values;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		*powerVector = integratePowerFunction(iBasisFunctionValues, targetFunctionValues);

		iBasisFunctionValues += NumericalIntegrationMethod::nSteps;
		++powerVector;
	}
}

template<class Basis>
double* GlobalCalculator<Basis>::_memoryBuffer = nullptr;

template<class Basis>
double* GlobalCalculator<Basis>::_massMatrix = nullptr;

template<class Basis>
double* GlobalCalculator<Basis>::_massVector = nullptr;

template<class Basis>
double* GlobalCalculator<Basis>::_values = nullptr;

template<class Basis>
LocalCoordinates3D* GlobalCalculator<Basis>::_localGradients = nullptr;

#define X(BasisName) template class GlobalCalculator<BasisName>;
BASES
#undef X;
