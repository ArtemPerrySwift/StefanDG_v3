#include "FaceDGCalculator.h"
#include "CoordinatesFunctions.h"
void FaceDGCalculator::computeFlowMatrix(const uint8_t faceIndex, const double normalDerivatives[N_BASIS_VALUES], double flowMatrix[N_LOCAL_MATRIX_ELEMENTS])
{
	const double* iFunctionValues = _valuesByFace[faceIndex];
	double* ijElement = flowMatrix;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		const double* jFunctionNormalDerivatives = normalDerivatives;
		for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
		{
			*ijElement = integrateFlowFunction(iFunctionValues, jFunctionNormalDerivatives);
			++ijElement;

			jFunctionNormalDerivatives += NumericalIntegrationMethod::nSteps;
		}
		iFunctionValues += NumericalIntegrationMethod::nSteps;
	}
}

const double(*FaceDGCalculator::getValuesByFace())[N_BASIS_VALUES]
{
	return _valuesByFace;
}

const LocalCoordinates3D(*FaceDGCalculator::getLocalgradientsByFace())[N_BASIS_VALUES]
{
	return _localGradientsByFace;
}

const double* const* FaceDGCalculator::getMassMatrixies()
{
	return _massMatrixByFace;
}

const double(*FaceDGCalculator::getCrossMassMatrixies())[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS]
{
	return _massCrossMatrixByFaces;
}

void FaceDGCalculator::init()
{
	_massCrossMatrixByFaces = new double[constants::tetrahedron::N_FACES][constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];
	_massMatrixByFace = new double*[constants::tetrahedron::N_FACES];
	_valuesByFace = new double[constants::tetrahedron::N_FACES][N_BASIS_VALUES];
	_localGradientsByFace = new LocalCoordinates3D[constants::tetrahedron::N_FACES][N_BASIS_VALUES];

	computeBasisDataByFace();
	computeMassMatrixies();
}

void FaceDGCalculator::finalize()
{
	delete[] _massCrossMatrixByFaces;
	delete[] _massMatrixByFace;
	delete[] _valuesByFace;
	delete[] _localGradientsByFace;

	_massCrossMatrixByFaces = nullptr;
	_massMatrixByFace = nullptr;
	_valuesByFace = nullptr;
	_localGradientsByFace = nullptr;
}

double (*FaceDGCalculator::_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS] = nullptr;
double** FaceDGCalculator::_massMatrixByFace = nullptr;

double (*FaceDGCalculator::_valuesByFace)[N_BASIS_VALUES] = nullptr;
LocalCoordinates3D(*FaceDGCalculator::_localGradientsByFace)[N_BASIS_VALUES] = nullptr;

void FaceDGCalculator::computeBasisDataByFace()
{
	LocalCoordinates3D localCoordinates3D[NumericalIntegrationMethod::nSteps];

	CoordinatesFunctions::translate0FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[0]);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[0]);

	CoordinatesFunctions::translate1FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[1]);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[1]);

	CoordinatesFunctions::translate2FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[2]);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[2]);

	CoordinatesFunctions::translate3FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[3]);
	Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[3]);

}

void FaceDGCalculator::computeMassMatrix(const double values[N_BASIS_VALUES], double massMatrix[N_LOCAL_MATRIX_ELEMENTS])
{
	double* di = massMatrix;
	const double* iFunctionValues = values;
	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		*di = integrateMassFunction(iFunctionValues, iFunctionValues);
		const double* jFunctionValues = iFunctionValues;
		double* ijElement = di;
		double* jiElement = di;
		for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
		{
			++ijElement;
			jiElement += Basis::N_FUNCTIONS;
			jFunctionValues += NumericalIntegrationMethod::nSteps;

			*ijElement = *jiElement = integrateMassFunction(iFunctionValues, jFunctionValues);

		}

		iFunctionValues += NumericalIntegrationMethod::nSteps;
		di += Basis::N_FUNCTIONS + 1;
	}
}

void FaceDGCalculator::computeCrossMassMatrix(const double values1[N_BASIS_VALUES], const double values2[N_BASIS_VALUES], double crossMassMatrix1[N_LOCAL_MATRIX_ELEMENTS], double crossMassMatrix2[N_LOCAL_MATRIX_ELEMENTS])
{
	const double* iFunctionValues = values1;
	double* ijElementMatrix1 = crossMassMatrix1;
	double* matrix2InitRowElementIt = crossMassMatrix2;

	for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
	{
		const double* jFunctionValues = values2;
		double* jiElementMatrix2 = matrix2InitRowElementIt;
		for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
		{
			*ijElementMatrix1 = *jiElementMatrix2 = integrateFlowFunction(iFunctionValues, jFunctionValues);
			++ijElementMatrix1;

			jiElementMatrix2 += Basis::N_FUNCTIONS;
			jFunctionValues += NumericalIntegrationMethod::nSteps;
		}
		iFunctionValues += NumericalIntegrationMethod::nSteps;
		++matrix2InitRowElementIt;
	}
}

void FaceDGCalculator::computeMassMatrixies()
{
	double** massMetrixIt = _massMatrixByFace;
	double (*iiMassMatrixPtr)[N_LOCAL_MATRIX_ELEMENTS] = *_massCrossMatrixByFaces;
	const double (*iFaceValuesIt)[N_BASIS_VALUES] = _valuesByFace;

	for (uint8_t i = 0; i < constants::tetrahedron::N_FACES; ++i)
	{
		*massMetrixIt = *iiMassMatrixPtr;
		computeMassMatrix(*iFaceValuesIt, *iiMassMatrixPtr);
		double (*ijMassMatrixIt)[N_LOCAL_MATRIX_ELEMENTS] = iiMassMatrixPtr;
		double (*jiMassMatrixIt)[N_LOCAL_MATRIX_ELEMENTS] = iiMassMatrixPtr;

		const double (*jFaceValuesIt)[N_BASIS_VALUES] = iFaceValuesIt;

		for (uint8_t j = i + 1; j < constants::tetrahedron::N_FACES; ++j)
		{
			++jFaceValuesIt;
			++ijMassMatrixIt;
			jiMassMatrixIt += constants::tetrahedron::N_FACES;

			computeCrossMassMatrix(*iFaceValuesIt, *jFaceValuesIt, *ijMassMatrixIt, *jiMassMatrixIt);
		}

		iiMassMatrixPtr += constants::tetrahedron::N_FACES + 1;
		++massMetrixIt;
	}
}

double FaceDGCalculator::integrateFlowFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionNormalDerivativeIt)
{
	const double* weightsIt = NumericalIntegrationMethod::weights;
	double sum = 0.0;
	for (uint8_t k = 0; k < NumericalIntegrationMethod::nSteps; ++k)
	{
		sum += (*weightsIt) * (*iBasisFunctionValueIt) * (*jBasisFunctionNormalDerivativeIt);
		++weightsIt;
		++iBasisFunctionValueIt;
		++jBasisFunctionNormalDerivativeIt;
	}

	return sum;
}

double FaceDGCalculator::integrateMassFunction(const double* iBasisFunctionValueIt, const double* jBasisFunctionValueIt)
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
