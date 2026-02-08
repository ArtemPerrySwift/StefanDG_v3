#include "FEM.h"


void FEM::calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod, const double functionsSetValues[], const uint8_t nFunctions, double bilinearMatrix[])
{
    const uint8_t pointsCount = integrationMethod.getStepsCount();

    const double* iFunctionSetValues = functionsSetValues;
    double* iBilinearMatrixRow = bilinearMatrix;

    for (uint8_t i = 0; i < nFunctions; ++i)
    {
        iBilinearMatrixRow[i] = integrationMethod.integrateMultiplication(iFunctionSetValues, iFunctionSetValues);

        const double* jFunctionSetValues = iFunctionSetValues + pointsCount;
        double* jBilinearMatrixRow = bilinearMatrix + nFunctions;

        for (uint8_t j = i + 1; j < nFunctions; ++j)
        {
            jBilinearMatrixRow[i] = iBilinearMatrixRow[j] = integrationMethod.integrateMultiplication(iFunctionSetValues, jFunctionSetValues);
            jFunctionSetValues += pointsCount;
            jBilinearMatrixRow += nFunctions;
        }

        iFunctionSetValues += pointsCount;
        iBilinearMatrixRow += nFunctions;
    }
}

void FEM::calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod, const double firstFunctionsSetValues[], const double secondFunctionsSetValues[], const uint8_t nFunctions, double bilinearMatrix[])
{
    const uint8_t pointsCount = integrationMethod.getStepsCount();

    const double* iFunctionFirstSetValues = firstFunctionsSetValues;
    //const double* jFunctionSecondSetValuesBuffer = secondFunctionsSetValues;
    double* ijBilinearMatrixElementPtr = bilinearMatrix;

    for (uint8_t i = 0; i < nFunctions; ++i)
    {
        ijBilinearMatrixElementPtr += i;
        const double* jFunctionSecondSetValues = secondFunctionsSetValues;
        *ijBilinearMatrixElementPtr = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSecondSetValues);
        double* jiBilinearMatrixElementPtr = ijBilinearMatrixElementPtr;

        for (uint8_t j = i + 1; j < nFunctions; ++j)
        {
            ++ijBilinearMatrixElementPtr;
            jiBilinearMatrixElementPtr += nFunctions;
            jFunctionSecondSetValues += pointsCount;

            *jiBilinearMatrixElementPtr = *ijBilinearMatrixElementPtr = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSecondSetValues);
        }

        ++ijBilinearMatrixElementPtr;
        iFunctionFirstSetValues += pointsCount;
        secondFunctionsSetValues += pointsCount;
    }
}

void FEM::calcBilinearMatrix(const IIntegrationMethod& integrationMethod, const double firstFunctionsSetValues[], const double secondFunctionsSetValues[], const uint8_t nFunctions, double bilinearMatrix[])
{
    const uint8_t pointsCount = integrationMethod.getStepsCount();

    const double* iFunctionFirstSetValues = firstFunctionsSetValues;

    for (uint8_t i = 0, k = 0; i < nFunctions; ++i)
    {
        const double* jFunctionSetValues = secondFunctionsSetValues;

        for (uint8_t j = 0; j < nFunctions; ++j, ++k)
        {
            bilinearMatrix[k] = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSetValues);
            jFunctionSetValues += pointsCount;
        }

        iFunctionFirstSetValues += pointsCount;
    }
}

void FEM::calcBilinearMatrix(const IIntegrationMethod& integrationMethod, const double firstFunctionsSetValues[], const uint8_t nFirstSetFunctions, const double secondFunctionsSetValues[], const uint8_t nSecondSetFunctions, double bilinearMatrix[])
{
    const uint8_t pointsCount = integrationMethod.getStepsCount();

    const double* iFunctionFirstSetValues = firstFunctionsSetValues;

    for (uint8_t i = 0, k = 0; i < nFirstSetFunctions; ++i)
    {
        const double* jFunctionSetValues = secondFunctionsSetValues;

        for (uint8_t j = 0; j < nSecondSetFunctions; ++j, ++k)
        {
            bilinearMatrix[k] = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSetValues);
            jFunctionSetValues += pointsCount;
        }

        iFunctionFirstSetValues += pointsCount;
    }
}

void FEM::calcLinearVector(const IIntegrationMethod& integrationMethod, const double functionsSetValues[], const double targetFunctionValues[], const uint8_t nFunctions, double linearVector[])
{
    const uint8_t pointsCount = integrationMethod.getStepsCount();
    const double* iFunctionSetValues = functionsSetValues;
    for (uint8_t i = 0; i < nFunctions; i++)
    {
        linearVector[i] = integrationMethod.integrateMultiplication(iFunctionSetValues, targetFunctionValues);
        iFunctionSetValues += pointsCount;
    }
}
