#pragma once
#include "Integration\IIntegrationMethod.h"

namespace FEM
{
    void calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod, const double functionsSetValues[], const uint8_t nFunctions, double bilinearMatrix[]);

    void calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod,
                                     const double firstFunctionsSetValues[],
                                     const double secondFunctionsSetValues[],
                                     const uint8_t nFunctions,
                                     double bilinearMatrix[]);

    void calcBilinearMatrix(const IIntegrationMethod& integrationMethod,
                            const double firstFunctionsSetValues[],
                            const double secondFunctionsSetValues[],
                            const uint8_t nFunctions,
                            double bilinearMatrix[]);

    void calcBilinearMatrix(const IIntegrationMethod& integrationMethod,
                            const double firstFunctionsSetValues[],
                            const uint8_t  nFirstSetFunctions,
                            const double secondFunctionsSetValues[],
                            const uint8_t  nSecondSetFunctions,
                            double bilinearMatrix[]);

    void calcLinearVector(const IIntegrationMethod& integrationMethod, const double functionsSetValues[], const double targetFunctionValues[], const uint8_t nFunctions, double linearVector[]);

    template<const uint8_t N_FUNCTIONS>
    void calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod,
                                     const double functionsSetValues[],
                                     double bilinearMatrix[N_FUNCTIONS * N_FUNCTIONS])
    {
        const uint8_t nPoints = integrationMethod.getStepsCount();

        const double* iFunctionSetValues = functionsSetValues;
        double* iBilinearMatrixRow = bilinearMatrix;

        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            iBilinearMatrixRow[i] = integrationMethod.integrateMultiplication(iFunctionSetValues, iFunctionSetValues);

            const double* jFunctionSetValues = iFunctionSetValues + nPoints;
            double* jBilinearMatrixRow = bilinearMatrix + N_FUNCTIONS;

            for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
            {
                jBilinearMatrixRow[i] = iBilinearMatrixRow[j] = integrationMethod.integrateMultiplication(iFunctionSetValues, jFunctionSetValues);
                jFunctionSetValues += nPoints;
                jBilinearMatrixRow += N_FUNCTIONS;
            }

            iFunctionSetValues += nPoints;
            iBilinearMatrixRow += N_FUNCTIONS;
        }
    }

    template<const uint8_t N_FUNCTIONS>
    void calcBilinearSymmetricMatrix(const IIntegrationMethod& integrationMethod,
                                     const double firstFunctionsSetValues[],
                                     const double secondFunctionsSetValues[],
                                     double bilinearMatrix[N_FUNCTIONS * N_FUNCTIONS])
    {
        const uint8_t nPoints = integrationMethod.getStepsCount();

        const double* iFunctionFirstSetValues = firstFunctionsSetValues;
        //const double* jFunctionSecondSetValuesBuffer = secondFunctionsSetValues;
        double* ijBilinearMatrixElementPtr = bilinearMatrix;

        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            ijBilinearMatrixElementPtr += i;
            const double* jFunctionSecondSetValues = secondFunctionsSetValues;
            *ijBilinearMatrixElementPtr = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSecondSetValues);
            double* jiBilinearMatrixElementPtr = ijBilinearMatrixElementPtr;

            for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
            {
                ++ijBilinearMatrixElementPtr;
                jiBilinearMatrixElementPtr += N_FUNCTIONS;
                jFunctionSecondSetValues += nPoints;

                *jiBilinearMatrixElementPtr = *ijBilinearMatrixElementPtr = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSecondSetValues);
            }

            ++ijBilinearMatrixElementPtr;
            iFunctionFirstSetValues += nPoints;
            secondFunctionsSetValues += nPoints;
        }
    }

    template<const uint8_t N_FUNCTIONS>
    void calcBilinearMatrix(const IIntegrationMethod& integrationMethod,
                            const double firstFunctionsSetValues[],
                            const double secondFunctionsSetValues[],
                            double bilinearMatrix[N_FUNCTIONS * N_FUNCTIONS])
    {
        const uint8_t pointsCount = integrationMethod.getStepsCount();

        const double* iFunctionFirstSetValues = firstFunctionsSetValues;

        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            const double* jFunctionSetValues = secondFunctionsSetValues;

            for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
            {
                *bilinearMatrix = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSetValues);
                jFunctionSetValues += pointsCount;
                ++bilinearMatrix;
            }

            iFunctionFirstSetValues += pointsCount;
        }
    }

    template<const uint8_t N_FIRST_SET_FUNCTIONS, const uint8_t N_SECOND_SET_FUNCTIONS>
    void calcBilinearMatrix(const IIntegrationMethod& integrationMethod,
        const double firstFunctionsSetValues[],
        const double secondFunctionsSetValues[],
        double bilinearMatrix[N_FIRST_SET_FUNCTIONS * N_SECOND_SET_FUNCTIONS])
    {
        const uint8_t pointsCount = integrationMethod.getStepsCount();

        const double* iFunctionFirstSetValues = firstFunctionsSetValues;

        for (uint8_t i = 0; i < N_FIRST_SET_FUNCTIONS; ++i)
        {
            const double* jFunctionSetValues = secondFunctionsSetValues;

            for (uint8_t j = 0; j < N_SECOND_SET_FUNCTIONS; ++j)
            {
                *bilinearMatrix = integrationMethod.integrateMultiplication(iFunctionFirstSetValues, jFunctionSetValues);
                
                jFunctionSetValues += pointsCount;
                ++*bilinearMatrix;
            }

            iFunctionFirstSetValues += pointsCount;
        }
    }

    template<const uint8_t N_FUNCTIONS>
    void calcLinearVector(const IIntegrationMethod& integrationMethod, 
                          const double functionsSetValues[], 
                          const double targetFunctionValues[],
                          double linearVector[N_FUNCTIONS])
    {
        const uint8_t pointsCount = integrationMethod.getStepsCount();
        const double* iFunctionSetValues = functionsSetValues;
        for (uint8_t i = 0; i < N_FUNCTIONS; i++)
        {
            *linearVector = integrationMethod.integrateMultiplication(iFunctionSetValues, targetFunctionValues);
            iFunctionSetValues += pointsCount;
            ++linearVector;
        }
    }
};

