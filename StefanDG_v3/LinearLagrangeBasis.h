#pragma once
#include <cstdint>
#include "LocalCoordinates3D.h"

class LinearLagrangeBasis
{
public:
	static const uint8_t N_FUNCTIONS = 4;
	static const uint8_t ORDER = 1;

	static void compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);
	static void compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	//static void compute(const LocalCoordinates* localPointIt, const size_t nPoints, double* gradientIt);

private:
	static double* compute0Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);
	static double* compute1Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);
	static double* compute2Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);
	static double* compute3Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);

	static LocalCoordinates3D* compute0Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	static LocalCoordinates3D* compute1Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	static LocalCoordinates3D* compute2Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	static LocalCoordinates3D* compute3Function(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
};

//const uint8_t LinearLagrangeBasis::N_FUNCTIONS = 4;
//const uint8_t LinearLagrangeBasis::ORDER = 1;

