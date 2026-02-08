#pragma once
#include "LocalCoordinates3D.h"
#include "LocalCoordinates2D.h"
#include <cstdint>
#include <cmath>

namespace NumericalIntegration
{
	namespace Tetrahedron
	{
		template <size_t ORDER> requires (ORDER == 2)
		class Gauss
		{
		public:
			Gauss() = default;
			static const LocalCoordinates3D localPoints[];
			static const double weights[];
			static const uint8_t nSteps;
		};

		template<>
		class Gauss<2>
		{
		public:
			static const uint8_t nSteps = 4;
			static const LocalCoordinates3D localPoints[nSteps];
			static const double weights[nSteps];
		private:
			Gauss<2>() = default;
			static const double coord1;
			static const double coord2;
			static const double weight;
		};

	}

	namespace Triangle
	{
		template <size_t ORDER> requires (ORDER == 2)
			class Gauss
		{
		public:
			Gauss() = default;
			static const LocalCoordinates2D localPoints[];
			static const double weights[];
			static const uint8_t nSteps;
		};

		template<>
		class Gauss<2>
		{
		public:
			static const uint8_t nSteps = 3;
			static const LocalCoordinates2D localPoints[nSteps];
			static const double weights[nSteps];
		private:
			Gauss<2>() = default;
			static const double coord1;
			static const double coord2;
			static const double weight;
		};

	}


}
