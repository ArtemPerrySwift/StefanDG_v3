#include "NumericalIntegration.h"

namespace NumericalIntegration
{
	namespace Tetrahedron
	{
		const double Gauss<2>::coord1 = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
		const double Gauss<2>::coord2 = (5.0 - sqrt(5.0)) / 20.0;

		const double Gauss<2>::weight = 1.0 / 24.0;

		const LocalCoordinates3D Gauss<2>::localPoints[] = { {coord1, coord2, coord2},
													 {coord2, coord1, coord2},
													 {coord2, coord2, coord1},
													 {coord2, coord2, coord2} };

		const double Gauss<2>::weights[] = { weight, weight, weight, weight };
	}

	namespace Triangle
	{
		const double Gauss<2>::coord1 = 1.0 / 6.0;
		const double Gauss<2>::coord2 = 2.0 / 3.0;
		const double Gauss<2>::weight = 1.0 / 6.0;

		const LocalCoordinates2D Gauss<2>::localPoints[] = { {coord1, coord1},
													 {coord2, coord1},
													 {coord1, coord2}};

		const double Gauss<2>::weights[] = { weight, weight, weight};
	}
}