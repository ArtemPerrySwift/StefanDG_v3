#include "ElementDGCalculator.h"
#include "LinearLagrangeBasis.h"
#include "FaceDGCalculator.h"

int main()
{
	ElementDGCalculator<LinearLagrangeBasis>::init();
	const double* massMatrix = ElementDGCalculator<LinearLagrangeBasis>::getMassMatrix();
	ElementDGCalculator<LinearLagrangeBasis>::finalize();
	FaceDGCalculator::init();
	const double (*_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][16] = FaceDGCalculator::getCrossMassMatrixies();
	FaceDGCalculator::finalize();
	return 0;
}