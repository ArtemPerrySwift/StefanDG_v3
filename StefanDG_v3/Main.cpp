#include "GlobalCalculator.h"
#include "LinearLagrangeBasis.h"
#include "FaceDGCalculator.h"

int main()
{
	GlobalCalculator<LinearLagrangeBasis>::init();
	const double* massMatrix = GlobalCalculator<LinearLagrangeBasis>::getMassMatrix();
	GlobalCalculator<LinearLagrangeBasis>::finalize();
	FaceDGCalculator::init();
	const double (*_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][16] = FaceDGCalculator::getCrossMassMatrixies();
	FaceDGCalculator::finalize();
	return 0;
}