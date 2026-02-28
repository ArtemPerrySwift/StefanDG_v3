#include "ElementDGCalculator.h"
#include "LinearLagrangeBasis.h"
#include "FaceDGCalculator.h"

int main()
{
	ElementDGCalculator<LinearLagrangeBasis>::init();
	const double* massMatrix = ElementDGCalculator<LinearLagrangeBasis>::getMassMatrix();
	ElementDGCalculator<LinearLagrangeBasis>::finalize();

	FaceDGCalculator<LinearLagrangeBasis>::init();
	const double (*_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][16] = FaceDGCalculator<LinearLagrangeBasis>::getCrossMassMatrixies();
	FaceDGCalculator<LinearLagrangeBasis>::finalize();

	return 0;
}