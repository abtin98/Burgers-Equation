#include "util.h"

enum FluxType {lax_friedrichs = 1, roe = 2, energy_conserving = 3};

class Parameters
{
public:
	double dt;
	double finalTime;
	double alpha;
	FluxType fluxType;
	Parameters();
	Parameters(double dt, double finalTime, double alpha, FluxType fluxType);
};
