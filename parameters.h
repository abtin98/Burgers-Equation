#include "util.h"

class Parameters
{
public:
	double dt;
	double finalTime;
	double alpha;
	enum FluxType {lax_friedrichs, roe, energy_conserving};
	Parameters();
	Parameters(double dt, double finalTime, double alpha, FluxType fluxType);
};
