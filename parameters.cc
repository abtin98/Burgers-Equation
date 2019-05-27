#include "parameters.h"

Parameters::Parameters()
{
	this->dt = 0.01;
	this->alpha = 0.5;
	this->finalTime = 3;
	this->fluxType = lax_friedrichs;
}

Parameters::Parameters(double dt, double finalTime, double alpha, FluxType fluxType )
{
	this->dt = dt;
	this->finalTime = finalTime;
	this->alpha = alpha;
	this->fluxType = fluxType;
}
