
#include "equations.h"
#include "parameters.h"

void Equations::compute_numerical_flux (const Tensor<1,1> & normal,
		 	 	 	 	 	 	 	 	const Vector<double> &uplus,
										const Vector<double> &uminus,
										const double alpha,
										std::array<double, 1> &numerical_flux, Parameters parameters)
{
	switch (parameters.FluxType)
	{
	case (lax_friedrichs): numerical_flux = (1./4.) * (uplus*uplus + uminus*uminus) - 1./2. * std::max(std::abs(uplus),std::abs(uminus)); break;
	case (roe): numerical_flux = 1./4. * (uplus*uplus + uminus*uminus) - 1./2. * std::abs(uplus+uminus) * (uminus - uplus); break;
	case (energy_conserving): numerical_flux = 1./4. * (uplus*uplus + uminus*uminus) - 1./12. * (uminus - uplus)*(uminus-uplus); break;
	}

}
