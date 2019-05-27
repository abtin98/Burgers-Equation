
#include "equations.h"

void Equations::compute_numerical_flux (const Tensor<1,1> & normal,
		 	 	 	 	 	 	 	 	const Vector<double> &uplus,
										const Vector<double> &uminus,
										const double alpha,
										std::array<double, 1> &numerical_flux,
										Parameters parameters)
{
	switch (parameters.fluxType)
	{
		case lax_friedrichs:
			{
			numerical_flux = (1./4.) * (uplus*uplus + uminus*uminus); //- 1./2. ;//* std::max(std::abs(uplus),std::abs(uminus));
			break;
			}
		case roe:
		{
			numerical_flux = 1./4. * (uplus*uplus + uminus*uminus); //- 1./2. ;//* std::abs(uplus+uminus) * (uminus - uplus);
			break;
		}
		case energy_conserving:
		{
			numerical_flux = 1./4. * (uplus*uplus + uminus*uminus) - 1./12. * (uminus - uplus)*(uminus-uplus);
			break;
		}
	}

}

void Equations::compute_flux(const Vector<double> &u, Vector<double> &flux)
{
	flux = u*u/2.;
}
