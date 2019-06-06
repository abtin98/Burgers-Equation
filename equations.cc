#include "equations.h"

template <int dim>
void Equations<dim>::compute_flux_vector(const double &U,Vector<double> &flux)
{
	for (unsigned int i = 0; i < dim; ++i)
	{
		flux[i] = U*U/2.;
	}
}

template <int dim>
void Equations<dim>::compute_numerical_normal_flux(const Tensor<1,dim>          &normal,
												   const double                 &Uplus,
												   const double                	&Uminus,
												   double 			         	&normal_flux)
{
	normal_flux = 0;
	switch (numFlux)
	{
	case lax_friedrichs:
		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2 * (Uplus > Uminus ? Uplus : Uminus) * (Uplus - Uminus)) * normal[i]; //check signs
		}
		break;
	case roe:
		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2 * std::fabs(Uplus+Uminus) * (Uplus - Uminus)) * normal[i]; //check signs
		}
		break;
	case energy_conserving:

		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./12 * (Uplus - Uminus)*(Uplus - Uminus)) * normal[i];
		}
		break;
	}
}

template class Equations<1>;
