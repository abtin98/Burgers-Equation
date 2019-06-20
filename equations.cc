#include "equations.h"

template <int dim>
void Equations<dim>::compute_flux_vector(const Vector<double> &U,std::vector<Vector<double>> &flux)
{
	for (unsigned int component = 0; component < dim; ++component)
	{
		for(unsigned int i = 0; i< U.size(); ++i)
		{
			flux[component][i] = U[i]*U[i]/2.; //In case of 1D Burgers', this turns into a simple vector.
		}

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
		for (unsigned int i = 0; i < dim; ++i) //This loop basically does the dot product of the numerical flux and the normal.
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2. * (std::fabs(Uplus) > std::fabs(Uminus) ? std::fabs(Uplus) : std::fabs(Uminus)) * (Uminus - Uplus)*normal[i]) * normal[i];
		}
		break;
	case roe:
		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2. * std::fabs(Uplus+Uminus) * (Uminus - Uplus) * normal[i]) * normal[i];
		}
		break;
	} //more numerical fluxes can be added in the switch statement
}

template class Equations<1>;
template class Equations<2>;
