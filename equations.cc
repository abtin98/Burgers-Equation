#include "equations.h"

template <int dim>
void Equations<dim>::compute_flux_vector(const Vector<double> &U,std::vector<Vector<double>> &flux)
{
	for (unsigned int component = 0; component < dim; ++component)
	{
		for(unsigned int i = 0; i< U.size(); ++i)
		{
			flux[component][i] = U[i]*U[i]/2.;
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
		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2. * (std::fabs(Uplus) > std::fabs(Uminus) ? std::fabs(Uplus) : std::fabs(Uminus)) * (Uminus - Uplus)*normal[i]) * normal[i]; //check signs
			//1./2. * (Uplus + Uminus) * normal[i];// + 1./2. * (Uplus - Uminus) * normal[i];
		}
		break;
	case roe:
		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./2. * std::fabs(Uplus+Uminus) * (Uminus - Uplus) * normal[i]) * normal[i]; //check signs
		}
		break;
	case energy_conserving:

		for (unsigned int i = 0; i < dim; ++i)
		{
			normal_flux += (1./4. * ((Uplus * Uplus) + (Uminus * Uminus)) - 1./12. * (Uplus - Uminus)*(Uplus - Uminus)) * normal[i];
		}
		break;

	case osher:
		for (unsigned int i = 0; i < dim; ++i)
		{
			if (Uplus < 0 && Uminus < 0)
				normal_flux += Uplus*Uplus/2. * normal[i];
			else if (Uplus > 0 && Uminus > 0)
				normal_flux += Uminus*Uminus/2. * normal[i];
			else if (Uplus <= 0 && Uminus >= 0)
				normal_flux += (Uplus*Uplus + Uminus*Uminus)/2. ;
			else if (Uplus >= 0 && Uminus <= 0 )
				normal_flux += 0.0 * normal[i];

		}
		break;
	}
}

template class Equations<1>;
template class Equations<2>;
