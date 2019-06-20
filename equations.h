#ifndef EQUATIONS
#define EQUATIONS

#include "util.h"

enum numerical_flux_type
	{
			lax_friedrichs,
			roe
	};

template <int dim>
class Equations
{
public:
	int n_components = 1; //Burgers' Equation has only one component

	numerical_flux_type numFlux = lax_friedrichs;

	void compute_flux_vector(const Vector<double> &U,std::vector<Vector<double>> &flux); //Given a Vector of U's, this calculates f(U)

	void compute_numerical_normal_flux(const Tensor<1,dim>             &normal,
            				   const double                    &Uplus,
					   const double            	   &Uminus,
					   double		           &normal_flux); //Gives the numerical flux at cell boundaries

	void compute_u_matrix(const Vector<double> &U,
		       	            FullMatrix<double> UMatrix); //Given U, this computes diag(U)
private:

};

#endif
