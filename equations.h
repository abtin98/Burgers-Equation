#ifndef EQUATIONS
#define EQUATIONS

#include "util.h"

enum numerical_flux_type
	{
			lax_friedrichs,
			roe,
			energy_conserving
	};

template <int dim>
class Equations
{
public:
	int n_components = 1;

	numerical_flux_type numFlux = energy_conserving;

	void compute_flux_vector(const Vector<double> &U,std::vector<Vector<double>> &flux); //[component][index]

	void compute_numerical_normal_flux(const Tensor<1,dim>             &normal,
            					const double                  &Uplus,
								const double                  &Uminus,
								double		                    &normal_flux);

	void compute_u_matrix(const Vector<double> &U,
						        FullMatrix<double> UMatrix);
private:

};

#endif
