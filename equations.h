#ifndef EQUATIONS
#define EQUATIONS

#include "util.h"
#include "parameters.h"

class Equations
{
public:
	void compute_numerical_flux (const Tensor<1,1> & normal,
								 const Vector<double> &uplus,
								 const Vector<double> &uminus,
								 const double alpha,
								 std::array<double, 1> &normal_flux,
								 Parameters parameters);
	void compute_flux (const Vector<double> &u, Vector<double> &flux);

private:

};
#endif
