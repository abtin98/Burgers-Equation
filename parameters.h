#ifndef PARAMETERS
#define PARAMETERS

#include "util.h"

class Parameters
{
public:
	double initial_time, final_time, delta_t;

	double delta_x;

	const int n_refinements = 10;

	const double left = 0.0;
	const double right = 10.0;

	enum quadrature_kind
	{
		gauss,
		gauss_lobatto
	};

	int polynomial_order_dg = 7;

	int quadrature_degree = 8;

	int face_quadrature_degree = 8;
};

#endif
