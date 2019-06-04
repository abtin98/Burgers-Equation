#ifndef PARAMETERS
#define PARAMETERS

#include "util.h"

template <int dim>
struct Parameters
{
	double initial_time, final_time, delta_t;

	Point<dim> corner_a, corner_b;

	double delta_x;

	enum quadrature_kind
	{
		gauss,
		gauss_lobatto
	};

	int polynomial_order_dg;
};

#endif
