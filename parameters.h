#ifndef PARAMETERS
#define PARAMETERS

#include "util.h"

struct Parameters
{
public:

	double initial_time{0.0}, final_time{3.0}, delta_t{0.001};

	double delta_x{0.1};

	 int n_refinements{3};

	 double left{0.0};
	 double right{2.0};

	int polynomial_order_dg{4};

	int quadrature_degree{8};

	int face_quadrature_degree{8};
};

#endif
