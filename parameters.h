#ifndef PARAMETERS
#define PARAMETERS

#include "util.h"

struct Parameters //this could be templated but it's not really necessary
{
public:

	double initial_time{0.0}, final_time{3.0}, delta_t{0.0001};

	 int n_refinements{5}; //number of times the grid needs to be refined. if 5, the grid is divided in half 5 times, resulting in 32 active cells.

	 double left{0.},right{2.}; //boundaries of grid. 

	int polynomial_order_dg{7};

	int quadrature_degree{8}; //number of quadrature points needs to match that of the nodes for collocation

	int face_quadrature_degree{8}; //in 1D, this number doesn't really matter since the faces reduce down to single points.
};

#endif
