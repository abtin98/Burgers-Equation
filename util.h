#ifndef UTIL
#define UTIL

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>
#include <cmath>
using namespace dealii;

Vector<double> RK4 (Vector<double> rhs);
void invert_mass_matrix (const FullMatrix<double> &M, FullMatrix<double> &Minv) //compute inverse of a diagonal mass matrix
{
	for (unsigned int i = 0; i < M.n_rows(); ++i)
	{
		Minv(i,i) = 1./M(i,i);
	}
}

#endif
