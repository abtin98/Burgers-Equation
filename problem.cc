#include "problem.h"
#include "util.h"

template <int dim>
Problem<dim>::Problem()
{

}

template<int dim>
void Problem<dim>::run()
{
	assemble_grid();
	assemble_system();
	perform_runge_kutta_45();
}

template<int dim>
void Problem<dim>::assemble_grid()
{

}

template<int dim>
void Problem<dim>::assemble_mass_and_stiffness_matrices(FullMatrix<double> &M, std::array<FullMatrix,dim> &S)
{

}

template<int dim>
void Problem<dim>::compute_inverse_mass_matrix(const FullMatrix<double> &M, FullMatrix<double> &M_inv)
{
	M_inv =
}

template<int dim>
void Problem<dim>::compute_rhs_vector()
{

}

template<int dim>
void Problem<dim>::assemble_system()
{

}

template<int dim>
void Problem<dim>::assemble_cell_term()
{

}

template<int dim>
void Problem<dim>::assemble_face_term()
{

}

template<int dim>
void Problem<dim>::perform_runge_kutta_45()
{

}
