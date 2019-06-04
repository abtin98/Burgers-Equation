#ifndef PROBLEM
#define PROBLEM

#include "util.h"
#include "equations.h"

template <int dim>
class Problem
{
public:
	Triangulation<dim> triangulation;
	DoFHandler<dim> dof_handler;

	Problem();
	void run();
private:
	void assemble_grid();
	void assemble_mass_and_stiffness_matrices(FullMatrix<double> &M, std::array<FullMatrix, dim> &S);
	void compute_inverse_mass_matrix(const FullMatrix<double> &M, FullMatrix<double> &M_inv);
	void compute_rhs_vector();
	void assemble_system();
	void assemble_cell_term();
	void assemble_face_term();
	void perform_runge_kutta_45();
};

#endif
