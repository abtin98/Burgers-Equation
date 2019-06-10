#ifndef PROBLEM
#define PROBLEM

#include "util.h"
#include "equations.h"
#include "parameters.h"

template <int dim>
class Problem
{
public:
	Problem(Parameters parameters);
	void run();

	const UpdateFlags update_flags = update_values | update_gradients | update_q_points | update_JxW_values,
					  face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors,
					  neighbour_face_update_flags = update_values;

private:
	Triangulation<dim> triangulation;
	const MappingQ1<dim> mapping;

	FE_DGQ<dim> fe;
	Vector<double> current_solution;
	Vector<double> old_solution;
	Vector<double> global_rhs;

	//we will compute these at the beginning and store them for later use.
	std::vector < FullMatrix<double> > inverse_mass_matrix;
	std::vector < std::vector < FullMatrix <double> > > stiffness_matrix; //a vector of n_active cells of dim stiffness matrices.


	DoFHandler<dim> dof_handler;

	const QGaussLobatto<dim> quadrature;
	const QGaussLobatto<dim-1> face_quadrature;

	Equations<dim> burgers_equation;

	FEValues<dim> fe_values;
	FEFaceValues<dim> fe_face_values;
	FEFaceValues<dim> fe_neighbour_face_values;

	Parameters class_parameters;

	//double compute_energy(const Vector<double> &u);

	void initialize_system();
	void assemble_grid();
	void compute_stiffness_and_inverse_mass_matrix();
	void compute_rhs_vector();
	void assemble_cell_term(int cell_index, const std::vector<types::global_dof_index> &dof_indices);
	void assemble_face_term(int 										cell_index,
							const unsigned int 						    face_no,
							const std::vector<types::global_dof_index> &dof_indices,
							const std::vector<types::global_dof_index> &neighbour_dof_indices);
	void perform_runge_kutta_45();
	void output_data(int n_iteration);
	void diagonalize_U (Vector<double> U, FullMatrix<double> &UMatrix);

	//debug functions
	void print_inverse_mass_matrix();
	void print_stiffness_matrix();
};

#endif
