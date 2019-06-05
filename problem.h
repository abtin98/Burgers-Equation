#ifndef PROBLEM
#define PROBLEM

#include "util.h"
#include "equations.h"

template <int dim>
class Problem
{
public:
	Problem();
	void run();
private:
	Triangulation<dim> triangulation;
	const MappingQ<dim> mapping;

	FE_DGQ<dim> fe;
	Vector<double> current_solution;
	Vector<double> old_solution;
	Vector<double> rhs;
	FullMatrix<double> mass_matrix;
	FullMatrix<double> inverse_mass_matrix;
	std::array<FullMatrix, dim> stiffness_matrix;
	std::array<Vector, dim> flux_vector;

	DoFHandler<dim> dof_handler;

	const QGaussLobatto<dim> quadrature;
	const QGaussLobatto<dim-1> face_quadrature;

	Parameters parameters;
	Equations<dim> burgers_equation;


	std::string variables = "x";
	std::map<std::string,double> constants;
	constants["pi"] = 3.14;
	std::string experession = "sin(2*pi*x) + 0.01";
	FunctionParser<dim> initial_condition(1);
	initial_condition.initialize(variables,
	              	  	  	  	 expression,
								 constants);

	double compute_energy(const Vector<double> &u);

	void assemble_grid();
	void compute_inverse_mass_matrix(const FullMatrix<double> &M, FullMatrix<double> &M_inv);
	void compute_rhs_vector();
	void assemble_system();
	void assemble_cell_term(const FEValues<dim> &fe_values, const std::vector<types::global_dof_index> &dofs_indices);
	void assemble_face_term(const unsigned int               face_no,
            				const FEFaceValues<dim>     &fe_v,
							const FEFaceValues<dim>     &fe_v_neighbor,
							const std::vector<types::global_dof_index> &dof_indices,
							const std::vector<types::global_dof_index> &neighbour_dof_indices);
	void perform_runge_kutta_45();
};

#endif
