#include "problem.h"
#include "util.h"

template <int dim>
Problem<dim>::Problem()
:
mapping(),
fe (parameters.polynomial_order_dg),
dof_handler(triangulation),
quadrature(parameters.quadrature_degree),
face_quadrature(parameters.face_quadrature_degree)
{
	const UpdateFlags update_flags = update_values | update_gradients | update_q_points | update_JxW_values,
					  face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors,
					  neighbour_face_update_flags = update_values;

	fe_values = new FEValues<dim> (mapping, fe, quadrature, update_flags);
	fe_face_values = new FEFaceValues<dim> (mapping, fe, face_quadrature, face_update_flags);
	fe_neighbour_face_values = new FEFaceValues<dim> (mapping, fe, face_quadrature, neighbour_face_update_flags);
}

template<int dim>
void Problem<dim>::run()
{

	assemble_grid();
	initialize_system();
	perform_runge_kutta_45();
}

template <int dim>
void Problem<dim>::initialize_system()
{
	//initial conditions
	FunctionParser<dim> initial_condition;
	std::string variables = "x";
	std::map<std::string,double> constants;
	constants["pi"] = numbers::PI;
	std::string expression = "sin(pi*x) + 0.01";
	initial_condition.initialize(variables,
	              	  	  	  	 expression,
								 constants);
	VectorTools::interpolate(mapping, dof_handler, initial_condition, current_solution);
	old_solution = current_solution;

}

template<int dim>
void Problem<dim>::assemble_grid()
{
	GridGenerator::hyper_cube(triangulation,parameters.left,parameters.right);

	std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator> > matched_pairs;
	GridTools::collect_periodic_faces(dof_handler, 0, 1, 0, matched_pairs );
		triangulation.add_periodicity(matched_pairs);
		DoFTools::make_periodicity_constraints <DoFHandler<dim>> (matched_pairs, );

	triangulation.refine_global(parameters.n_refinements);
	std::cout << "number of active cells: " << triangulation.n_active_cells() << std::endl;
	GridTools::collect_periodic_faces(dof_handler, 0, 1, 0, );
	triangulation.add_periodicity();
	DoFTools::make_periodicity_constraints();

	dof_handler.distribute_dofs (fe);
}

template <int dim>
void Problem<dim>::compute_stiffness_and_inverse_mass_matrix()
{
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_cell_q_points = quadrature.size();
	const unsigned int n_face_q_points = face_quadrature.size();
	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
	std::vector<types::global_dof_index> dof_indices_neighbour(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	for (; cell != endc; ++cell)
	{
		fe_values.reinit(cell);
		int cell_index = cell->index();
		cell->get_dof_indices (dof_indices);
		FullMatrix<double> mass_matrix (dofs_per_cell, dofs_per_cell);
		std::vector<FullMatrix<double>> stiffness_matrix_cell;
		stiffness_matrix_cell.resize(dim, FullMatrix<double>(dofs_per_cell));

		inverse_mass_matrix.resize(triangulation.n_active_cells(), FullMatrix<double>(dofs_per_cell));
		stiffness_matrix.resize(triangulation.n_active_cells(), (dim, FullMatrix<double>(dofs_per_cell)));
		for (unsigned int q_index = 0; q_index < n_cell_q_points; ++q_index)
		{
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
				{
					mass_matrix(i,j) += fe_values.shape_value (i, q_index) * fe_values.shape_value (j, q_index) * fe_values.JxW(q_index);
					for (int component = 0; component < dim; ++component)
					{
						stiffness_matrix_cell[component](i,j)+= fe_values.shape_grad_component(i, q_index, component) * fe_values.shape_value(j,q_index) * fe_values.JxW(q_index);
					}
				}
			}
		}
		inverse_mass_matrix[cell_index].invert(mass_matrix);
		stiffness_matrix[cell_index] = stiffness_matrix_cell;
	}
}


template<int dim>
void Problem<dim>::compute_rhs_vector()
{
	//transferred from Problem<dim>::assemble_system. need to pass arguments to assemble_system.
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_cell_q_points = quadrature.size();
	const unsigned int n_face_q_points = face_quadrature.size();
	//assemble_face_term();
	//assemble_cell_term();

	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
	std::vector<types::global_dof_index> dof_indices_neighbour(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	for (; cell != endc; ++cell)
	{
		fe_values.reinit(cell);
		cell->get_dof_indices (dof_indices);
		int cell_index = cell->index();
		assemble_cell_term(cell_index, fe_values, dof_indices); // not sure what arguments to pass...

		for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
		{
			fe_face_values.reinit(cell, face_no); //what about the neighbour cell? does it need to be reinit'ed?
			typename DoFHandler<dim>::active_cell_iterator neighbour_cell = cell->neighbour(face_no); //not sure if this variable exists
			neighbour_cell->get_dof_indices(dof_indices_neighbour);
			assemble_face_term(face_no, fe_face_values, fe_neighbour_face_values,dof_indices, dof_indices_neighbour); // not sure what arguments to pass
			if (cell->at_boundary(face_no))
			{
				// this is not necessary because it'll be periodic.
			}
		}
	}
}

template<int dim>
void Problem<dim>::assemble_cell_term(int cell_index, const std::vector<types::global_dof_index> &dof_indices) //are these local or global dofs?
{
	const unsigned int n_q_points = fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
	Vector<double> cell_rhs_intermediate (dofs_per_cell), cell_rhs (dofs_per_cell);
	cell_rhs = 0;
	std::vector<Vector<double>> flux_vector;
	flux_vector.resize(dim, Vector<double> (dofs_per_cell));

	std::vector<double> independent_local_dof_values;
	Vector<double> U (dofs_per_cell);

	for (unsigned int i = 0; i< dofs_per_cell; ++i)
	{
		independent_local_dof_values[i] = current_solution(dof_indices[i]);
	}

	for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
	{
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			U(i) += independent_local_dof_values[i] * fe_values.shape_value(i,q_index);
		}
	}

	burgers_equation.compute_flux_vector(U,flux_vector);

	for (int component = 0; component < dim; ++component)
	{
		stiffness_matrix[cell_index][component].vmult(cell_rhs_intermediate, flux_vector[component], true); //S_x * F_x + S_y * F_y + S_z + F_z
	}

	inverse_mass_matrix[cell_index].vmult(cell_rhs, cell_rhs_intermediate, true);

	for (int i = 0; i < dofs_per_cell; i++)
	{
		global_rhs(dof_indices[i]) += cell_rhs(i); //adds everything to global rhs_vector to prepare for time integration.
	}
}

template<int dim>
void Problem<dim>::assemble_face_term(int cell_index,
									  const unsigned int          face_no,
									  const std::vector<types::global_dof_index> &dof_indices,
									  const std::vector<types::global_dof_index> &neighbour_dof_indices ) //where does face_no come in????
{

	const unsigned int n_q_points = fe_face_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_face_values.dofs_per_cell;

	Vector<double> cell_rhs_intermediate (dofs_per_cell), cell_rhs (dofs_per_cell);
	cell_rhs = 0;
	cell_rhs_intermediate = 0;
	std::vector < Vector < double > > flux_vector;

	flux_vector.resize(dim, Vector<double> (dofs_per_cell));

	double normal_flux = 0;

	Vector<double> Uplus, Uminus;

	Tensor<1,dim> normal;

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		Uplus[i] = current_solution(dof_indices[i]); //or should it be old solution?
		Uminus[i] = current_solution(neighbour_dof_indices[i]);
	}



	for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
	{
		for (unsigned int i = 0 ; i < dofs_per_cell; ++i)
		{
			burgers_equation.compute_numerical_normal_flux(fe_face_values.normal_vector(q_index),Uplus[q_index],Uminus[q_index],normal_flux); //how to incorporate the arguments here?
			cell_rhs_intermediate(i) += normal_flux * fe_face_values.shape_values(i,q_index) * fe_face_values.JxW(q_index);
		}
	}

	inverse_mass_matrix[cell_index].vmult(cell_rhs, cell_rhs_intermediate);
	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		global_rhs(dof_indices[i]) -= cell_rhs(i);
	}
}


/*
 * Not sure how to do the whole indexing for assembling face and cell terms... should I have one big global matrix + vector?...
 * How to propagate the assemble_face_term info from each face to the cell_rhs vector properly without causing interference?
 */

template<int dim>
void Problem<dim>::perform_runge_kutta_45()
{

	for (int n_iteration = 0; n_iteration < (parameters.final_time - parameters.initial_time)/parameters.delta_t ; ++n_iteration)
	{
		compute_rhs_vector();
		old_solution = current_solution;
		current_solution = old_solution + global_rhs * parameters.delta_t;
	}
}

template <int dim>
double Problem<dim>::compute_energy(const Vector<double> &u) //need to loop over cells, current function won't work.
{
	Vector<double> product;
	double energy {0};
	//mass_matrix.vmult(product, u);

	for (int i = 0; i < u.size(); ++i)
	{
		energy += product[i]*u[i];
	}
	return energy;
}

//template <int dim>
//Problem<dim>::~Problem ()
//{
//std::cout << "Destructing DGBase..." << std::endl;
//delete_fe_values();
//}
//
//template <int dim>
//void Problem<dim>::delete_fe_values ()
//{
//	if (fe_values          != NULL) delete fe_values;
//	if (fe_face_values      != NULL) delete fe_face_values;
//	if (fe_neighbour_face_values  != NULL) delete fe_neighbour_face_values;
//	fe_values          = NULL;
//	fe_face_values      = NULL;
//	fe_neighbour_face_values   = NULL;
//}


template class Problem<1>;
