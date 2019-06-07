#include "problem.h"
#include "util.h"

template <int dim>
Problem<dim>::Problem(Parameters parameters)
:
//parameters(),
mapping(),
fe (parameters.polynomial_order_dg),
dof_handler(triangulation),
quadrature(parameters.quadrature_degree),
face_quadrature(parameters.face_quadrature_degree),
fe_values(mapping, fe, quadrature, update_flags),
fe_face_values(mapping, fe, face_quadrature, face_update_flags),
fe_neighbour_face_values(mapping, fe, face_quadrature, neighbour_face_update_flags)
{
	//fe_values = new FEValues<dim> (mapping, fe, quadrature, update_flags);
	//fe_face_values = new FEFaceValues<dim> (mapping, fe, face_quadrature, face_update_flags);
	//fe_neighbour_face_values = new FEFaceValues<dim> (mapping, fe, face_quadrature, neighbour_face_update_flags);
	class_parameters = parameters;
}

template<int dim>
void Problem<dim>::run()
{

	assemble_grid();
	initialize_system();
	compute_stiffness_and_inverse_mass_matrix();
	perform_runge_kutta_45();
}

template <int dim>
void Problem<dim>::initialize_system()
{
	//initial conditions
	FunctionParser<dim> initial_condition;
	std::string variables = "x,y";
	std::map<std::string,double> constants;
	constants["pi"] = numbers::PI;
	std::string expression = "sin(pi*x) + 0.01";
	initial_condition.initialize(variables,
	              	  	  	  	 expression,
								 constants);
	current_solution.reinit(dof_handler.n_dofs());
	old_solution.reinit(dof_handler.n_dofs());
	global_rhs.reinit(dof_handler.n_dofs());
	VectorTools::interpolate(mapping, dof_handler, initial_condition, current_solution);
	old_solution = current_solution;
	std::cout << "initial condition implemented" << std::endl;
}

template<int dim>
void Problem<dim>::assemble_grid()
{
	GridGenerator::hyper_cube(triangulation,class_parameters.left,class_parameters.right, true);

	std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator> > matched_pairs; //can do either Triangulation or DoFHandler
	GridTools::collect_periodic_faces(triangulation, 0, 1, 0, matched_pairs );
	GridTools::collect_periodic_faces(triangulation, 2, 3, 1, matched_pairs ); // for 2d
	triangulation.add_periodicity(matched_pairs);
		//DoFTools::make_periodicity_constraints <DoFHandler<dim>> (matched_pairs, );

	triangulation.refine_global(class_parameters.n_refinements);
	std::cout << "number of active cells: " << triangulation.n_active_cells() << std::endl;

	dof_handler.distribute_dofs (fe);
}

template <int dim>
void Problem<dim>::compute_stiffness_and_inverse_mass_matrix()
{
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_cell_q_points = quadrature.size();
	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

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
		stiffness_matrix.resize(triangulation.n_active_cells(), std::vector<FullMatrix<double>>(dim, FullMatrix<double>(dofs_per_cell)));
		for (unsigned int q_index = 0; q_index < n_cell_q_points; ++q_index)
		{
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
				{
					mass_matrix(i,j) += fe_values.shape_value (i, q_index) * fe_values.shape_value (j, q_index) * fe_values.JxW(q_index);
					for (int component = 0; component < dim; ++component)
					{
						stiffness_matrix_cell[component][i][j]+= (double) fe_values.shape_grad(i, q_index)[component] * fe_values.shape_value(j,q_index) * fe_values.JxW(q_index);
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

	//const unsigned int n_cell_q_points = quadrature.size();
	//const unsigned int n_face_q_points = face_quadrature.size();
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

		assemble_cell_term(cell_index, dof_indices);
		for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
		{
			if (!cell->face(face_no)->at_boundary())
			{
			fe_face_values.reinit(cell, face_no); //what about the neighbour cell? does it need to be reinit'ed?
			typename DoFHandler<dim>::active_cell_iterator neighbour_cell = cell->neighbor(face_no); //not sure if this variable exists
			neighbour_cell->get_dof_indices(dof_indices_neighbour);
			assemble_face_term(cell_index, face_no,dof_indices, dof_indices_neighbour); // not sure what arguments to pass
			//if (cell->at_boundary(face_no))
			//{
				// this is not necessary because it'll be periodic.
			//}
			}

			else
			{
				std::cout << "face " << face_no << " of cell " << cell_index << " is at boundary" << std::endl;
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
	independent_local_dof_values.resize(dofs_per_cell);
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

	for (unsigned int i = 0; i < dofs_per_cell; i++)
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

	Vector<double> Uplus(dofs_per_cell), Uminus(dofs_per_cell);

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
			cell_rhs_intermediate(i) += normal_flux * fe_face_values.shape_value(i,q_index) * fe_face_values.JxW(q_index);
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

	for (int n_iteration = 0; n_iteration < (class_parameters.final_time - class_parameters.initial_time)/class_parameters.delta_t ; ++n_iteration)
	{
		std::cout << "---------------------" << std::endl
				  << "current iteration #" << n_iteration << std::endl
				  << "number of active cells: " << triangulation.n_active_cells() << std::endl
				  << "number of total degrees of freedom: " << dof_handler.n_dofs() << std:: endl
				  << "dofs per cell: " << fe_values.dofs_per_cell << std::endl
				  << "--------------------" << std::endl;

		compute_rhs_vector();
		old_solution = current_solution;

		//for (unsigned int i = 0; i < global_rhs.size(); i++)
			//global_rhs(i) = global_rhs(i) * parameters.delta_t;

		global_rhs *= class_parameters.delta_t;
		old_solution += global_rhs;
		current_solution = old_solution;
		if (n_iteration < 30)
			output_data(n_iteration);
		//current_solution = old_solution + global_rhs;
	}
}

template<int dim>
void Problem<dim>::output_data(int n_iteration)
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector(current_solution, "solution");
	data_out.build_patches();
	std::ofstream output ("solution-" + Utilities::int_to_string(n_iteration,4) +".vtk");
	data_out.write_vtk (output);
}

//template <int dim>
//double Problem<dim>::compute_energy(const Vector<double> &u) //need to loop over cells, current function won't work.
//{
//	Vector<double> product;
//	double energy {0};
//	//mass_matrix.vmult(product, u);
//
//	for (int i = 0; i < u.size(); ++i)
//	{
//		energy += product[i]*u[i];
//	}
//	return energy;
//}

template class Problem<1>;
template class Problem<2>;
