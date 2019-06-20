#include "problem.h"
#include "util.h"

template <int dim>
Problem<dim>::Problem(Parameters parameters)
:
mapping(),
fe (parameters.polynomial_order_dg),
dof_handler(triangulation),
quadrature(parameters.quadrature_degree),
face_quadrature(parameters.face_quadrature_degree),
fe_values(mapping, fe, quadrature, update_flags),
fe_face_values(mapping, fe, face_quadrature, face_update_flags),
fe_neighbour_face_values(mapping, fe, face_quadrature, neighbour_face_update_flags)
{
	class_parameters = parameters;
}

template<int dim>
void Problem<dim>::run() //this function puts everything together.
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
	std::string variables = "x";
	std::map<std::string,double> constants;
	constants["pi"] = numbers::PI;
	std::string expression = "sin(pi*(x)) + 0.01"; //same function that gassner and ranocha use as well
	initial_condition.initialize(variables,
	              	  	  	  	 expression,
								 constants);
	current_solution.reinit(dof_handler.n_dofs()); //resize the solution vectors to the number of degrees of freedom.
	old_solution.reinit(dof_handler.n_dofs());
	global_rhs.reinit(dof_handler.n_dofs());
	VectorTools::interpolate(mapping, dof_handler, initial_condition, current_solution);
	old_solution = current_solution;
	std::cout << "initial condition implemented" << std::endl;
}

template<int dim>
void Problem<dim>::assemble_grid()
{
	GridGenerator::hyper_cube(triangulation,class_parameters.left,class_parameters.right, true); //makes a hypercube of dim-dimensions

	//std::vector<GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator> > matched_pairs; //can do either Triangulation or DoFHandler
	//GridTools::collect_periodic_faces(triangulation, 0, 1, 0, matched_pairs );
	//GridTools::collect_periodic_faces(triangulation, 2, 3, 1, matched_pairs ); // for 2d
	//triangulation.add_periodicity(matched_pairs);

	triangulation.refine_global(class_parameters.n_refinements);
	std::cout << "number of active cells: " << triangulation.n_active_cells() << std::endl;

	dof_handler.distribute_dofs (fe);
}

template <int dim>
void Problem<dim>::compute_stiffness_and_inverse_mass_matrix() //this function computes and stores the needed matrices only once.
{
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_cell_q_points = quadrature.size();
	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end(); //these are for looping through the cells.

	for (; cell != endc; ++cell) //we now loop over the cells, quadrature points, and degrees of freedom to assemble the mass and stiffness matrices cell by cell.
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
						stiffness_matrix_cell[component][i][j]+= (double) fe_values.shape_grad(j, q_index)[component] * fe_values.shape_value(i,q_index) * fe_values.JxW(q_index);
					}
				}
			}
		}
		inverse_mass_matrix[cell_index].invert(mass_matrix);
		stiffness_matrix[cell_index] = stiffness_matrix_cell;
	}
}

template <int dim>
void Problem<dim>::diagonalize_U (Vector<double> &U, FullMatrix<double> &UMatrix)
{
	for (unsigned int i = 0; i < U.size(); i++)
	{
		UMatrix(i,i) = U(i);
	}
}


template<int dim>
void Problem<dim>::compute_rhs_vector()
{
	const unsigned int dofs_per_cell = fe.dofs_per_cell;

	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
	std::vector<types::global_dof_index> dof_indices_neighbour(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	//we loop over the cells, quadrature points and degrees of freedom and assemble the cell and face contributions to the rhs
	for (; cell != endc; ++cell)
	{

		fe_values.reinit(cell);
		cell->get_dof_indices (dof_indices); //here the dof_indices come from fe_values
		int cell_index = cell->index();
		assemble_cell_term(cell_index, dof_indices);
		for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
		{
			if (!cell->face(face_no)->at_boundary()) //internal face
			{
			fe_face_values.reinit(cell, face_no);

			typename DoFHandler<dim>::active_cell_iterator neighbour_cell = cell->neighbor(face_no);
			neighbour_cell->get_dof_indices(dof_indices_neighbour);
			int neighbour_cell_index = neighbour_cell->index();
			fe_neighbour_face_values.reinit(neighbour_cell,(face_no == 1) ? 0 : 1);
			assemble_face_term(cell_index, face_no,dof_indices, dof_indices_neighbour);
			}

			else //if we're dealing with boundaries, we set the periodic boundary conditions
			{
				if (cell_index == 0 && face_no == 0)
				{
					fe_face_values.reinit(cell, face_no);
					typename DoFHandler<dim>::active_cell_iterator neighbour_cell = dof_handler.begin_active();
					for (unsigned int i = 0 ; i < std::pow(2,class_parameters.n_refinements) - 1; ++i)
					{
						++neighbour_cell;
					}
					neighbour_cell->get_dof_indices(dof_indices_neighbour);
					int neighbour_cell_index = neighbour_cell->index();
					fe_neighbour_face_values.reinit(neighbour_cell,(face_no == 1) ? 0 : 1);

				}
				else if (cell_index == std::pow(2,class_parameters.n_refinements) - 1 && face_no == 1)
				{
					fe_face_values.reinit(cell, face_no);
					typename DoFHandler<dim>::active_cell_iterator neighbour_cell = dof_handler.begin_active();
					neighbour_cell->get_dof_indices(dof_indices_neighbour);
					int neighbour_cell_index = neighbour_cell->index();
					fe_neighbour_face_values.reinit(neighbour_cell,(face_no == 1) ? 0 : 1); 
				}

				assemble_face_term(cell_index,face_no,dof_indices,dof_indices_neighbour);
			}
		}
	}
}

template<int dim>
void Problem<dim>::assemble_cell_term(int cell_index, const std::vector<types::global_dof_index> &dof_indices) //This function calculates the volume integrals in DG
	//in other words, it calculates -1/3 * S^T * f(u) + 1/6 * (U * (S - S^T) * u), which is the weak form of Burgers' split form. 
{
	const unsigned int n_q_points = fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
	Vector<double> cell_rhs_intermediate (dofs_per_cell), cell_rhs (dofs_per_cell);
	cell_rhs = 0;
	cell_rhs_intermediate = 0;
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

	FullMatrix<double> UMatrix (dofs_per_cell,dofs_per_cell);
	UMatrix = 0;

	Vector<double> int1 (dofs_per_cell), int2 (dofs_per_cell);
	int1 = 0;
	int2 = 0;

	diagonalize_U(U,UMatrix);

	for (unsigned int component = 0; component < dim; ++component)
	{
		stiffness_matrix[cell_index][component].vmult(int1,U, true);
		stiffness_matrix[cell_index][component].Tvmult(int2, U, true);
	}

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		int1(i) = int1(i) -  int2(i);
	}


	UMatrix.vmult(cell_rhs_intermediate, int1);


	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		cell_rhs(i) += cell_rhs_intermediate(i)/6.0 ;
	}

	cell_rhs_intermediate = 0;

	burgers_equation.compute_flux_vector(U,flux_vector);

	for (int component = 0; component < dim; ++component)
	{
		stiffness_matrix[cell_index][component].Tvmult(cell_rhs_intermediate, flux_vector[component], true); //S_x * F_x + S_y * F_y + S_z + F_z
	}

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		cell_rhs(i) -= (cell_rhs_intermediate(i))* 2.0/3.0; //second split form of flux added to rhs.
	}

	Vector <double> rhs (dofs_per_cell);
	rhs = 0;
	inverse_mass_matrix[cell_index].vmult(rhs, cell_rhs, true);

	for (unsigned int i = 0; i < dofs_per_cell; i++)
	{
		global_rhs(dof_indices[i]) -= rhs(i); //adds everything to global rhs_vector to prepare for time integration.
	}
}

template<int dim>
void Problem<dim>::assemble_face_term(int cell_index,
									  const unsigned int          face_no,
									  const std::vector<types::global_dof_index> &dof_indices,
									  const std::vector<types::global_dof_index> &neighbour_dof_indices )
       //this function calculates the surface integrals involving the numerical flux.	
{

	const unsigned int n_q_points = fe_face_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_face_values.dofs_per_cell;

	Vector<double> cell_rhs_intermediate (dofs_per_cell), cell_rhs (dofs_per_cell);
	cell_rhs = 0;
	cell_rhs_intermediate = 0;
	std::vector < Vector < double > > flux_vector;

	flux_vector.resize(dim, Vector<double> (dofs_per_cell));

	double normal_flux;

	std::vector<double> independent_local_dof_values(dofs_per_cell), independent_neighbour_dof_values(dofs_per_cell),
						Uplus(n_q_points), Uminus(n_q_points);

	Tensor<1,dim> normal;

	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		independent_local_dof_values[i] = current_solution(dof_indices[i]); 
		independent_neighbour_dof_values[i] = current_solution(neighbour_dof_indices[i]);
	}

	for (unsigned int q = 0; q < n_q_points; ++q)
	{
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			Uplus[q] += independent_local_dof_values[i] * fe_face_values.shape_value(i,q) ;
			Uminus[q] += independent_neighbour_dof_values[i] * fe_neighbour_face_values.shape_value(i,q);
		}
	}

	for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
	{
		for (unsigned int i = 0 ; i < dofs_per_cell; ++i)
		{
			normal_flux = 0;
			burgers_equation.compute_numerical_normal_flux(fe_face_values.normal_vector(q_index),Uplus[q_index],Uminus[q_index],normal_flux); 

			cell_rhs_intermediate(i) += normal_flux * fe_face_values.shape_value(i,q_index) * fe_face_values.JxW(q_index);
			//added Uplus in numerical flux for strong form.
		}
	}

	inverse_mass_matrix[cell_index].vmult(cell_rhs, cell_rhs_intermediate);
	for (unsigned int i = 0; i < dofs_per_cell; ++i)
	{
		global_rhs(dof_indices[i]) -= cell_rhs(i);
	}
}

template<int dim>
void Problem<dim>::perform_runge_kutta_45() //this function does the time integration and outputs solution files and an energy plot at the end.
{
	std::ofstream myfile ("energy_llf.gpl" , std::ios::trunc);

	for (int n_iteration = 0; n_iteration < (class_parameters.final_time - class_parameters.initial_time)/class_parameters.delta_t ; ++n_iteration)
	{
//		std::cout << "---------------------" << std::endl
//				  << "current iteration #" << n_iteration << std::endl
//				  << "number of active cells: " << triangulation.n_active_cells() << std::endl
//				  << "number of total degrees of freedom: " << dof_handler.n_dofs() << std:: endl
//				  << "dofs per cell: " << fe_values.dofs_per_cell << std::endl
//				  << "--------------------" << std::endl;

		if (n_iteration % 100 == 0) //output data every 100 iterations. this number can be changed
					output_data(n_iteration);
		global_rhs = 0;
		compute_rhs_vector();
		Vector<double> v1 (dof_handler.n_dofs()), v2 (dof_handler.n_dofs()), u(dof_handler.n_dofs());
		u = current_solution;
		//old_solution = current_solution;
		v1 = 0;
		v1 += u;
		global_rhs *= class_parameters.delta_t;
		v1 += global_rhs;

		current_solution = v1;
		global_rhs = 0;



		compute_rhs_vector();

		for (unsigned int i = 0; i<dof_handler.n_dofs(); ++i)
		{
			current_solution(i) = 1./2. * (u(i) + v1(i) + class_parameters.delta_t * global_rhs(i));
		}

		old_solution = current_solution;

		double energy = compute_energy();

		myfile << n_iteration * class_parameters.delta_t << " " << energy << std::endl;
	}
	myfile.close();
}

template<int dim>
void Problem<dim>::output_data(int n_iteration)
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector(current_solution, "solution");
	data_out.build_patches();
	std::ofstream output ("solution-" + Utilities::int_to_string(n_iteration,4) +".gpl");
	data_out.write_gnuplot (output);

//	std::cout << "iteration " << n_iteration << std::endl;
//	std::cout << "----------" << std::endl;
//	for (unsigned int i = 0 ; i< dof_handler.n_dofs(); ++i)
//	{
//		std::cout << current_solution(i) << std::endl;
//	}
}


//debug functions from here on.

template <int dim>
void Problem<dim>::print_inverse_mass_matrix()
{
	for (unsigned int cell_no = 0; cell_no < triangulation.n_active_cells(); ++cell_no)
	{
		std::cout << "cell number " << cell_no << std::endl;
		for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
			{
				for (unsigned int j = 0 ; j < fe_values.dofs_per_cell; ++j)
				{
					std::cout << inverse_mass_matrix[cell_no][i][j] << " ";
				}
				std::cout << std::endl;
			}
		std::cout << "----------" << std::endl;
	}

}

template <int dim>
void Problem<dim>::print_stiffness_matrix()
{
	for (unsigned int cell_no = 0; cell_no < triangulation.n_active_cells(); ++cell_no)
		{
			std::cout << "cell number " << cell_no << std::endl;
			for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
				{
					for (unsigned int j = 0 ; j < fe_values.dofs_per_cell; ++j)
					{
						std::cout << stiffness_matrix[cell_no][0][i][j] << " ";
					}
					std::cout << std::endl;
				}
			std::cout << "----------" << std::endl;
		}
}

template <int dim>
void Problem<dim>::print_diff_matrix()
{
	FullMatrix<double> D(fe_values.dofs_per_cell,fe_values.dofs_per_cell);
	inverse_mass_matrix[0].mmult(D,stiffness_matrix[0][0]);
	for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
					{
						for (unsigned int j = 0 ; j < fe_values.dofs_per_cell; ++j)
						{
							std::cout << D(i,j) << " ";
						}
						std::cout << std::endl;
					}
				std::cout << "----------" << std::endl;

}

template <int dim>
double Problem<dim>::compute_energy() //need to loop over cells, current function won't work.
{
	const unsigned int n_q_points = fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_values.dofs_per_cell;

	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
	double E = 0;

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	for (; cell != endc; ++cell)
	{
		fe_values.reinit(cell);
		cell->get_dof_indices (dof_indices);
		int cell_index = cell->index();
		std::vector<double> independent_local_dof_values;
			independent_local_dof_values.resize(dofs_per_cell);
			Vector<double> U (dofs_per_cell);
			U = 0;

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
		for (unsigned int i = 0 ; i< dofs_per_cell; ++i)
		{
			E += U(i) * U(i) * 1./inverse_mass_matrix[cell_index][i][i];
		}
	}


	return E;
}

template class Problem<1>;
template class Problem<2>;
