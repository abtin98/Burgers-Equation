#include "problem.h"
#include "util.h"

template <int dim>
Problem<dim>::Problem() :   fe (3), dof_handler (triangulation) //fe(1) indicates polynomial of degree 1.
{
}

template <int dim>
void Problem<dim>::make_grid(int n_refinements)
{
	GridGenerator::hyper_cube (triangulation,0,1);
	triangulation.refine_global (n_refinements);
	std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}

template <int dim>
void Problem<dim>::setup_system()
{
	dof_handler.distribute_dofs(fe);
	std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler,dsp);
	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void Problem<dim>::assemble_system()
{
	QGaussLobatto<dim> quadrature_formula(4);
	QGaussLobatto<dim-1> face_quadrature_formula(2);

	const UpdateFlags update_flags = update_values
									| update_gradients
									| update_JxW_values,
									face_update_flags = update_values,
									neighbor_face_update_flags = update_values;

	FEValues<dim> fe_values (fe, quadrature_formula, update_flags);
	FEFaceValues<dim> fe_face_values (fe,face_quadrature_formula,face_update_flags);
	FEFaceValues<dim> fe_neighbor_face_values (fe, face_quadrature_formula, neighbor_face_update_flags);
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	FullMatrix<double> cell_mass_matrix(dofs_per_cell,dofs_per_cell);
	FullMatrix<double> cell_stiffness_matrix(dofs_per_cell,dofs_per_cell);
	Vector<double> cell_rhs (dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
	for (const auto &cell: dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		cell_mass_matrix = 0;
		cell_stiffness_matrix = 0;
		cell_rhs = 0;
		for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
		{
			for (unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
				{
					cell_mass_matrix(i,j) += (fe_values.shape_value(i,q_index) * fe_values.shape_value(j,q_index) * fe_values.JxW(q_index));
					//cell_stiffness_matrix(i,j) += (fe_values.shape_value(j,q_index) * fe_values.shape_grad(i,q_index) * fe_values.JxW(q_index));
				}
			}


			for (unsigned int i = 0; i < dofs_per_cell; ++ i)
			{
				cell_rhs(i) += (fe_values.shape_value(i,q_index)*1.5*fe_values.JxW(q_index));
			}
		}

		//print mass matrix
					for (unsigned int i = 0; i < dofs_per_cell; ++i)
						{
							for (unsigned int j = 0; j < dofs_per_cell; ++j)
							{
								std::cout << cell_mass_matrix(i,j) << " ";
							}
							std::cout<<std::endl;
						}

		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			for  (unsigned int j = 0; j < dofs_per_cell; ++j)
			{
				system_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_mass_matrix(i,j));
			}
		}
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}

	}
	std::map<types::global_dof_index,double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler,0,Functions::ZeroFunction<dim>(),boundary_values);
	MatrixTools::apply_boundary_values(boundary_values,system_matrix,solution,system_rhs);
	//system_rhs(dof_handler.n_dofs()-1) = 0;

}

template <int dim>
void Problem<dim>::solve()
{
	SolverControl solver_control (1000, 1e-12);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix,solution,system_rhs,PreconditionIdentity());
	for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
	{
		std::cout << solution(i) << std::endl;
	}
}

template <int dim>
void Problem<dim>::output()
{
	DataOut<1> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");
	data_out.build_patches();
	std::ofstream output ("solution.vtu");
	data_out.write_vtu(output);
}

template <int dim>
void Problem<dim>::run()
{
	int n_refinements = 0;
	make_grid(n_refinements);
	setup_system();
	assemble_system();
	//solve();
	//output();
}
