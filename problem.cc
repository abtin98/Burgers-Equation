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

}

template<int dim>
void Problem<dim>::run()
{
	assemble_grid();


	perform_runge_kutta_45();
}

template<int dim>
void Problem<dim>::assemble_grid()
{
	GridGenerator::hyper_cube(triangulation,parameters.left,parameters.right);
	triangulation.refine_global(parameters.n_refinements);
	std::cout << "number of active cells: " << triangulation.n_active_cells() << std::endl;
	GridTools::collect_periodic_faces(dof_handler, 0, 1, 0, );
	triangulation.add_periodicity();
	DoFTools::make_periodicity_constraints();

	dof_handler.distribute_dofs (fe);
}

template<int dim>
void Problem<dim>::compute_inverse_mass_matrix(const FullMatrix<double> &M, FullMatrix<double> &M_inv)
{
	M_inv.invert(M);
}

template<int dim>
void Problem<dim>::compute_rhs_vector()
{
	assemble_system();
	Vector<double> sum_of_flux_and_num_flux;
	sum_of_flux_and_num_flux = 0;
	for (int i = 0; i<dim ;++i)
	{
		stiffness_matrix[i].vmult(sum_of_flux_and_num_flux,flux_vector[i],true);
	}

	sum_of_flux_and_num_flux -= rhs;

	inverse_mass_matrix.vmult(rhs, sum_of_flux_and_num_flux);
}

template<int dim>
void Problem<dim>::assemble_system()
{
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_cell_q_points = quadrature.size();
	const unsigned int n_face_q_points = face_quadrature.size();
	//assemble_face_term();
	//assemble_cell_term();

	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
	std::vector<types::global_dof_index> dof_indices_neighbour(dofs_per_cell);

	const UpdateFlags update_flags = update_values | update_gradients | update_q_points | update_JxW_values,
					  face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors,
					  neighbour_face_update_flags = update_values;
	FEValues<dim> fe_values (mapping, fe, quadrature, update_flags);
	FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature, face_update_flags);
	FEFaceValues<dim> fe_neighbour_face_values (mapping, fe, face_quadrature, neighbour_face_update_flags);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	for (; cell != endc; ++cell)
	{
		fe_values.reinit(cell);
		cell->get_dof_indices (dof_indices);

		if (cell == dof_handler.begin_active())
			assemble_cell_term(fe_values, dof_indices); // not sure what arguments to pass... assemble cell term only once b/c mass and stiffness matrices don't change from cell to cell.

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
void Problem<dim>::assemble_cell_term(const FEValues<dim> &fe_values, const std::vector<types::global_dof_index> &dofs_indices)
{
	const unsigned int n_q_points = fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_values.dofs_per_cell;

	mass_matrix.reinit(dofs_per_cell,dofs_per_cell);
	mass_matrix = 0;
	for (int i = 0; i<dim ; ++i) // initialize each stiffness matrix
	{
		stiffness_matrix[i].reinit(dofs_per_cell,dofs_per_cell);
		stiffness_matrix[i] = 0;
	}



	for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) //this loop creates the mass and dim x stiffness matrices.
	{
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < dofs_per_cell; ++ j)
			{
				mass_matrix += fe_values.shape_value (i, q_index) * fe_values.shape_value (j, q_index) * fe_values.JxW(q_index);
				for (int di = 0; di < dim ; ++di)
				{
					stiffness_matrix[di](i,j) += fe_values.shape_grad_component(i, q_index, di) * fe_values.shape_value(j,q_index) * fe_values.JxW(q_index);
				}
			}
		}
	}
	compute_inverse_mass_matrix(mass_matrix,inverse_mass_matrix);
}

template<int dim>
void Problem<dim>::assemble_face_term(const unsigned int          face_no,
									  const FEFaceValues<dim>     &fe_face_values,
									  const FEFaceValues<dim>     &fe_neighbour_face_values,
									  const std::vector<types::global_dof_index> &dof_indices,
									  const std::vector<types::global_dof_index> &neighbour_dof_indices )
{

	const unsigned int n_q_points = fe_face_values.n_quadrature_points;
	const unsigned int dofs_per_cell = fe_face_values.dofs_per_cell;
	rhs.reinit(dofs_per_cell);
	rhs = 0;
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
			rhs(i) += normal_flux * fe_face_values.shape_values(i,q_index) * fe_face_values.JxW(q_index);
		}
	}
}


/*
 * Not sure how to do the whole indexing for assembling face and cell terms... should I have one big global matrix + vector?...
 * How to propagate the assemble_face_term info from each face to the cell_rhs vector properly without causing intereference?
 */

template<int dim>
void Problem<dim>::perform_runge_kutta_45()
{

	for (int n_iteration = 0; n_iteration < (parameters.final_time - parameters.initial_time)/parameters.delta_t ; ++n_iteration)
	{
		compute_rhs_vector();
		old_solution = current_solution;
		current_solution = old_solution + rhs * parameters.delta_t;
	}
}

template <int dim>
double Problem<dim>::compute_energy(const Vector<double> &u) //need to loop over cells, current function won't work.
{
	Vector<double> product;
	double energy {0};
	mass_matrix.vmult(product, u);

	for (int i = 0; i < u.size(); ++i)
	{
		energy += product[i]*u[i];
	}
	return energy;
}

template class Problem<1>;
