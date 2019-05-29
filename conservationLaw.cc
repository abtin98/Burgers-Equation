#include "util.h"
#include "conservationLaw.h"
#include "burgersEquation.h"
#include "parameters.h"

  template <int dim>
  ConservationLaw<dim>::ConservationLaw (const char *input_filename)
    :
    mapping (),
    fe (FE_Q<dim>(1), BurgersEquation<dim>::n_components),
    dof_handler (triangulation),
    quadrature (2),
    face_quadrature (2),
    verbose_cout (std::cout, false)
  {
    ParameterHandler prm;
    Parameters::AllParameters<dim>::declare_parameters (prm);
    prm.parse_input (input_filename);
    parameters.parse_parameters (prm);
    verbose_cout.set_condition (parameters.output == Parameters::Solver::verbose);
  }
  template <int dim>
  void ConservationLaw<dim>::setup_system ()
  {
    DynamicSparsityPattern dsp (dof_handler.n_dofs(),
                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    system_matrix.reinit (dsp);
  }
  template <int dim>
  void ConservationLaw<dim>::assemble_system ()
  {
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
    std::vector<types::global_dof_index> dof_indices_neighbor (dofs_per_cell);
    const UpdateFlags update_flags               = update_values
                                                   | update_gradients
                                                   | update_q_points
                                                   | update_JxW_values,
                                                   face_update_flags          = update_values
                                                       | update_q_points
                                                       | update_JxW_values
                                                       | update_normal_vectors,
                                                       neighbor_face_update_flags = update_values;
    FEValues<dim>        fe_v                  (mapping, fe, quadrature,
                                                update_flags);
    FEFaceValues<dim>    fe_v_face             (mapping, fe, face_quadrature,
                                                face_update_flags);
    FESubfaceValues<dim> fe_v_subface          (mapping, fe, face_quadrature,
                                                face_update_flags);
    FEFaceValues<dim>    fe_v_face_neighbor    (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);
    FESubfaceValues<dim> fe_v_subface_neighbor (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_v.reinit (cell);
        cell->get_dof_indices (dof_indices);
        assemble_cell_term(fe_v, dof_indices);
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if (cell->at_boundary(face_no))
            {
              fe_v_face.reinit (cell, face_no);
              assemble_face_term (face_no, fe_v_face,
                                  fe_v_face,
                                  dof_indices,
                                  std::vector<types::global_dof_index>(),
                                  true,
                                  cell->face(face_no)->boundary_id(),
                                  cell->face(face_no)->diameter());
            }
          else
            {
              if (cell->neighbor(face_no)->has_children())
                {
                  const unsigned int neighbor2=
                    cell->neighbor_of_neighbor(face_no);
                  for (unsigned int subface_no=0;
                       subface_no < cell->face(face_no)->n_children();
                       ++subface_no)
                    {
                      const typename DoFHandler<dim>::active_cell_iterator
                      neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);
                      Assert (neighbor_child->face(neighbor2) ==
                              cell->face(face_no)->child(subface_no),
                              ExcInternalError());
                      Assert (neighbor_child->has_children() == false,
                              ExcInternalError());
                      fe_v_subface.reinit (cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);
                      neighbor_child->get_dof_indices (dof_indices_neighbor);
                      assemble_face_term (face_no, fe_v_subface,
                                          fe_v_face_neighbor,
                                          dof_indices,
                                          dof_indices_neighbor,
                                          false,
                                          numbers::invalid_unsigned_int,
                                          neighbor_child->face(neighbor2)->diameter());
                    }
                }
              else if (cell->neighbor(face_no)->level() != cell->level())
                {
                  const typename DoFHandler<dim>::cell_iterator
                  neighbor = cell->neighbor(face_no);
                  Assert(neighbor->level() == cell->level()-1,
                         ExcInternalError());
                  neighbor->get_dof_indices (dof_indices_neighbor);
                  const std::pair<unsigned int, unsigned int>
                  faceno_subfaceno = cell->neighbor_of_coarser_neighbor(face_no);
                  const unsigned int neighbor_face_no    = faceno_subfaceno.first,
                                     neighbor_subface_no = faceno_subfaceno.second;
                  Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                               neighbor_subface_no)
                          == cell,
                          ExcInternalError());
                  fe_v_face.reinit (cell, face_no);
                  fe_v_subface_neighbor.reinit (neighbor,
                                                neighbor_face_no,
                                                neighbor_subface_no);
                  assemble_face_term (face_no, fe_v_face,
                                      fe_v_subface_neighbor,
                                      dof_indices,
                                      dof_indices_neighbor,
                                      false,
                                      numbers::invalid_unsigned_int,
                                      cell->face(face_no)->diameter());
                }
            }
      }
  }
  template <int dim>
  void
  ConservationLaw<dim>::
  assemble_cell_term (const FEValues<dim>             &fe_v,
                      const std::vector<types::global_dof_index> &dof_indices)
  {
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int n_q_points    = fe_v.n_quadrature_points;
    Table<2,Sacado::Fad::DFad<double> >
    W (n_q_points, BurgersEquation<dim>::n_components);
    Table<2,double>
    W_old (n_q_points, BurgersEquation<dim>::n_components);
    Table<3,Sacado::Fad::DFad<double> >
    grad_W (n_q_points, BurgersEquation<dim>::n_components, dim);
    Table<3,double>
    grad_W_old(n_q_points, BurgersEquation<dim>::n_components, dim);
    std::vector<double> residual_derivatives (dofs_per_cell);
    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i].diff (i, dofs_per_cell);
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<BurgersEquation<dim>::n_components; ++c)
        {
          W[q][c]       = 0;
          W_old[q][c]   = 0;
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_W[q][c][d] = 0;
              grad_W_old[q][c][d] = 0;
            }
        }
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
          W[q][c] += independent_local_dof_values[i] *
                     fe_v.shape_value_component(i, q, c);
          W_old[q][c] += old_solution(dof_indices[i]) *
                         fe_v.shape_value_component(i, q, c);
          for (unsigned int d = 0; d < dim; d++)
            {
              grad_W[q][c][d] += independent_local_dof_values[i] *
                                 fe_v.shape_grad_component(i, q, c)[d];
              grad_W_old[q][c][d] += old_solution(dof_indices[i]) *
                                     fe_v.shape_grad_component(i, q, c)[d];
            }
        }
    std::vector <
    std::array <std::array <Sacado::Fad::DFad<double>, dim>, BurgersEquation<dim>::n_components >
    > flux(n_q_points);
    std::vector <
    std::array <std::array <double, dim>, BurgersEquation<dim>::n_components >
    > flux_old(n_q_points);
    std::vector < std::array< Sacado::Fad::DFad<double>, BurgersEquation<dim>::n_components> > forcing(n_q_points);
    std::vector < std::array< double, BurgersEquation<dim>::n_components> > forcing_old(n_q_points);
    for (unsigned int q=0; q<n_q_points; ++q)
      {
        BurgersEquation<dim>::compute_flux_matrix (W_old[q], flux_old[q]);
        BurgersEquation<dim>::compute_forcing_vector (W_old[q], forcing_old[q]);
        BurgersEquation<dim>::compute_flux_matrix (W[q], flux[q]);
        BurgersEquation<dim>::compute_forcing_vector (W[q], forcing[q]);
      }
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      {
        Sacado::Fad::DFad<double> R_i = 0;
        const unsigned int
        component_i = fe_v.get_fe().system_to_component_index(i).first;
        for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
          {
            if (parameters.is_stationary == false)
              R_i += 1.0 / parameters.time_step *
                     (W[point][component_i] - W_old[point][component_i]) *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point);
            for (unsigned int d=0; d<dim; d++)
              R_i -= ( parameters.theta * flux[point][component_i][d] +
                       (1.0-parameters.theta) * flux_old[point][component_i][d] ) *
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.JxW(point);
            for (unsigned int d=0; d<dim; d++)
              R_i += 1.0*std::pow(fe_v.get_cell()->diameter(),
                                  parameters.diffusion_power) *
                     ( parameters.theta * grad_W[point][component_i][d] +
                       (1.0-parameters.theta) * grad_W_old[point][component_i][d] ) *
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.JxW(point);
            R_i -= ( parameters.theta  * forcing[point][component_i] +
                     (1.0 - parameters.theta) * forcing_old[point][component_i] ) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
          }
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          residual_derivatives[k] = R_i.fastAccessDx(k);
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
        right_hand_side(dof_indices[i]) -= R_i.val();
      }
  }
  template <int dim>
  void
  ConservationLaw<dim>::assemble_face_term(const unsigned int           face_no,
                                           const FEFaceValuesBase<dim> &fe_v,
                                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                                           const std::vector<types::global_dof_index>   &dof_indices,
                                           const std::vector<types::global_dof_index>   &dof_indices_neighbor,
                                           const bool                   external_face,
                                           const unsigned int           boundary_id,
                                           const double                 face_diameter)
  {
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
    std::vector<Sacado::Fad::DFad<double> >
    independent_local_dof_values (dofs_per_cell),
                                 independent_neighbor_dof_values (external_face == false ?
                                     dofs_per_cell :
                                     0);
    const unsigned int n_independent_variables = (external_face == false ?
                                                  2 * dofs_per_cell :
                                                  dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; i++)
      {
        independent_local_dof_values[i] = current_solution(dof_indices[i]);
        independent_local_dof_values[i].diff(i, n_independent_variables);
      }
    if (external_face == false)
      for (unsigned int i = 0; i < dofs_per_cell; i++)
        {
          independent_neighbor_dof_values[i]
            = current_solution(dof_indices_neighbor[i]);
          independent_neighbor_dof_values[i]
          .diff(i+dofs_per_cell, n_independent_variables);
        }
    Table<2,Sacado::Fad::DFad<double> >
    Wplus (n_q_points, BurgersEquation<dim>::n_components),
          Wminus (n_q_points, BurgersEquation<dim>::n_components);
    Table<2,double>
    Wplus_old(n_q_points, BurgersEquation<dim>::n_components),
              Wminus_old(n_q_points, BurgersEquation<dim>::n_components);
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
          Wplus[q][component_i] +=  independent_local_dof_values[i] *
                                    fe_v.shape_value_component(i, q, component_i);
          Wplus_old[q][component_i] +=  old_solution(dof_indices[i]) *
                                        fe_v.shape_value_component(i, q, component_i);
        }
    if (external_face == false)
      {
        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int component_i = fe_v_neighbor.get_fe().
                                               system_to_component_index(i).first;
              Wminus[q][component_i] += independent_neighbor_dof_values[i] *
                                        fe_v_neighbor.shape_value_component(i, q, component_i);
              Wminus_old[q][component_i] += old_solution(dof_indices_neighbor[i])*
                                            fe_v_neighbor.shape_value_component(i, q, component_i);
            }
      }
    else
      {
        Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
                ExcIndexRange (boundary_id, 0,
                               Parameters::AllParameters<dim>::max_n_boundaries));
        std::vector<Vector<double> >
        boundary_values(n_q_points, Vector<double>(BurgersEquation<dim>::n_components));
        parameters.boundary_conditions[boundary_id]
        .values.vector_value_list(fe_v.get_quadrature_points(),
                                  boundary_values);
        for (unsigned int q = 0; q < n_q_points; q++)
          {
            BurgersEquation<dim>::compute_Wminus (parameters.boundary_conditions[boundary_id].kind,
                                                 fe_v.normal_vector(q),
                                                 Wplus[q],
                                                 boundary_values[q],
                                                 Wminus[q]);
            BurgersEquation<dim>::compute_Wminus (parameters.boundary_conditions[boundary_id].kind,
                                                 fe_v.normal_vector(q),
                                                 Wplus_old[q],
                                                 boundary_values[q],
                                                 Wminus_old[q]);
          }
      }
    std::vector< std::array < Sacado::Fad::DFad<double>, BurgersEquation<dim>::n_components> >  normal_fluxes(n_q_points);
    std::vector< std::array < double, BurgersEquation<dim>::n_components> >  normal_fluxes_old(n_q_points);
    double alpha;
    switch (parameters.stabilization_kind)
      {
      case Parameters::Flux::constant:
        alpha = parameters.stabilization_value;
        break;
      case Parameters::Flux::mesh_dependent:
        alpha = face_diameter/(2.0*parameters.time_step);
        break;
      default:
        Assert (false, ExcNotImplemented());
        alpha = 1;
      }
    for (unsigned int q=0; q<n_q_points; ++q)
      {
        BurgersEquation<dim>::numerical_normal_flux(fe_v.normal_vector(q),
                                                   Wplus[q], Wminus[q], alpha,
                                                   normal_fluxes[q]);
        BurgersEquation<dim>::numerical_normal_flux(fe_v.normal_vector(q),
                                                   Wplus_old[q], Wminus_old[q], alpha,
                                                   normal_fluxes_old[q]);
      }
    std::vector<double> residual_derivatives (dofs_per_cell);
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
        {
          Sacado::Fad::DFad<double> R_i = 0;
          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const unsigned int
              component_i = fe_v.get_fe().system_to_component_index(i).first;
              R_i += ( parameters.theta * normal_fluxes[point][component_i] +
                       (1.0 - parameters.theta) * normal_fluxes_old[point][component_i] ) *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point);
            }
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            residual_derivatives[k] = R_i.fastAccessDx(k);
          system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
          if (external_face == false)
            {
              for (unsigned int k=0; k<dofs_per_cell; ++k)
                residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell+k);
              system_matrix.add (dof_indices[i], dof_indices_neighbor,
                                 residual_derivatives);
            }
          right_hand_side(dof_indices[i]) -= R_i.val();
        }
  }
  template <int dim>
  std::pair<unsigned int, double>
  ConservationLaw<dim>::solve (Vector<double> &newton_update)
  {
    switch (parameters.solver)
      {
      case Parameters::Solver::direct:
      {
        SolverControl solver_control (1,0);
        TrilinosWrappers::SolverDirect::AdditionalData data (
          parameters.output == Parameters::Solver::verbose);
        TrilinosWrappers::SolverDirect direct (solver_control, data);
        direct.solve (system_matrix, newton_update, right_hand_side);
        return std::pair<unsigned int, double> (solver_control.last_step(),
                                                solver_control.last_value());
      }
      case Parameters::Solver::gmres:
      {
        Epetra_Vector x(View, system_matrix.trilinos_matrix().DomainMap(),
                        newton_update.begin());
        Epetra_Vector b(View, system_matrix.trilinos_matrix().RangeMap(),
                        right_hand_side.begin());
        AztecOO solver;
        solver.SetAztecOption(AZ_output,
                              (parameters.output ==
                               Parameters::Solver::quiet
                               ?
                               AZ_none
                               :
                               AZ_all));
        solver.SetAztecOption(AZ_solver, AZ_gmres);
        solver.SetRHS(&b);
        solver.SetLHS(&x);
        solver.SetAztecOption(AZ_precond,         AZ_dom_decomp);
        solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
        solver.SetAztecOption(AZ_overlap,         0);
        solver.SetAztecOption(AZ_reorder,         0);
        solver.SetAztecParam(AZ_drop,      parameters.ilut_drop);
        solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill);
        solver.SetAztecParam(AZ_athresh,   parameters.ilut_atol);
        solver.SetAztecParam(AZ_rthresh,   parameters.ilut_rtol);
        solver.SetUserMatrix(const_cast<Epetra_CrsMatrix *>
                             (&system_matrix.trilinos_matrix()));
        solver.Iterate(parameters.max_iterations, parameters.linear_residual);
        return std::pair<unsigned int, double> (solver.NumIters(),
                                                solver.TrueResidual());
      }
      }
    Assert (false, ExcNotImplemented());
    return std::pair<unsigned int, double> (0,0);
  }
  template <int dim>
  void ConservationLaw<dim>::output_results () const
  {
    typename BurgersEquation<dim>::Postprocessor
    postprocessor (parameters.schlieren_plot);
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (current_solution,
                              BurgersEquation<dim>::component_names (),
                              DataOut<dim>::type_dof_data,
                              BurgersEquation<dim>::component_interpretation ());
    data_out.add_data_vector (current_solution, postprocessor);
    data_out.build_patches ();
    static unsigned int output_file_number = 0;
    std::string filename = "solution-" +
                           Utilities::int_to_string (output_file_number, 3) +
                           ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
    ++output_file_number;
  }
  template <int dim>
  void ConservationLaw<dim>::run ()
  {
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(parameters.mesh_filename.c_str());
      Assert (input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));
      grid_in.read_ucd(input_file);
    }
    dof_handler.clear();
    dof_handler.distribute_dofs (fe);
    old_solution.reinit (dof_handler.n_dofs());
    current_solution.reinit (dof_handler.n_dofs());
    predictor.reinit (dof_handler.n_dofs());
    right_hand_side.reinit (dof_handler.n_dofs());
    setup_system();
    VectorTools::interpolate(dof_handler,
                             parameters.initial_conditions, old_solution);
    current_solution = old_solution;
    predictor = old_solution;
    if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
        {
          Vector<double> refinement_indicators (triangulation.n_active_cells());
          compute_refinement_indicators(refinement_indicators);
          refine_grid(refinement_indicators);
          setup_system();
          VectorTools::interpolate(dof_handler,
                                   parameters.initial_conditions, old_solution);
          current_solution = old_solution;
          predictor = old_solution;
        }
    output_results ();
    Vector<double> newton_update (dof_handler.n_dofs());
    double time = 0;
    double next_output = time + parameters.output_step;
    predictor = old_solution;
    while (time < parameters.final_time)
      {
        std::cout << "T=" << time << std::endl
                  << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << std::endl;
        std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                  << "   _____________________________________" << std::endl;
        unsigned int nonlin_iter = 0;
        current_solution = predictor;
        while (true)
          {
            system_matrix = 0;
            right_hand_side = 0;
            assemble_system ();
            const double res_norm = right_hand_side.l2_norm();
            if (std::fabs(res_norm) < 1e-10)
              {
                std::printf("   %-16.3e (converged)\n\n", res_norm);
                break;
              }
            else
              {
                newton_update = 0;
                std::pair<unsigned int, double> convergence
                  = solve (newton_update);
                current_solution += newton_update;
                std::printf("   %-16.3e %04d        %-5.2e\n",
                            res_norm, convergence.first, convergence.second);
              }
            ++nonlin_iter;
            AssertThrow (nonlin_iter <= 10,
                         ExcMessage ("No convergence in nonlinear solver"));
          }
        time += parameters.time_step;
        if (parameters.output_step < 0)
          output_results ();
        else if (time >= next_output)
          {
            output_results ();
            next_output += parameters.output_step;
          }
        predictor = current_solution;
        predictor.sadd (2.0, -1.0, old_solution);
        old_solution = current_solution;
      }
  }
}
