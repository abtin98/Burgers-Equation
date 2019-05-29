#include "util.h"

template <int dim>
  class ConservationLaw
  {
  public:
    ConservationLaw (const char *input_filename);
    void run ();
  private:
    void setup_system ();
    void assemble_system ();
    void assemble_cell_term (const FEValues<dim>             &fe_v,
                             const std::vector<types::global_dof_index> &dofs);
    void assemble_face_term (const unsigned int               face_no,
                             const FEFaceValuesBase<dim>     &fe_v,
                             const FEFaceValuesBase<dim>     &fe_v_neighbor,
                             const std::vector<types::global_dof_index> &dofs,
                             const std::vector<types::global_dof_index> &dofs_neighbor,
                             const bool                       external_face,
                             const unsigned int               boundary_id,
                             const double                     face_diameter);
    std::pair<unsigned int, double> solve (Vector<double> &solution);
    void compute_refinement_indicators (Vector<double> &indicator) const;
    void refine_grid (const Vector<double> &indicator);
    void output_results () const;
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    const FESystem<dim>  fe;
    DoFHandler<dim>      dof_handler;
    const QGauss<dim>    quadrature;
    const QGauss<dim-1>  face_quadrature;
    Vector<double>       old_solution;
    Vector<double>       current_solution;
    Vector<double>       predictor;
    Vector<double>       right_hand_side;
    TrilinosWrappers::SparseMatrix system_matrix;
    Parameters::AllParameters<dim>  parameters;
    ConditionalOStream              verbose_cout;
  };
