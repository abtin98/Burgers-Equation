#ifndef BURGERS_EQN
#define BURGERS_EQN
#include "util.h"

template <int dim>
  class BurgersEquation
  {
  public:
    static const unsigned int n_components             = dim;
    static const unsigned int first_u_component = 0;
    template <typename InputVector>
    static
    void compute_flux_matrix (const InputVector &W,
                              std::array <std::array
                              <typename InputVector::value_type, dim>,
                              n_components > &flux);
    template <typename InputVector>
    static
    void numerical_normal_flux (const Tensor<1,dim>                &normal,
                                const InputVector                  &Wplus,
                                const InputVector                  &Wminus,
                                const double                        alpha,
                                std::array
                                <typename InputVector::value_type, n_components>
                                &normal_flux);
    enum BoundaryKind
    {
      periodic_boundary,
	  inflow_boundary
    };
    template <typename DataVector>
    static
    void
    compute_Wminus (const BoundaryKind  (&boundary_kind)[n_components],
                    const Tensor<1,dim>  &normal_vector,
                    const DataVector     &Wplus,
                    const Vector<double> &boundary_values,
                    const DataVector     &Wminus);
  };
#endif
