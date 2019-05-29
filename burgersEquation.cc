#include "util.h"
#include "burgersEquation.h"

template <int dim>
    template <typename InputVector> //removed static
    void BurgersEquation<dim>::compute_flux_matrix (const InputVector &W,
                              std::array <std::array
                              <typename InputVector::value_type, dim>,
                              BurgersEquation<dim>::n_components > &flux)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
    	  flux[d][d] = W[d]*W[d]/2.;
        }
    }
template<int dim>
    template <typename InputVector> //removed static
    void BurgersEquation<dim>::numerical_normal_flux (const Tensor<1,dim>                &normal,
                                const InputVector                  &Wplus,
                                const InputVector                  &Wminus,
                                const double                        alpha,
                                std::array
                                <typename InputVector::value_type, n_components>
                                &normal_flux)
    {
      std::array
      <std::array <typename InputVector::value_type, dim>,
      n_components > iflux, oflux;
      compute_flux_matrix (Wplus, iflux);
      compute_flux_matrix (Wminus, oflux);
      for (unsigned int di=0; di<n_components; ++di)
        {
          normal_flux[di] = 0;
          normal_flux[di] = (Wplus[di]*Wplus[di] + Wminus[di]*Wminus[di]/4.) - (Wplus[di] - Wminus[di])*(Wplus[di] - Wminus[di])/12.;
          //for (unsigned int d=0; d<dim; ++d)
           // normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];
          //normal_flux[di] += 0.5*alpha*(Wplus[di] - Wminus[di]);
        }
    }
template <int dim>
    template <typename DataVector> //removed static
    void
    BurgersEquation<dim>::compute_Wminus (const BoundaryKind  (&boundary_kind)[n_components],
                    const Tensor<1,dim>  &normal_vector,
                    const DataVector     &Wplus,
                    const Vector<double> &boundary_values,
                    const DataVector     &Wminus)
    {
      for (unsigned int c = 0; c < n_components; c++)
        switch (boundary_kind[c])
          {
          case inflow_boundary:
          {
            Wminus[c] = boundary_values(c);
            break;
          }
          case periodic_boundary:
          {
            //Wminus[c] = Wplus[c]; not sure what to do here.
            break;
          }
          default:
            Assert (false, ExcNotImplemented());
          }
    }
