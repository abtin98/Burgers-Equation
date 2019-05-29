#ifndef PARAMETERS
#define PARAMETERS


#include "util.h"
#include "burgersEquation.h"
namespace Parameters
  {
    struct Solver
    {
      enum SolverType { gmres, direct };
      SolverType solver;
      enum  OutputType { quiet, verbose };
      OutputType output;
      double linear_residual;
      int max_iterations;
      double ilut_fill;
      double ilut_atol;
      double ilut_rtol;
      double ilut_drop;
      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };

      static const unsigned int max_n_boundaries = 10;

  template <int dim>
    struct AllParameters
	{
	  struct BoundaryConditions
	        {
	          typename BurgersEquation<dim>::BoundaryKind kind[BurgersEquation<dim>::n_components];
	          FunctionParser<dim> values;
	          BoundaryConditions ();
	        };
      AllParameters ();
      double time_step, final_time;
      double theta;
      std::string mesh_filename;
      FunctionParser<dim> initial_conditions;
      BoundaryConditions  boundary_conditions[max_n_boundaries];
      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };

  }

#endif
