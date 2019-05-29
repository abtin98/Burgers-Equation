#include "parameters.h"
namespace Parameters
  {
    void Solver::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("linear solver");
      {
        prm.declare_entry("output", "quiet",
                          Patterns::Selection("quiet|verbose"),
                          "State whether output from solver runs should be printed. "
                          "Choices are <quiet|verbose>.");
        prm.declare_entry("method", "gmres",
                          Patterns::Selection("gmres|direct"),
                          "The kind of solver for the linear system. "
                          "Choices are <gmres|direct>.");
        prm.declare_entry("residual", "1e-10",
                          Patterns::Double(),
                          "Linear solver residual");
        prm.declare_entry("max iters", "300",
                          Patterns::Integer(),
                          "Maximum solver iterations");
        prm.declare_entry("ilut fill", "2",
                          Patterns::Double(),
                          "Ilut preconditioner fill");
        prm.declare_entry("ilut absolute tolerance", "1e-9",
                          Patterns::Double(),
                          "Ilut preconditioner tolerance");
        prm.declare_entry("ilut relative tolerance", "1.1",
                          Patterns::Double(),
                          "Ilut relative tolerance");
        prm.declare_entry("ilut drop tolerance", "1e-10",
                          Patterns::Double(),
                          "Ilut drop tolerance");
      }
      prm.leave_subsection();
    }
    void Solver::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("linear solver");
      {
        const std::string op = prm.get("output");
        if (op == "verbose")
          output = verbose;
        if (op == "quiet")
          output = quiet;
        const std::string sv = prm.get("method");
        if (sv == "direct")
          solver = direct;
        else if (sv == "gmres")
          solver = gmres;
        linear_residual = prm.get_double("residual");
        max_iterations  = prm.get_integer("max iters");
        ilut_fill       = prm.get_double("ilut fill");
        ilut_atol       = prm.get_double("ilut absolute tolerance");
        ilut_rtol       = prm.get_double("ilut relative tolerance");
        ilut_drop       = prm.get_double("ilut drop tolerance");
      }
      prm.leave_subsection();
    }

      static const unsigned int max_n_boundaries = 10;

  template <int dim>
    AllParameters<dim>::BoundaryConditions::BoundaryConditions ()
      :
      values (BurgersEquation<dim>::n_components)
    {
      for (unsigned int c=0; c<BurgersEquation<dim>::n_components; ++c)
        kind[c] = BurgersEquation<dim>::periodic_boundary; //not sure what to do here.
    }
    template <int dim>
    AllParameters<dim>::AllParameters ()
      :
      time_step(.1),
      final_time(3.),
      theta(.5),
      initial_conditions (BurgersEquation<dim>::n_components)
    {}
    template <int dim>
    void
    AllParameters<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry("mesh", "grid.inp",
                        Patterns::Anything(),
                        "intput file name");
      prm.enter_subsection("time stepping");
      {
        prm.declare_entry("time step", "0.1",
                          Patterns::Double(0),
                          "simulation time step");
        prm.declare_entry("final time", "10.0",
                          Patterns::Double(0),
                          "simulation end time");
        prm.declare_entry("theta scheme value", "1.0",
                          Patterns::Double(0,1),
                          "value for theta that interpolated between explicit "
                          "Euler (theta=0), Crank-Nicolson (theta=0.5), and "
                          "implicit Euler (theta=1).");
      }
      prm.leave_subsection();
      for (unsigned int b=0; b<max_n_boundaries; ++b)
        {
          prm.enter_subsection("boundary_" +
                               Utilities::int_to_string(b));
          {
            for (unsigned int di=0; di<BurgersEquation<dim>::n_components; ++di)
              {
                prm.declare_entry("w_" + Utilities::int_to_string(di),
                                  "periodic",
                                  Patterns::Selection("inflow|periodic"),
                                  "<inflow|periodic>");
                prm.declare_entry("w_" + Utilities::int_to_string(di) +
                                  " value", "0.0",
                                  Patterns::Anything(),
                                  "expression in x,y,z");
              }
          }
          prm.leave_subsection();
        }
      prm.enter_subsection("initial condition");
      {
        for (unsigned int di=0; di<BurgersEquation<dim>::n_components; ++di)
          prm.declare_entry("w_" + Utilities::int_to_string(di) + " value",
                            "0.0",
                            Patterns::Anything(),
                            "expression in x,y,z");
      }
      prm.leave_subsection();
      Parameters::Solver::declare_parameters (prm);
    }
    template <int dim>
    void
    AllParameters<dim>::parse_parameters (ParameterHandler &prm)
    {
      mesh_filename = prm.get("mesh");
      prm.enter_subsection("time stepping");
      {
        time_step = prm.get_double("time step");
        if (time_step == 0)
          {
            time_step = 1.0;
            final_time = 1.0;
          }
        final_time = prm.get_double("final time");
        theta = prm.get_double("theta scheme value");
      }
      prm.leave_subsection();
      for (unsigned int boundary_id=0; boundary_id<max_n_boundaries;
           ++boundary_id)
        {
          prm.enter_subsection("boundary_" +
                               Utilities::int_to_string(boundary_id));
          {
            std::vector<std::string>
            expressions(BurgersEquation<dim>::n_components, "0.0");
            for (unsigned int di=0; di<BurgersEquation<dim>::n_components; ++di)
              {
                const std::string boundary_type
                  = prm.get("w_" + Utilities::int_to_string(di));
                if (di < dim)
                  boundary_conditions[boundary_id].kind[di]
                    = BurgersEquation<dim>::periodic_boundary;
                else if (boundary_type == "inflow")
                  boundary_conditions[boundary_id].kind[di]
                    = BurgersEquation<dim>::inflow_boundary;
                else
                  AssertThrow (false, ExcNotImplemented());
                expressions[di] = prm.get("w_" + Utilities::int_to_string(di) +
                                          " value");
              }
            boundary_conditions[boundary_id].values
            .initialize (FunctionParser<dim>::default_variable_names(),
                         expressions,
                         std::map<std::string, double>());
          }
          prm.leave_subsection();
        }
      prm.enter_subsection("initial condition");
      {
        std::vector<std::string> expressions (BurgersEquation<dim>::n_components,
                                              "0.0");
        for (unsigned int di = 0; di < BurgersEquation<dim>::n_components; di++)
          expressions[di] = prm.get("w_" + Utilities::int_to_string(di) +
                                    " value");
        initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                       expressions,
                                       std::map<std::string, double>());
      }
      prm.leave_subsection();
      Parameters::Solver::parse_parameters (prm);
    }
  }

