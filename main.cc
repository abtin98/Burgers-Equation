/*int main()
{
	Problem problem(parameters);
	problem.setup();
	problem.run();
	problem.output();
	return 0;
}*/
/*

#include "problem.h"
#include "util.h"

int main ()
{
	deallog.depth_console (2);
	Problem<1> burgersEqn;
	burgersEqn.run();
}
*/

#include "util.h"
//#include "burgersEquation.h"
//#include "parameters.h"
#include "conservationLaw.h"

int main (int argc, char *argv[])
{
  try
    {
	  const unsigned int dim = 2;
	  Triangulation<dim> triangulation;
	  GridGenerator::hyper_cube(triangulation,0.,10.);
	  triangulation.refine_global(5);
	  std::ofstream out ("mesh.ucd");
	  GridOut grid_out;
	  grid_out.write_ucd(triangulation,out);
	  std::cout << "grid has been written!" << std::endl;


      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
      ConservationLaw<dim> cons (argv[1]);
      cons.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
  return 0;
}
