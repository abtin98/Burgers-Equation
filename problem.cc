#include "problem.h"

Problem::Problem() :   fe (1), dof_handler (triangulation) //fe(1) indicates polynomial of degree 1.
{
}

void Problem::make_grid()
{
	GridGenerator::hyper_cube (triangulation,0,10);
	triangulation.refine_global (5); //can change number
	std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}

void Problem::setup_system()
{
	std::cout << "system is set up" << std::endl;
}

void Problem::assemble_system()
{
	std::cout << "yeet";
}

void Problem::solve()
{
	std::cout << "system is solved" << std::endl;
}

void Problem::output()
{
	std::cout << "Output saved." << std::endl;
}

void Problem::run()
{
	make_grid();
	setup_system();
	assemble_system();
	solve();
	output();
}
