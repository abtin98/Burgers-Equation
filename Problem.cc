#include "Problem.h"

void Problem::make_grid()
{
	GridGenerator::hyper_cube (triangulation);
	triangulation.refine_global (4); //can change number
	std::cout << "Grid generated!" << std::endl;
}

void Problem::assemble_system()
{
	std::cout << "yeet";
}

void Problem::run()
{
	make_grid();
	setup_system();
	assemble_system();
	solve();
	output();
}
