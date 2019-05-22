#include "Problem.h"

void Problem::make_grid()
{

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
