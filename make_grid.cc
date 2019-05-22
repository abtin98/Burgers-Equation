#include "util.h"

void make_grid (Triangulation<1> &triangulation, int refine_global)
{
  	GridGenerator::hyper_cube (triangulation);
  	triangulation.refine_global (refine_global);
  	std::cout << "Grid generated!" << std::endl;
}