/*int main()
{
	Problem problem(parameters);
	problem.setup();
	problem.run();
	problem.output();
	return 0;
}*/

#include "util.h"
#include "make_grid.h"

int main ()
{
  int refine_global = 4;
  Triangulation<1> triangulation;
  make_grid(triangulation, refine_global);
}
