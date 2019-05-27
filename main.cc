/*int main()
{
	Problem problem(parameters);
	problem.setup();
	problem.run();
	problem.output();
	return 0;
}*/

#include "problem.h"
#include "util.h"

int main ()
{
	deallog.depth_console (2);
	Problem burgersEqn;
	burgersEqn.run();
}
