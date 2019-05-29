#include "util.h"

//Vector<double> RK4 (Vector<double> rhs)
//{
//	for (unsigned int i = 0; i < 4; i++)
//	{
//
//	}
//}
void invert_mass_matrix (const FullMatrix<double> &M, FullMatrix<double> &Minv)
{
	for (unsigned int i = 0; i < M.n_rows(); ++i)
	{
		Minv(i,i) = 1./M(i,i);
	}
}
