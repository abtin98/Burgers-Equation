#include "util.h"

template <int dim>
class Problem
{
	public:
		Problem();
		void run();
	private:
		void make_grid(int n_refinements);
		void setup_system();
		void assemble_system();
		void solve();
		void output();

		Triangulation<dim> triangulation;
		FE_DGQ<dim> fe;
		DoFHandler<dim> dof_handler;

		SparsityPattern sparsity_pattern;
		SparseMatrix<double> system_matrix;

		Vector<double> solution;
		Vector<double> system_rhs;

};
